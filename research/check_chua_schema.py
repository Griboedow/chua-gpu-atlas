"""Verify renamed records and arrays against the preceding Git publication."""
from pathlib import Path
import argparse,ast,io,json,re,subprocess,sys
from urllib.parse import urlsplit,unquote
import numpy as np

root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'research'))
from publication_schema import current_names
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=root/'.publication/schema-check.json')
args=parser.parse_args()
for commit,name in [('0fc47ea','results/chua-layer-prediction/predictions-frozen.json'),
                    ('112b01f','results/chua-followup/sheets/forecasts-frozen.json')]:
    original=subprocess.check_output(['git','show',commit+':'+name],cwd=root).decode('utf-8')
    assert json.loads(current_names(original))==json.loads((root/name).read_text(encoding='utf-8')),(commit,name)
base='37dd180'
names=subprocess.check_output(['git','ls-tree','-r','--name-only',base],cwd=root,text=True).splitlines()
selected=[n for n in names if Path(n).suffix in {'.json','.npz','.py','.md'}]
data=subprocess.check_output(['git','cat-file','--batch'],cwd=root,input=''.join(f'{base}:{n}\n' for n in selected).encode())
position=0;counts={'json':0,'arrays':0,'python':0,'math_blocks':0};paths=[]
for name in selected:
    end=data.index(b'\n',position);size=int(data[position:end].split()[-1]);old=data[end+1:end+1+size];position=end+size+2
    path=root/current_names(name);assert path.is_file(),str(path)
    if str(path.relative_to(root)).replace('\\','/')!=name:paths.append([name,path.relative_to(root).as_posix()])
    if path.suffix=='.json':
        if name=='results/chua-layer-atlas/manifest.json':continue
        expected=json.loads(current_names(old.decode('utf-8')));actual=json.loads(path.read_text(encoding='utf-8'))
        assert expected==actual,('record changed beyond names',name)
        counts['json']+=1
    elif path.suffix=='.npz':
        with np.load(io.BytesIO(old),allow_pickle=False) as a,np.load(path,allow_pickle=False) as b:
            assert {current_names(k) for k in a.files}==set(b.files),name
            for k in a.files:
                x=a[k];y=b[current_names(k)]
                assert x.dtype==y.dtype and np.array_equal(x,y,equal_nan=x.dtype.kind in 'fc'),(name,k)
                counts['arrays']+=1
    elif path.suffix=='.py':ast.parse(path.read_text(encoding='utf-8'));counts['python']+=1
    elif path.suffix=='.md':
        a=re.findall(r'```math\n[\s\S]*?```',old.decode().replace('\r\n','\n'))
        b=re.findall(r'```math\n[\s\S]*?```',path.read_text(encoding='utf-8'))
        assert a==b,('math changed',name);counts['math_blocks']+=len(a)

def anchors(p):
    s=p.read_text(encoding='utf-8');seen={};out=set(re.findall(r'<a id="([^"]+)"',s))
    for line in s.splitlines():
        if not re.match(r'^#{1,6} ',line):continue
        h=re.sub(r'[^\w\-\s]','',re.sub(r'^#+\s+','',line).lower()).replace(' ','-')
        n=seen.get(h,0);seen[h]=n+1;out.add(h+(f'-{n}' if n else ''))
    return out
links=0
for name in names:
    if not name.endswith('.md'):continue
    p=root/current_names(name);s=p.read_text(encoding='utf-8')
    for url in re.findall(r'\]\(([^\n)]+)\)',s):
        u=urlsplit(url)
        if u.scheme or url.startswith('//'):continue
        target=(p.parent/unquote(u.path)).resolve() if u.path else p
        assert target.exists(),(name,url)
        if u.fragment and target.suffix=='.md':assert unquote(u.fragment) in anchors(target),(name,url)
        links+=1
counts['local_links']=links
print(json.dumps(counts));print('All records and arrays match the prior publication after identifier renaming.')
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps({'base_commit':base,'counts':counts,'renamed_paths':paths,'passed':True},indent=2)+'\n',encoding='utf-8')
