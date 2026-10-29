// Optional browser smoke test for the detailed scientific report.
const path=require('node:path');
const fs=require('node:fs');
const {pathToFileURL}=require('node:url');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
(async()=>{
 const output=path.resolve('outputs/detailed');
 const browser=await chromium.launch({headless:true,channel:'msedge',args:['--enable-unsafe-swiftshader']});
 try {
  const page=await browser.newPage({viewport:{width:1500,height:1050},deviceScaleFactor:1});
  const errors=[];page.on('pageerror',e=>errors.push(String(e)));
  await page.goto(pathToFileURL(path.join(output,'index.html')).href,{waitUntil:'load',timeout:60000});
  await page.waitForFunction(()=>document.querySelectorAll('.js-plotly-plot').length===10 &&
    [...document.querySelectorAll('.js-plotly-plot')].every(p=>p._fullLayout),{timeout:60000});
  await page.evaluate(async()=>Promise.all([...document.querySelectorAll('.js-plotly-plot')].map(p=>Plotly.Plots.resize(p))));
  await page.locator('#region-classic').screenshot({path:path.join(output,'region-preview.png')});
  await page.locator('#window-zoom').screenshot({path:path.join(output,'zoom-preview.png')});
  const region=page.locator('#region-classic');
  await region.locator('.updatemenu-item-text').filter({hasText:/^Any seed$/}).click();
  await region.locator('.updatemenu-item-text').filter({hasText:/^Both outer regions$/}).click();
  await page.waitForFunction(()=>document.getElementById('region-classic').data[3].visible===true);
  const slice=page.locator('#window-slices');
  await slice.scrollIntoViewIfNeeded();
  // Use Plotly's public restyle API with the actual slider payload.
  const frame=await slice.evaluate(async el=>{
    const step=el.layout.sliders[0].steps[0];
    await Plotly.restyle(el,...step.args);
    return {sliceCount:el.layout.sliders[0].steps.length,rows:el.data[0].z.length};
  });
  await page.locator('#basin').screenshot({path:path.join(output,'basin-preview.png')});
  const overflow=await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth);
  const result={plotCount:await page.locator('.js-plotly-plot').count(),frame,overflow,pageErrors:errors};
  fs.writeFileSync(path.join(output,'browser-check.json'),JSON.stringify(result,null,2));
  console.log(JSON.stringify(result));
  if(errors.length||overflow||frame.sliceCount!==41||frame.rows!==41)process.exitCode=1;
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
