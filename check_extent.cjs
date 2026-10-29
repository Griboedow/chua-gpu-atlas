const path=require('node:path');
const fs=require('node:fs');
const {pathToFileURL}=require('node:url');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
(async()=>{
 const output=path.resolve('outputs/extent');
 const browser=await chromium.launch({headless:true,channel:'msedge',args:['--enable-unsafe-swiftshader']});
 try {
  const page=await browser.newPage({viewport:{width:1500,height:1050},deviceScaleFactor:1});
  const errors=[];page.on('pageerror',e=>errors.push(String(e)));
  await page.goto(pathToFileURL(path.join(output,'index.html')).href,{waitUntil:'load',timeout:60000});
  await page.waitForFunction(()=>document.querySelectorAll('.js-plotly-plot').length===3 &&
    [...document.querySelectorAll('.js-plotly-plot')].every(p=>p._fullLayout));
  await page.evaluate(async()=>Promise.all([...document.querySelectorAll('.js-plotly-plot')].map(p=>Plotly.Plots.resize(p))));
  await page.locator('#extent-linear').screenshot({path:path.join(output,'linear-preview.png')});
  await page.locator('#extent-log').screenshot({path:path.join(output,'log-preview.png')});
  await page.locator('#extent-log .legend .traces').filter({hasText:'Scouted points'}).locator('.legendtoggle').click();
  await page.waitForFunction(()=>document.getElementById('extent-log').data[0].visible===true);
  const result={plots:await page.locator('.js-plotly-plot').count(),pageErrors:errors,
    horizontalOverflow:await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth)};
  fs.writeFileSync(path.join(output,'browser-check.json'),JSON.stringify(result,null,2));
  console.log(JSON.stringify(result));
  if(errors.length||result.horizontalOverflow)process.exitCode=1;
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
