// Optional development-only browser smoke test. Not required to run the Python project.
const path = require('node:path');
const fs = require('node:fs');
const {pathToFileURL} = require('node:url');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'playwright');

(async () => {
  const output = path.resolve(process.argv[2] || 'outputs');
  const browser = await chromium.launch({headless: true, channel: 'msedge',
    args: ['--enable-unsafe-swiftshader']});
  try {
    const page = await browser.newPage({viewport: {width: 1500, height: 1050}, deviceScaleFactor: 1});
    const errors = [];
    page.on('pageerror', e => errors.push(String(e)));
    await page.goto(pathToFileURL(path.join(output, 'index.html')).href, {waitUntil: 'load', timeout: 60000});
    await page.waitForFunction(() => document.querySelectorAll('.js-plotly-plot').length === 6 &&
      [...document.querySelectorAll('.js-plotly-plot')].every(el => el._fullLayout), {timeout: 60000});
    await page.evaluate(async () => { await Promise.all([...document.querySelectorAll('.js-plotly-plot')].map(p=>Plotly.Plots.resize(p))); });
    const dimensions = await page.locator('.js-plotly-plot').evaluateAll(els => els.map(el => ({client:el.clientWidth,plot:el._fullLayout.width,
      ...(el._fullLayout.scene ? {domain:el._fullLayout.scene.domain,size:el._fullLayout._size,
        ranges:['xaxis','yaxis','zaxis'].map(k=>el._fullLayout.scene[k].range)}:{})})));
    const plots = await page.locator('.js-plotly-plot').count();
    const canvases = await page.locator('canvas').count();
    await page.locator('section').first().scrollIntoViewIfNeeded();
    await page.screenshot({path: path.join(output, 'browser-preview.png')});
    const region = page.locator('.js-plotly-plot').nth(1);
    await region.getByText('LLE points', {exact: true}).click();
    await page.waitForFunction(() => document.querySelectorAll('.js-plotly-plot')[1].data[2].visible === true);
    const switched = await region.evaluate(el => el.data.map(t => t.visible));
    await region.getByText('Any seed', {exact: true}).click();
    await region.screenshot({path: path.join(output, 'region-preview.png')});
    const cameras = await page.locator('.js-plotly-plot').evaluateAll(els => els.filter(el=>el._fullLayout.scene).map(el=>el._fullLayout.scene._scene.getCamera()));
    const overflow = await page.evaluate(() => document.documentElement.scrollWidth > innerWidth);
    const result = {plots, canvases, dimensions, cameras, switched, horizontalOverflow: overflow, pageErrors: errors};
    fs.writeFileSync(path.join(output, 'browser-check.json'), JSON.stringify(result, null, 2));
    console.log(JSON.stringify(result));
    if (plots !== 6 || canvases < 4 || errors.length || overflow) process.exitCode = 1;
  } finally { await browser.close(); }
})().catch(error => {console.error(error); process.exitCode = 1;});
