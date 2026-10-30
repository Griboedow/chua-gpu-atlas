// Exercise the exported viewer controller without a browser or WebGL renderer.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');

async function main() {
  const folder = process.argv[2];
  if (!folder) throw new Error('Usage: node research/check_chua_slope_volume_viewer.cjs OUTPUT_DIRECTORY');
  const html = fs.readFileSync(path.join(folder, 'volumes.html'), 'utf8');
  const start = html.lastIndexOf('const model=');
  const end = html.indexOf('\nrefresh();', start) + '\nrefresh();'.length;
  assert(start >= 0 && end > start);
  const graph = {data: Array.from({length: 12}, () => ({})), layout: {}};
  const elements = {volumes: graph, scope: {value: 'interior'}, mode: {value: 'regular'},
    'show-plane': {checked: false}, 'alpha-slice': {value: '20'}, 'beta-slice': {value: '20'},
    'alpha-value': {}, 'beta-value': {}, status: {dataset: {}}, reset: {}};
  for (const element of Object.values(elements)) {
    element.events = {};
    element.addEventListener = (event, fn) => { element.events[event] = fn; };
  }
  const context = vm.createContext({document: {getElementById: id => elements[id]}, Plotly: {
    restyle: async (g, update, indices = g.data.map((_, i) => i)) => {
      for (const [key, values] of Object.entries(update)) {
        indices.forEach((index, i) => { g.data[index][key] = values[i % values.length]; });
      }
    },
    relayout: async (g, update) => { Object.assign(g.layout, update); }
  }});
  await vm.runInContext(html.slice(start, end), context);
  const {model, traces, refresh} = vm.runInContext('({model, traces, refresh})', context);
  function verify() {
    const scope = elements.scope.value, mode = elements.mode.value;
    for (const kind of ['alpha', 'beta']) {
      const t = traces[kind], d = model[kind], level = Number(elements[kind + '-slice'].value);
      assert.equal(graph.data[t.map].z, d.labels[level]);
      assert.equal(graph.data[t.map].visible, true);
      for (const s of ['full', 'interior']) {
        assert.equal(graph.data[t[s + '2']].visible, s === scope && mode !== 'chaos');
        assert.equal(graph.data[t[s + '4']].visible, s === scope && mode !== 'regular');
      }
      const plane = graph.data[t.plane];
      assert.equal(plane.visible, elements['show-plane'].checked);
      const xi = d.m0.map((v, i) => i).filter(i => scope === 'full' || d.m0[i] <= -1.06);
      const yi = d.m1.map((v, i) => i).filter(i => scope === 'full' || (d.m1[i] >= -.90 && d.m1[i] <= -.20));
      assert.equal(plane.x.length, xi.length); assert.equal(plane.y.length, yi.length);
      for (const y of [0, yi.length - 1]) for (const x of [0, xi.length - 1]) {
        assert.equal(plane.surfacecolor[y][x], d.labels[level][yi[y]][xi[x]]);
        assert.equal(plane.z[y][x], d.control[level]);
      }
      assert.equal(elements[kind + '-value'].textContent, d.control[level].toFixed(5));
    }
  }
  verify(); let combinations = 0;
  for (const scope of ['full', 'interior']) for (const mode of ['regular', 'chaos', 'both']) {
    for (const show of [false, true]) for (const level of [0, 20, 40]) {
      elements.scope.value = scope; elements.mode.value = mode; elements['show-plane'].checked = show;
      elements['alpha-slice'].value = String(level); elements['beta-slice'].value = String(40 - level);
      await refresh(); verify(); combinations++;
    }
  }
  elements['alpha-slice'].value = '0'; const first = refresh();
  elements['alpha-slice'].value = '37'; elements['beta-slice'].value = '3';
  await Promise.all([first, refresh()]); verify();
  elements.reset.events.click();
  assert.equal(graph.layout['scene.camera'].up.z, 1);
  assert.equal(graph.layout['scene2.camera'].up.z, 1);
  const result = {status: 'passed',control_combinations: combinations,latest_slider_update: 'passed',camera_reset: 'passed',
    browser_render_check: 'Not performed: local browser use is disabled in Work settings.',
    scope: 'Exported JavaScript executed with DOM and Plotly state stubs. Checks controller state and data selection, not actual WebGL rendering.'};
  fs.writeFileSync(path.join(folder, 'viewer-checks.json'), JSON.stringify(result, null, 2) + '\n');
  console.log(JSON.stringify(result));
}
main().catch(error => { console.error(error); process.exitCode = 1; });
