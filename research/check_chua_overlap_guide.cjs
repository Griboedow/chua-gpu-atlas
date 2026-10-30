// Check the generated offline guide's real script and its embedded sections.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const root = path.resolve(__dirname, '..');
const html = fs.readFileSync(path.join(root, 'results/chua-overlap-geometry/guide.html'), 'utf8');
assert(!/\b(?:STORY|ORBITS|SECTIONS|PREDICTION|VOLUMES|DATA)\b/.test(html));
assert(!/<(?:script|link)[^>]+(?:src|href)=["']https?:/i.test(html));
assert.equal((html.match(/<img[^>]+src="data:image\/png;base64,/g) || []).length, 5);
const script = html.match(/<script>([\s\S]*?)<\/script>/)[1];
let change;
const elements = {
  alpha: {value: '8.86', addEventListener(event, callback) {assert.equal(event, 'change'); change = callback;}},
  description: {textContent: ''},
  view: {innerHTML: ''},
};
vm.runInNewContext(script, {document: {getElementById(id) {assert(id in elements); return elements[id];}}}, {timeout: 3000});
assert(elements.description.textContent.includes('at least three'));
for (const alpha of ['8.84', '8.85', '8.86', '8.865']) {
  elements.alpha.value = alpha; change();
  assert(elements.view.innerHTML.includes('<svg'));
  assert(elements.view.innerHTML.includes('</svg>'));
  assert(!elements.view.innerHTML.includes('NaN'));
  assert(elements.description.textContent.length > 70);
}
console.log('Offline guide passed: 5 embedded images, 4 selectable SVG sections, real change handler, no remote runtime dependency.');
