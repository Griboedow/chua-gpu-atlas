"""Additional exact state-sum controls and outward-rounded certificate checks.

This audit does not reconstruct or certify orbits. It tests combinatorial edge
cases and rechecks pair separation using the serialized enclosing Arb balls.
"""
from search_chua_links import ROOT, OUT, read, save, digest
from flint import arb, ctx
from spherogram import Link
from chua_link_bracket import frontier_jones
from chua_link_diagrams import jones_A, unlink_jones
from chua_link_invariants import split_union_jones


def controls():
    cases = []
    for count in (1, 2, 3):
        link = Link([])
        link.unlinked_unknot_components = count
        cases.append((f'free_unlink_{count}', link, unlink_jones(count)))
    trefoil = Link('3_1')
    for count in (1, 2):
        link = Link(trefoil.PD_code())
        link.unlinked_unknot_components = count
        expected = split_union_jones([jones_A(trefoil)] + [{'0': 1}] * count)
        cases.append((f'trefoil_with_{count}_free_components', link, expected))
    first, second = trefoil.PD_code(), trefoil.mirror().PD_code()
    offset = max(v for crossing in first for v in crossing) + 1
    split = Link(first + [[v + offset for v in crossing] for crossing in second])
    cases.append(('split_trefoil_and_mirror', split,
                  split_union_jones([jones_A(trefoil), jones_A(trefoil.mirror())])))
    # All four positions of the two repeated labels test immediate loop closure.
    ri_codes = [[0, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 0], [1, 0, 0, 1]]
    for i, pd in enumerate(ri_codes):
        for free in (0, 1):
            link = Link([pd])
            link.unlinked_unknot_components = free
            cases.append((f'Reidemeister_I_{i}_free_{free}', link, unlink_jones(free + 1)))
    rows = []
    for name, link, expected in cases:
        grouped, diagnostic = frontier_jones(link)
        exhaustive = jones_A(link)
        assert grouped == exhaustive == expected, (name, grouped, exhaustive, expected)
        rows.append(dict(name=name, pd_code=link.PD_code(),
                         free_components=link.unlinked_unknot_components,
                         frontier_polynomial=grouped, exhaustive_polynomial=exhaustive,
                         expected_polynomial=expected, diagnostic=diagnostic, passed=True))
    assert len(rows) == 14
    return rows


def interval_separation():
    ctx.dps = 180
    rows = []
    for digits in (140, 160):
        data = read(OUT / f'certificates-d{digits}.json')
        assert data['passed']
        assert data['shifted_root_control']['rejected']
        assert data['shifted_root_control']['reason'].startswith('KrawczykRejection:')
        for witness in data['witnesses']:
            records = witness['records']
            assert len(records) == 2 and all(r['primitive_period'] for r in records)
            start = [arb(v['ball']) for v in records[0]['root']]
            margins = []
            for point in records[1]['upward_section_points']:
                differences = [start[i] - arb(point[i]['ball']) for i in range(2)]
                separated = [v for v in differences if not v.contains(0)]
                assert separated, (digits, witness['name'], point)
                margin = max(abs(v).lower() for v in separated)
                assert margin > 0
                margins.append(str(margin))
            rows.append(dict(digits=digits, pair=witness['name'],
                             second_orbit_upward_crossings=len(margins),
                             separating_coordinate_gap_lower_bounds=margins,
                             enclosure_source='Outward-rounded serialized ball strings',
                             passed=True))
    return rows


def main():
    polynomials = controls()
    separation = interval_separation()
    paths = [ROOT/'research'/name for name in [
        'audit_chua_link_algorithms.py', 'chua_link_bracket.py',
        'chua_link_diagrams.py', 'chua_link_invariants.py',
        'certify_chua_link_witnesses.py']]
    paths += [OUT/f'certificates-d{d}.json' for d in (140, 160)]
    save(OUT/'independent-audit.json', dict(
        passed=True, polynomial_control_count=len(polynomials),
        polynomial_controls=polynomials, interval_pair_checks=separation,
        source_sha256={p.relative_to(ROOT).as_posix(): digest(p) for p in paths},
        hash_normalization='CRLF converted to LF for Python and JSON files.',
        scope='Exact state-sum edge cases and pair separation from saved outward-rounded balls. '
              'No new orbit existence calculation and no ambient-isotopy certification.'))
    print(f'{len(polynomials)} exact polynomial controls and {len(separation)} interval pair checks passed')


if __name__ == '__main__':
    main()
