"""Exact sufficient tests for non-split links represented by oriented diagrams.

All meridians map to the same variable t.  A split link with at least two
components has zero single-variable Alexander polynomial, hence every maximal
reduced Fox minor vanishes.  A nonzero evaluated minor proves non-splitness;
zero is inconclusive.  The tests concern the supplied diagram, not its numerical
extraction from an ODE trajectory.  Individual knottedness is not joint linking.
"""
from __future__ import annotations

from bisect import bisect_right
from collections import Counter
import hashlib
import json
from pathlib import Path

import sympy as sp


def component_count(link):
    return len(link.link_components) + int(link.unlinked_unknot_components)


def fox_matrix(link, t):
    """Return the oriented Wirtinger matrix, or None for a free split component.

    The outgoing under-arc labels each row.  Arc numbering wraps independently
    around each component.  An over-only component is an unknotted split
    component: it can be lifted above the rest of the diagram.
    """
    if link.unlinked_unknot_components:
        return None
    sequences = link.link_components
    under_positions = [[i for i, entry in enumerate(seq)
                        if entry.is_under_crossing()] for seq in sequences]
    if any(not positions for positions in under_positions):
        return None
    n = len(link.crossings)
    if not n:
        return sp.zeros(0, 0)
    assert sum(map(len, under_positions)) == n
    matrix = sp.zeros(n, n)
    over_arcs = {}
    offset = 0
    for seq, positions in zip(sequences, under_positions):
        for i, entry in enumerate(seq):
            if entry.is_over_crossing():
                arc = (bisect_right(positions, i) - 1) % len(positions)
                over_arcs[id(entry.crossing)] = offset + arc
        offset += len(positions)
    assert len(over_arcs) == n
    offset = 0
    for seq, positions in zip(sequences, under_positions):
        size = len(positions)
        for local, pos in enumerate(positions):
            crossing = seq[pos].crossing
            row = offset + local
            incoming = offset + (local - 1) % size
            matrix[row, over_arcs[id(crossing)]] += 1 - t
            if crossing.sign > 0:
                matrix[row, incoming] += t
                matrix[row, row] -= 1
            else:
                matrix[row, incoming] -= 1
                matrix[row, row] += t
        offset += size
    assert all(sum(matrix.row(i)) == 0 for i in range(n))
    return matrix


def alexander_minor(link, t=-1):
    """Exact integer cofactor evaluation, not a normalized Alexander polynomial.

    Delete the final row and column.  Nonzero implies a non-split link when
    there are at least two components.  At t=2 units can change the numerical
    value between equivalent diagrams, but not whether it vanishes.
    """
    if not isinstance(t, (int, sp.Integer)) or t == 0:
        raise ValueError('t must be a nonzero integer')
    count = component_count(link)
    if count < 1:
        raise ValueError('The empty link has no Alexander minor in this API')
    if not link.crossings:
        return 1 if count == 1 else 0
    matrix = fox_matrix(link, sp.Integer(t))
    if matrix is None:
        return 0
    return int(matrix[:-1, :-1].det(method='domain-ge'))


def fox_nonsplit_test(link, evaluations=(-1, 2)):
    values = {str(int(t)): alexander_minor(link, t) for t in evaluations}
    count = component_count(link)
    n = len(link.crossings)
    return dict(components=count, crossings=n, evaluations=values,
                omitted_row=n - 1 if n else None,
                omitted_column=n - 1 if n else None,
                non_split_by_alexander=count > 1 and any(values.values()),
                scope='Exact arithmetic on the supplied oriented diagram; zero is inconclusive')


def multiply_polynomials(left, right):
    """Multiply Laurent polynomials encoded by integer-exponent string keys."""
    out = Counter()
    for a, x in left.items():
        for b, y in right.items():
            out[int(a) + int(b)] += int(x) * int(y)
    return {str(e): int(v) for e, v in sorted(out.items()) if v}


def split_union_jones(component_polynomials):
    """Jones polynomial in A for the fully split union of the given knots.

    V(K1 disjoint-union ... disjoint-union Kr) = delta**(r-1) product V(Ki),
    delta = -A**2 - A**-2.  Inequality excludes complete splitting; for r>2
    it need not exclude a split partition into nontrivial sublinks.
    """
    polynomials = list(component_polynomials)
    if not polynomials:
        raise ValueError('At least one component polynomial is required')
    out = {'0': 1}
    for polynomial in polynomials:
        out = multiply_polynomials(out, polynomial)
    for _ in range(len(polynomials) - 1):
        out = multiply_polynomials(out, {'2': -1, '-2': -1})
    return out


def diagrammatic_split_count(link):
    """Number of visibly split pieces; one does not prove non-splitness.

    Includes isolated unknotted components removed by diagram simplification.
    The count can increase after Reidemeister simplification of a split link.
    """
    pieces = len(link.split_link_diagram()) if link.crossings else 0
    return pieces + int(link.unlinked_unknot_components)


def controls(output=None):
    """Replay exact standard-link and split-knotted controls, without Sage."""
    from spherogram import Link
    from chua_knot_diagrams import fox_coefficients
    from chua_link_diagrams import jones_A, unlink_jones

    rows = []
    for name, table, expected in [('hopf', 'L2a1', True),
                                  ('whitehead', 'L5a1', True),
                                  ('borromean', 'L6a4', True),
                                  ('trefoil', '3_1', False)]:
        link = Link(table)
        row = dict(name=name, table_name=table, pd_code=link.PD_code(),
                   free_components=link.unlinked_unknot_components,
                   **fox_nonsplit_test(link))
        assert row['non_split_by_alexander'] == expected, row
        if table == '3_1':
            coefficients = fox_coefficients(link)
            t = sp.Symbol('t')
            matrix = fox_matrix(link, t)
            polynomial = sp.Poly(matrix[:-1, :-1].det(), t)
            lo = min(e[0] for e, c in polynomial.terms())
            values = [int(polynomial.nth(i)) for i in range(lo, polynomial.degree() + 1)]
            if values[-1] < 0:
                values = [-v for v in values]
            assert values == coefficients
            row['single_component_comparison'] = coefficients
        row['passed'] = True
        rows.append(row)

    trefoil = Link('3_1')
    first = trefoil.PD_code()
    second = trefoil.mirror().PD_code()
    offset = max(v for crossing in first for v in crossing) + 1
    split_pair = Link(first + [[v + offset for v in crossing] for crossing in second])
    split_polynomial = split_union_jones([jones_A(trefoil), jones_A(trefoil.mirror())])
    assert jones_A(split_pair) == split_polynomial
    assert split_polynomial != unlink_jones(2)
    assert diagrammatic_split_count(split_pair) == 2
    row = dict(name='split_trefoil_and_mirror', pd_code=split_pair.PD_code(),
               free_components=0, **fox_nonsplit_test(split_pair),
               jones_A=jones_A(split_pair), split_union_jones_A=split_polynomial,
               diagrammatic_split_pieces=diagrammatic_split_count(split_pair))
    assert not row['non_split_by_alexander'] and not any(row['evaluations'].values())
    row['passed'] = True
    rows.append(row)
    for size in (1, 2, 3):
        link = Link([])
        link.unlinked_unknot_components = size
        row = dict(name=f'unlink_{size}', pd_code=[], free_components=size,
                   **fox_nonsplit_test(link))
        assert not row['non_split_by_alexander']
        assert set(row['evaluations'].values()) == ({1} if size == 1 else {0})
        assert diagrammatic_split_count(link) == size
        row['passed'] = True
        rows.append(row)
    result = dict(passed=all(row['passed'] for row in rows), rows=rows,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
                  scope='Exact combinatorial controls. Nonzero Alexander minor is sufficient, never necessary, for non-splitness.')
    if output is not None:
        path = Path(output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    return result


if __name__ == '__main__':
    target = Path(__file__).resolve().parents[1] / 'results/chua-link-search/invariant-controls.json'
    result = controls(target)
    print(json.dumps({'passed': result['passed'], 'controls': len(result['rows']),
                      'output': str(target)}))
