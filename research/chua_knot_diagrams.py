"""Polygonal knot diagrams and exact Alexander polynomials.

Crossings use floating-point geometry. Polynomial arithmetic is exact once the
diagram is fixed; this does not certify the topology of the underlying ODE orbit.
"""
from __future__ import annotations

import numpy as np
from numba import njit
from spherogram import Link
import sympy as sp


@njit(cache=True)
def intersections(points):
    """Crossings of a closed polygon under projection to its first two axes."""
    n = len(points)
    rows = []
    for i in range(n):
        a = points[i]; b = points[(i + 1) % n]; u = b - a
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            c = points[j]; d = points[(j + 1) % n]; v = d - c
            if max(a[0], b[0]) < min(c[0], d[0]) or max(c[0], d[0]) < min(a[0], b[0]):
                continue
            if max(a[1], b[1]) < min(c[1], d[1]) or max(c[1], d[1]) < min(a[1], b[1]):
                continue
            det = u[0] * v[1] - u[1] * v[0]
            if abs(det) < 1e-30:
                continue
            w = c - a
            s = (w[0] * v[1] - w[1] * v[0]) / det
            t = (w[0] * u[1] - w[1] * u[0]) / det
            if 0 < s < 1 and 0 < t < 1:
                height = a[2] + s * u[2] - c[2] - t * v[2]
                rows.append((float(i), float(j), s, t, height, det))
    return rows


def camera(seed=0):
    """Deterministic proper rotation, independent of the orbit data."""
    q, r = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def diagram(xyz, seed=0):
    xyz = np.asarray(xyz, dtype=float)
    if np.linalg.norm(xyz[-1] - xyz[0]) < 1e-14:
        xyz = xyz[:-1]
    # Positive diagonal scaling is an ambient isotopy; it improves conditioning.
    scales = np.ptp(xyz, axis=0)
    points = (xyz - np.mean(xyz, axis=0)) / scales @ camera(seed)
    rows = intersections(points)
    events = sorted((r[k] + r[k+2], i, k) for i, r in enumerate(rows) for k in (0, 1))
    labels = {(i, k): e for e, (_, i, k) in enumerate(events)}
    pd = []; n = len(events)
    for i, r in enumerate(rows):
        under, over = (1, 0) if r[4] > 0 else (0, 1)
        u = labels[i, under]; o = labels[i, over]
        cross_under_over = r[5] if under == 0 else -r[5]
        pd.append([(u-1) % n, o, u, (o-1) % n] if cross_under_over < 0
                  else [(u-1) % n, (o-1) % n, u, o])
    link = Link(pd) if pd else Link([])
    if pd:
        assert len(link.link_components) == 1
        assert link.is_planar()
    return link, dict(seed=seed, crossings=len(rows), pd_code=pd,
        minimum_crossing_height=min((abs(r[4]) for r in rows), default=None),
        minimum_vertex_fraction=min((min(r[2],1-r[2],r[3],1-r[3]) for r in rows), default=None),
        minimum_event_spacing=min((events[i+1][0]-events[i][0] for i in range(n-1)),default=None)), points, rows


def polynomial_coefficients(link):
    """Alexander coefficients, lowest degree zero, positive leading coefficient."""
    if not link.crossings:
        return [1]
    v = sp.Matrix(link.seifert_matrix()); t = sp.Symbol('t')
    poly = sp.Poly((t*v-v.T).det(method='domain-ge'), t)
    terms = poly.terms(); low = min(e[0] for e, _ in terms)
    top = max(e[0] for e, _ in terms)
    coefficients = [int(poly.nth(i)) for i in range(low, top+1)]
    if coefficients[-1] < 0:
        coefficients = [-c for c in coefficients]
    assert abs(sum(coefficients)) == 1
    assert coefficients == coefficients[::-1]
    return coefficients


def canonical_gauss(link):
    """Diagram code modulo starting point and traversal orientation, not isotopy.

    Crossing signs retain chirality. Equal codes identify the same oriented
    planar crossing data; unequal codes do not establish inequivalent knots.
    """
    if not link.crossings:
        return ()
    assert len(link.link_components) == 1
    seq = [(id(c.crossing), int(c.is_over_crossing()), c.crossing.sign)
           for c in link.link_components[0]]
    versions = []
    for oriented in (seq, seq[::-1]):
        for shift in range(len(seq)):
            labels = {}; code = []
            for key, over, sign in oriented[shift:] + oriented[:shift]:
                if key not in labels: labels[key] = len(labels)
                code.append((labels[key], over, sign))
            versions.append(tuple(code))
    return min(versions)


def fox_coefficients(link):
    """Independent Alexander matrix from oriented Wirtinger relations."""
    if not link.crossings:
        return [1]
    seq = link.link_components[0]; n = len(link.crossings)
    under = [i for i,c in enumerate(seq) if c.is_under_crossing()]
    arc = lambda i: (int(np.searchsorted(under, i, side='right'))-1) % n
    t=sp.Symbol('t'); matrix=sp.zeros(n,n)
    for row, pos in enumerate(under):
        crossing=seq[pos].crossing
        over=next(i for i,c in enumerate(seq) if c.crossing is crossing and c.is_over_crossing())
        matrix[row,arc(over)] += 1-t
        if crossing.sign>0:
            matrix[row,(row-1)%n] += t; matrix[row,row] -= 1
        else:
            matrix[row,(row-1)%n] -= 1; matrix[row,row] += t
    poly=sp.Poly(matrix[:n-1,:n-1].det(method='domain-ge'),t)
    lo=min(e[0] for e,c in poly.terms());hi=poly.degree()
    values=[int(poly.nth(i)) for i in range(lo,hi+1)]
    if values[-1]<0:values=[-x for x in values]
    return values


def identify_summand(link, max_diagrams=3000):
    """Match a diagram to torus references using explicit type-III moves.

    This is a sufficient identification test, not a general knot recognizer.
    The path is retained for replay. Mirror images remain distinguished.
    """
    from collections import deque
    from spherogram.links import simplify
    targets={}
    for p,q in ((2,3),(2,5),(3,5)):
        k=Link(f'T({p},{q})')
        for mirrored in (False,True):
            name=f'T({p},{q})'+(' mirror' if mirrored else '')
            targets[canonical_gauss(k.mirror() if mirrored else k)]=name
    first=link.copy();queue=deque([(first,[])]);seen={canonical_gauss(first)}
    while queue and len(seen)<=max_diagrams:
        current,path=queue.popleft();code=canonical_gauss(current)
        if code in targets:
            return dict(type=targets[code],start_pd=link.PD_code(),type_III_moves=path,
                        final_pd=current.PD_code(),visited_diagrams=len(seen))
        for index in range(len(simplify.possible_type_III_moves(current))):
            k=current.copy();triple=simplify.possible_type_III_moves(k)[index]
            descriptor=[(k.crossings.index(c),e) for c,e in triple]
            simplify.reidemeister_III(k,triple);k._rebuild()
            code=canonical_gauss(k)
            if code not in seen:
                seen.add(code);queue.append((k,path+[descriptor]))
    return dict(type='unresolved',start_pd=link.PD_code(),visited_diagrams=len(seen))


def analyze(xyz, seed=0, simplify='global'):
    link, record, points, crossings = diagram(xyz, seed)
    if link.crossings:
        import random
        random.seed(20261101 + seed)
        link.simplify(simplify)
    record['simplified_crossings'] = len(link.crossings)
    record['simplified_pd_code'] = link.PD_code()
    record['alexander_coefficients'] = polynomial_coefficients(link)
    record['fox_alexander_coefficients'] = fox_coefficients(link)
    assert record['alexander_coefficients'] == record['fox_alexander_coefficients']
    record['determinant'] = abs(sum(c*(-1)**i for i,c in enumerate(record['alexander_coefficients'])))
    record['unknot_by_reidemeister_reduction'] = not bool(link.crossings)
    record['summands'] = [identify_summand(k) for k in link.deconnect_sum()] if link.crossings else []
    record['knot_type'] = ' # '.join(sorted(s['type'] for s in record['summands'])) or 'unknot'
    return record, link
