"""An attracting cycle and explicit equilibrium neighborhoods avoiding its basin.

The central region is exactly linear. A real eigenbasis gives an expanding
coordinate and a contracting rotating plane. Certified exit boxes enter the
outer equilibria's invariant Lyapunov ellipsoids. Together with an enclosed
periodic sink this establishes hiddenness, without finite-time basin sampling.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import solve_continuous_lyapunov
from flint import arb, arb_mat, acb, ctx
from chua_interval_map import col, bounds
from certify_chua_coexisting_cycles import ROOT, OUT, ReturnMap, certify


def positive_definite(M):
    assert all(M[i,j].overlaps(M[j,i]) for i in range(M.nrows()) for j in range(M.ncols()))
    minors = [arb_mat([[M[i,j] for j in range(k)] for i in range(k)]).det() for k in range(1,4)]
    assert all(d > 0 for d in minors), "Positive definiteness rejected"
    return [bounds(d) for d in minors]


def lyapunov(A):
    An = np.array([[float(A[i,j].mid()) for j in range(3)] for i in range(3)])
    Pn = solve_continuous_lyapunov(An.T, -np.eye(3))
    # Decimal entries designate a fixed rational symmetric matrix, not an
    # assumed exact solution of the floating-point Lyapunov equation.
    Pn = (Pn+Pn.T)/2
    entries = [[format(Pn[i,j], '.14g') for j in range(3)] for i in range(3)]
    P = arb_mat([[arb(x) for x in row] for row in entries])
    Q = -(A.transpose()*P+P*A)
    certificate = dict(P=entries, P_minors=positive_definite(P), Q_minors=positive_definite(Q))
    return P, certificate


def eigenbasis(A, alpha, beta, m0):
    a1 = 1+alpha*(1+m0)
    a2 = beta+alpha*m0
    a3 = alpha*beta*(1+m0)
    poly = lambda x: x**3+a1*x*x+a2*x+a3
    deriv = lambda x: 3*x*x+2*a1*x+a2
    t = arb("2.712")
    for _ in range(12):
        t = (t-poly(t)/deriv(t)).mid()
    box = arb(t, arb("1e-50"))
    N = t-poly(t)/deriv(box)
    assert box.contains_interior(N) and box > 0
    lam = N
    sigma = (-a1-lam)/2
    omega = (-a3/lam-sigma*sigma).sqrt()
    assert sigma < 0 and omega > 0
    mu = acb(sigma, omega)
    v = col([1, 1+m0+lam/alpha, -beta*(1+m0+lam/alpha)/lam])
    yc = 1+m0+mu/alpha
    zc = -beta*yc/mu
    S = arb_mat([[1, 1, 0], [v[1,0], yc.real, yc.imag], [v[2,0], zc.real, zc.imag]])
    assert not S.det().contains(0)
    return S, dict(unstable=bounds(lam), stable_real=bounds(sigma), stable_imaginary=bounds(omega),
                   eigenbasis=[[bounds(S[i,j]) for j in range(3)] for i in range(3)])


def outer_path(f, A, c, start, time):
    work = [(arb(0), time, 0)]
    count = 0
    while work:
        left, right, depth = work.pop()
        state, _ = f.flow(A,c,start,left.union(right))
        speed = (A*(state-c))[0,0]
        if state[0,0] > 1 or (left.is_zero() and speed > 0):
            count += 1
            continue
        assert depth < 22, "Outer trajectory crosses the breakpoint or cannot be resolved"
        mid = (left+right)/2
        work.extend([(left,mid,depth+1),(mid,right,depth+1)])
    return count


def run(digits):
    ctx.dps = digits
    parameters = (arb(42)/5, arb(12), -arb(6)/5, -arb(1)/20)
    seed = json.loads((OUT/"hidden-cycle-seed.json").read_text())
    orbit = certify(seed, digits, parameters=parameters)
    f = ReturnMap(seed,digits,parameters)
    cycle_box=[arb(v['ball']) for v in orbit['root']]
    _,_,cycle_arcs=f.evaluate(cycle_box)
    amplitude=arb(0)
    for arc in cycle_arcs:
        A,c,_=f.matrices(None,arc['region'])
        for k in range(256):
            t=(arb(k)/256).union(arb(k+1)/256)*arc['duration']
            state,_=f.flow(A,c,arc['start'],t)
            amplitude=amplitude.max(abs(state[0,0]))
    assert amplitude.abs_upper()<12
    A0,c0,_ = f.matrices(None,0)
    A1,c1,_ = f.matrices(None,1)
    S, eig = eigenbasis(A0,*parameters[:3])
    P, lyap = lyapunov(A1)
    margin = c1[0,0]-1
    ellipsoid_level = margin*margin/(2*P.inv()[0,0])
    assert ellipsoid_level > 0
    # ||stable_coordinates||_2 <= delta contracts until the unstable coordinate
    # reaches +/- 1/2. The rectangular exit enclosure below contains the disk.
    delta = arb("1e-5")
    assert arb(1)/2+delta < 1
    start = S*col([arb(1)/2, arb(0,delta), arb(0,delta)])
    T = f.crossing(A0,c0,start,1,.25555)
    exit_state,_ = f.flow(A0,c0,start,T)
    exit_state[0,0] = 1
    central = f.verify_arc(A0,c0,dict(duration=T,start=start,region=0,target=1))
    flight = arb(120)
    path_cells = outer_path(f,A1,c1,exit_state,flight)
    end,_ = f.flow(A1,c1,exit_state,flight)
    offset = end-c1
    V = (offset.transpose()*P*offset)[0,0]
    assert V < ellipsoid_level, "Exit box not captured by the invariant ellipsoid"
    # Explicit Euclidean balls contained in the certified coordinate domain
    # and in each outer invariant ellipsoid. Frobenius/row norms are bounds.
    inv = S.inv()
    inv_frob = sum(inv[i,j]**2 for i in range(3) for j in range(3)).sqrt()
    origin_radius = (delta/(2*inv_frob)).lower()
    p_norm = max(sum(abs(P[i,j]) for j in range(3)).abs_upper() for i in range(3))
    outer_radius = (ellipsoid_level/(2*p_norm)).sqrt().lower()
    # Reject a tube wide enough to contain states outside the central region.
    bad = S*col([arb(1)/2,arb(0,arb(2)),arb(0,arb(2))])
    assert not (bad[0,0] > 0 and bad[0,0] < 1)
    report = dict(passed=True,working_digits=digits,
        parameters=dict(alpha="42/5",beta="12",m0="-6/5",m1="-1/20",gamma="0"),
        periodic_orbit=orbit, periodic_x_amplitude_upper=bounds(amplitude.abs_upper()),origin_spectrum=eig, outer_lyapunov=lyap,
        outer_equilibrium=[bounds(x) for x in c1.entries()],
        stable_coordinate_radius=str(delta), unstable_exit_coordinate="1/2",
        central_exit_time=bounds(T), central_path=central, outer_path_cells=path_cells,
        outer_flight_time="120", captured_V=bounds(V), ellipsoid_level=bounds(ellipsoid_level),
        excluded_ball_radii=dict(origin=bounds(origin_radius), each_outer_equilibrium=bounds(outer_radius)),
        wide_exit_control_rejected=True,
        scope="One hidden attracting periodic orbit and three explicit equilibrium neighborhoods disjoint from its basin. This is an ideal gamma=0 model; no maximum attractor count or chaotic-attractor claim.")
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--digits',type=int,default=80)
    args=parser.parse_args()
    report=run(args.digits)
    (OUT/f'hidden-cycle-certificate-{args.digits}.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Hidden cycle certified:',report['periodic_orbit']['period']['ball'],flush=True)
    print('Excluded equilibrium balls:',report['excluded_ball_radii'],flush=True)


if __name__ == '__main__':
    main()
