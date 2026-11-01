"""Compose enclosed return jets in coordinates following the saddle branch."""

import numpy as np
from flint import arb, arb_mat
from chua_interval_map import col, eye, upper
from certify_chua_saddle import submatrix, sum_matrices


def connect(
    f,
    a,
    v,
    beta,
    beta0,
    p0,
    pb0,
    pb,
    pbb,
    V,
    W,
    iterations,
    center_records=None,
    verify_path=False,
):
    delta_beta = beta - beta0
    p = p0 + pb0 * delta_beta + pbb * (delta_beta * delta_beta / 2)
    state = p + V * col([a, v])
    D = arb_mat([[V[i, 0], V[i, 1], pb[i, 0]] for i in range(2)])
    H = [arb_mat(3, 3) for i in range(2)]
    for i in range(2):
        H[i][2, 2] = pbb[i, 0]
    moving=getattr(f,'moving_frame',False)
    Q=arb_mat(V) if moving else eye(2)
    if moving:
        Qi=Q.inv();D=Qi*D
        H=[sum_matrices([H[j]*Qi[i,j] for j in range(2)]) for i in range(2)]
    deviations = [a - a.mid(), v - v.mid(), beta - beta.mid()]
    records = []
    for k in range(iterations):
        try:
            out, Jlocal, Hlocal, arcs = f.evaluate(
                list(state.entries()), beta, verify_path=verify_path
            )
        except ValueError as exc:
            raise ValueError(f'Return {k+1}/{iterations}, state {state}: {exc}') from exc
        if moving:
            # Fixed rational frames for this enclosure, chosen by floating QR.
            # Their numerical accuracy is irrelevant: all inverses and changes
            # of coordinates below are enclosed in Arb. No intervals are lost.
            product=arb_mat([[Jlocal[i,j] for j in range(2)] for i in range(2)])*Q
            numeric=np.array([[float(product[i,j].mid()) for j in range(2)] for i in range(2)])
            qn,_=np.linalg.qr(numeric)
            Qnew=arb_mat([[arb(format(x,'.17g')) for x in row] for row in qn]);Qi=Qnew.inv()
            R=arb_mat([[Q[0,0],Q[0,1],0],[Q[1,0],Q[1,1],0],[0,0,1]])
            converted=[R.transpose()*h*R for h in Hlocal]
            Hlocal=[sum_matrices([converted[j]*Qi[i,j] for j in range(2)]) for i in range(2)]
            Jlocal=Qi*Jlocal*R
        else:Qnew=eye(2)
        T = arb_mat(
            [[D[i, j] if i < 2 else int(j == 2) for j in range(3)] for i in range(3)]
        )
        newH = [
            T.transpose() * Hlocal[i] * T
            + sum_matrices([H[j] * Jlocal[i, j] for j in range(2)])
            for i in range(2)
        ]
        D = Jlocal * T
        H = newH
        state = out
        physical_H=[sum_matrices([H[j]*Qnew[i,j] for j in range(2)]) for i in range(2)]
        if center_records is not None:
            center = center_records[k]
            enclosure = center["state"]
            for j in range(3):
                enclosure += col([center["D"][i, j] for i in range(2)]) * deviations[j]
                for m in range(3):
                    enclosure += col([physical_H[i][j, m] for i in range(2)]) * (
                        deviations[j] * deviations[m] / 2
                    )
            state = col([state[i, 0].intersection(enclosure[i, 0]) for i in range(2)])
            center_D=Qnew.inv()*center['D']
            for i in range(2):
                for j in range(3):
                    derivative = center_D[i, j]
                    for m in range(3):
                        derivative += H[i][j, m] * deviations[m]
                    D[i, j] = D[i, j].intersection(derivative)
        records.append(
            dict(
                state=state,
                D=Qnew*D,
                H=physical_H,
                path_cells=[
                    z.get("path_certificate", {}).get("accepted_cells") for z in arcs
                ],
            )
        )
        Q=Qnew
    value = W * (state - p)
    D = Q*D
    H = [arb_mat(h) for h in physical_H]
    for i in range(2):
        D[i, 2] -= pb[i, 0]
        H[i][2, 2] -= pbb[i, 0]
    jac = W * D
    hess = [sum_matrices([H[i] * W[k, i] for i in range(2)]) for k in range(2)]
    return value, jac, hess, records
