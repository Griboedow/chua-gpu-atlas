"""Actual-orbit proximity audit for source193; no inference from FTLE alone."""
from pathlib import Path
import json
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree

from rayleigh_pair import ROOT, field_jac, save


def reference_cloud(states, period, copies=True):
    states = states[:-1]
    variants = [states, -states, states[:, [2,3,0,1]], -states[:, [2,3,0,1]]] if copies else [states]
    vertices = np.stack(variants)
    return vertices, cKDTree(vertices.reshape(-1, 4)), period


def distances(points, reference):
    vertices, tree, period = reference
    n = vertices.shape[1]
    _, nearest = tree.query(points)
    copies, indices = nearest//n, nearest%n
    best = np.full(len(points), np.inf)
    phases = np.zeros(len(points))
    for offset in (-1, 0):
        left = (indices+offset)%n
        a = vertices[copies, left]
        direction = vertices[copies, (left+1)%n]-a
        fraction = np.clip(np.sum((points-a)*direction, axis=1)/np.sum(direction*direction, axis=1), 0., 1.)
        d = np.linalg.norm(points-a-fraction[:, None]*direction, axis=1)
        replace = d < best
        best[replace] = d[replace]
        phases[replace] = ((left[replace]+fraction[replace])/n*period)%period
    return best, copies, phases


def intervals(times, mask):
    padded = np.r_[False,mask,False]
    starts = np.flatnonzero(np.diff(padded.astype(int)) == 1)
    stops = np.flatnonzero(np.diff(padded.astype(int)) == -1)-1
    return [[float(times[a]),float(times[b])] for a,b in zip(starts,stops)]


def main():
    started = time.perf_counter()
    deadline = started+150
    scout = np.load(ROOT/'scout.npz')
    refined = np.load(ROOT/'scout-refined.npz')
    row = 193
    match = np.flatnonzero(refined['source_indices'] == row)
    if len(match) != 1:
        raise RuntimeError('Expected one source193 refinement row')
    rr = int(match[0]);p = scout['parameters'][row].copy()
    initial = scout['initial_states'][row].copy()
    scout_time = float(scout['elapsed_sigma'][row])
    refined_time = scout_time+float(refined['elapsed_sigma'][rr])
    checkpoints = [scout_time,refined_time]
    quad_data = np.load(ROOT/'quadrature-cycle.npz')
    quad_meta = json.loads((ROOT/'quadrature-cycle.json').read_text())
    if not np.allclose(p,quad_data['parameters'],rtol=1e-14,atol=0):
        raise RuntimeError('Scout and quadrature reference parameters differ')
    quad = reference_cloud(quad_data['states'], quad_meta['period_sigma'])
    weak = json.loads((ROOT/'weak-coupling-modal-controls.json').read_text())
    if not np.allclose(p,[weak['parameters']['mu'],weak['parameters']['kappa']],rtol=1e-14,atol=0):
        raise RuntimeError('Scout and modal reference parameters differ')
    references = {'quadrature':quad}
    for control in weak['controls']:
        amplitude = control['section_velocity'];sign = 1 if control['mode']=='symmetric' else -1
        seed = np.array([0.,amplitude,0.,sign*amplitude]);period=control['period_sigma']
        sol = solve_ivp(lambda t,y:field_jac(y,p)[0],(0,period),seed,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.02,dense_output=True)
        references[control['mode']] = reference_cloud(sol.sol(np.linspace(0,period,6001)).T,period)

    times, states = [np.array([0.])], [initial[None,:]]
    y=initial.copy();t=0.;stop=30000.;termination='completed to 30000 sigma'
    while t < stop:
        if time.perf_counter() > deadline:
            termination='CPU budget reached; saved completed chunks'
            break
        end=min(stop,t+1000)
        sample=np.arange(t+1,end+.1,1.)
        sample=np.unique(np.r_[sample,[c for c in checkpoints if t<c<=end]])
        sol=solve_ivp(lambda tau,z:field_jac(z,p)[0],(t,end),y,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.1,t_eval=sample)
        if not sol.success:
            termination=sol.message;break
        times.append(sol.t);states.append(sol.y.T);y=sol.y[:,-1];t=end
        if int(t)%5000==0:print({'sigma':t,'wall_seconds':time.perf_counter()-started},flush=True)
    times=np.concatenate(times);states=np.concatenate(states)
    arrays={};summaries={}
    for name,reference in references.items():
        d,copy,phase=distances(states,reference)
        arrays[name+'_distance']=d;arrays[name+'_copy']=copy;arrays[name+'_phase_sigma']=phase
        summaries[name]={'minimum_distance':float(min(d)),'minimum_time_sigma':float(times[np.argmin(d)]),
                         'final_distance':float(d[-1]),'below_0.01_intervals_sigma':intervals(times,d<.01),
                         'below_0.05_intervals_sigma':intervals(times,d<.05),'below_0.1_intervals_sigma':intervals(times,d<.1)}
    gpu_points=np.vstack([initial,scout['final_states'][row],refined['final_states'][rr]])
    gpu_summary=[]
    for i,label in enumerate(['original start','scout endpoint','refined endpoint']):
        item={'label':label,'sigma':[0.,scout_time,refined_time][i],'state':gpu_points[i].tolist()}
        for name,reference in references.items():item[name+'_distance']=float(distances(gpu_points[i:i+1],reference)[0][0])
        at=np.flatnonzero(abs(times-item['sigma'])<1e-8)
        item['cpu_gpu_state_difference']=float(np.linalg.norm(states[at[0]]-gpu_points[i])) if len(at) else None
        gpu_summary.append(item)
    closest=arrays['quadrature_distance']
    coarse=reference_cloud(quad_data['states'][::2],quad_meta['period_sigma'])
    coarse_distance=distances(states,coarse)[0]
    fit_mask=(times>=1000)&(times<=12000)&(closest<.03)
    growth_fit={'status':'insufficient samples in prescribed finite-distance interval'}
    if np.count_nonzero(fit_mask)>100:
        coefficients=np.polyfit(times[fit_mask],np.log(closest[fit_mask]),1)
        residual=np.log(closest[fit_mask])-np.polyval(coefficients,times[fit_mask])
        variance=np.sum((np.log(closest[fit_mask])-np.mean(np.log(closest[fit_mask])))**2)
        growth_fit={'status':'finite geometric departure fit','time_interval_sigma':[1000.,12000.],
                    'interval_selection':'Post hoc chosen after inspecting the proximity/departure data; exploratory agreement, not a preregistered independent test.',
                    'fit_method':'Unweighted ordinary least squares of natural log phase-minimized Euclidean polygon distance versus canonical sigma, using all saved samples in the stated time and distance filters.',
                    'distance_upper_cutoff':.03,'sample_count':int(np.count_nonzero(fit_mask)),
                    'log_distance_growth_rate':float(coefficients[0]),
                    'R_squared':float(1-np.sum(residual**2)/variance),
                    'reference_unstable_floquet_exponent':float(max(quad_data['floquet_exponents'])),
                    'scope':'Corroborates finite saddle departure; does not establish an attracting set or infinite-time shadowing.'}
    save(ROOT/'saddle-shadowing',dict(parameters=p,source_row=row,refined_row=rr,time_unit='canonical sigma',
         input_files=['scout.npz','scout-refined.npz','quadrature-cycle.npz','weak-coupling-modal-controls.json'],
         initial_state=initial,times=times,states=states,**arrays,proximity_summary=summaries,gpu_checkpoint_comparisons=gpu_summary,
         actual_end_sigma=float(times[-1]),requested_end_sigma=stop,termination=termination,
         method='Independent DOP853 state integration, rtol1e-11,atol1e-13,max_step0.1; samples every1sigma plus exact GPU checkpoint times.',
         distance_definition='Unweighted Euclidean distance in (x1,v1,x2,v2), minimized over adjacent polygon segments near the nearest vertex of each sampled exact cycle. Quadrature/modal references include exchange and global-sign copies.',
         quadrature_polygon_resolution_comparison_max_abs=float(max(abs(coarse_distance-closest))),
         symmetry_copy_order=['identity','global sign','exchange oscillators','exchange plus global sign'],
         geometric_departure_fit=growth_fit,
         distance_scope='Geometric phase-minimized proximity, not a theorem guaranteeing shadowing for infinite time. Interval endpoints have roughly1sigma sampling resolution.',
         reference_periods={name:reference[2] for name,reference in references.items()},wall_seconds=time.perf_counter()-started,
         script='research/rayleigh_saddle_shadowing.py',scope='Check actual proximity and departure; no attractor or novelty claim.'))
    print('Summaries',summaries,flush=True);print('GPU checkpoints',gpu_summary,flush=True)
    print('Termination',termination,'wall',time.perf_counter()-started,flush=True)


if __name__=='__main__':main()
