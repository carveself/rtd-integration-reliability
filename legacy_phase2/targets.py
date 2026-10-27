"""Goal stopping rules; no reference solution is visible to these algorithms."""
import numpy as np

def coefficient_of_variation(y):return float(np.std(y,ddof=1)/np.mean(y))


def cv_envelope(y,eps):
    """Conditional norm bound: valid only IF eps bound the actual per-item errors."""
    y=np.asarray(y,float);eps=np.asarray(eps,float)
    if y.ndim!=1 or eps.shape!=y.shape or len(y)<2 or np.any(eps<0):raise ValueError('Invalid samples/errors')
    m=np.mean(y);s=np.std(y,ddof=1);dm=np.mean(eps);ds=np.linalg.norm(eps)/np.sqrt(len(y)-1)
    if m<=dm:return 0.,float('inf')
    return max(0.,s-ds)/(m+dm),(s+ds)/(m-dm)


def current_arrays(states):return np.array([s.value for s in states]),np.array([s.error for s in states])


def goal_cv(states,cfg):
    trajectory=[]
    for iteration in range(20000):
        y,e=current_arrays(states);c=coefficient_of_variation(y);lo,hi=cv_envelope(y,e)
        radius=max(c-lo,hi-c)
        near=hi<cfg['cv_near_zero']
        accepted=radius<=cfg['cv_absolute_budget'] if near else radius<=cfg['cv_relative_budget']*lo
        trajectory.append(dict(iteration=iteration,evaluations=sum(s.evaluations for s in states),estimate=c,lower=lo,upper=hi))
        if accepted:return dict(estimate=c,lower=lo,upper=hi,mode='absolute' if near else 'relative',trajectory=trajectory)
        # Largest contribution to the L2 perturbation bound; deterministic tie handling.
        states[int(np.argmax(e*e))].refine()
    raise RuntimeError('Goal-CV refinement limit')


def goal_decision(states,cfg):
    trajectory=[]
    for iteration in range(20000):
        y,e=current_arrays(states);lower=np.maximum(0.,y-e);upper=y+e;chosen=int(np.argmax(y))
        bound=float(1-lower[chosen]/np.max(upper))
        possible=np.flatnonzero(upper>=np.max(lower)).tolist()
        trajectory.append(dict(iteration=iteration,evaluations=sum(s.evaluations for s in states),chosen=chosen,regret_estimate=bound,possible=possible))
        if bound<=cfg['decision_regret_budget']:
            return dict(chosen=chosen,regret_estimate=bound,possible=possible,trajectory=trajectory)
        targets=set(possible+[chosen]);target=max(sorted(targets),key=lambda i:e[i])
        states[target].refine()
    raise RuntimeError('Goal-decision refinement limit')


def budget_lockstep(states,cfg,task):
    """Same target envelopes as goal, with uniform per-item relative tolerances.

    Strong scheduling baseline, introduced as a labeled post-review supplement.
    No development or reference targets are consulted by the stopping rule.
    """
    levels=[*cfg['baseline_rtol_grid'],1e-7,1e-8,1e-9,1e-10]
    for stage,rtol in enumerate([None,*levels]):
        if rtol is not None:
            for s in states:s.finish(rtol)
        y,e=current_arrays(states)
        if task=='cv':
            value=coefficient_of_variation(y);lo,hi=cv_envelope(y,e)
            radius=max(value-lo,hi-value);near=hi<cfg['cv_near_zero']
            accepted=radius<=cfg['cv_absolute_budget'] if near else radius<=cfg['cv_relative_budget']*lo
            if accepted:return dict(estimate=value,lower=lo,upper=hi,stage=stage)
        else:
            chosen=int(np.argmax(y));lower=np.maximum(0.,y-e);upper=y+e
            bound=float(1-lower[chosen]/max(upper))
            if bound<=cfg['decision_regret_budget']:
                return dict(chosen=chosen,regret_estimate=bound,stage=stage)
    raise RuntimeError('Budget-lockstep target not met')
