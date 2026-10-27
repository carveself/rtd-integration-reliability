"""Posthoc replay of all 308 OLD cases at identical phase2 target budgets.

This is a regression/development audit, never fresh validation. Old reference
values are read only after stopping, exactly as in the original phase2 study.
"""
from pathlib import Path
import sys,json,csv,time
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'vendor'),str(ROOT/'legacy_phase2')]
from models import model
from resonance import discover
from seeded import SeededIntegral
from targets import goal_cv,goal_decision,coefficient_of_variation
from run_study import write_csv


def main():
    old=json.loads((ROOT/'legacy_phase2/protocol.json').read_text())
    new=json.loads((ROOT/'protocol.json').read_text());cfg={**new,**old}
    refs={(r['context'],int(r['item'])):float(r['J_A_m2'])
          for r in csv.DictReader((ROOT/'legacy_phase2/references.csv').open())}
    contexts=[]
    for split in ['development','holdout','confirmation']:
        for c in old[split]:
            heights=np.random.default_rng(c['seed']).normal(.28,old['height_sigma_eV'],old['samples'])
            contexts.append((c,'cv',[model(c,height=float(h)) for h in heights]))
    for c in old['designs']:contexts.append((c,'decision',[model(c,well=w) for w in c['wells_nm']]))
    rows=[];details=[]
    for c,task,devices in contexts:
        for method in ('uniform','spectral','pole'):
            start=time.perf_counter();states=[];residuals=0;setup=0.
            for d in devices:
                disc=None
                if method!='uniform':
                    disc=discover(d,c['bias_V'],c['temperature_K'],refine=method=='pole')
                    residuals+=disc['residual_calls'];setup+=disc['setup_seconds']
                states.append(SeededIntegral(d,c['bias_V'],c['temperature_K'],cfg,disc,method))
            outcome=goal_cv(states,cfg) if task=='cv' else goal_decision(states,cfg)
            seconds=time.perf_counter()-start
            # Consult old reference only now, after the algorithm has stopped.
            truth=np.array([refs[c['id'],i] for i in range(len(devices))])
            y=np.array([s.value for s in states]);e=np.array([s.error for s in states]);err=abs(y-truth)
            if task=='cv':
                value=coefficient_of_variation(y);target=coefficient_of_variation(truth)
                target_error=abs(value-target)/target;passed=target_error<=cfg['cv_relative_budget']
            else:
                target_error=1-truth[outcome['chosen']]/max(truth);passed=target_error<=cfg['decision_regret_budget']
            rows.append({'context':c['id'],'task':task,'method':method,'items':len(states),'energy_evaluations':sum(s.evaluations for s in states),
                         'pole_residual_calls':residuals,'setup_seconds':setup,'seconds':seconds,
                         'undercovered_items':int(sum(err>e)),'maximum_underestimate_ratio':float(max(err/np.maximum(e,1e-30))),
                         'target_error':float(target_error),'target_pass':bool(passed)})
            details.extend({'context':c['id'],'task':task,'method':method,'item':i,'J_A_m2':float(y[i]),
                            'estimated_error_A_m2':float(e[i]),'actual_error_A_m2':float(err[i]),
                            'underestimate_ratio':float(err[i]/max(e[i],1e-30))} for i in range(len(y)))
            print(c['id'],method,'undercoverage',sum(err>e),'target pass',passed,flush=True)
    out=ROOT/'results_legacy';out.mkdir(exist_ok=True)
    write_csv(out/'benchmark.csv',rows);write_csv(out/'sample_errors.csv',details)

if __name__=='__main__':main()
