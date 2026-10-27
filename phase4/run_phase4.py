"""Fresh validation and posthoc legacy regression; no reference-driven stopping."""
from pathlib import Path
import sys,csv,json,time,hashlib,warnings
import numpy as np
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'phase4'
sys.path[:0]=[str(P/'src'),str(ROOT/'src'),str(ROOT/'vendor'),str(ROOT/'legacy_phase2'),str(ROOT)]
from selective import SelectiveIntegral
from resonance import make_device,discover,bounds
from seeded import SeededIntegral
from models import model
from targets import goal_cv,goal_decision,coefficient_of_variation
from run_study import independent_reference,write_csv
from highprecision import MPChain,integrate_window
from reference import direct_banded
from negf import supply,current_prefactor
from scipy.integrate import quad,IntegrationWarning

METHODS=['uniform','spectral','pole','selective']


def candidate(dev,case,cfg,method,finish=True):
    start=time.perf_counter();disc=None
    if method=='selective':st=SelectiveIntegral(dev,case['bias'],case['temperature'],cfg)
    else:
        if method!='uniform':disc=discover(dev,case['bias'],case['temperature'],refine=method=='pole')
        st=SeededIntegral(dev,case['bias'],case['temperature'],cfg,disc,method)
    if finish:st.finish(cfg['rtol'])
    return st,disc,time.perf_counter()-start


def hybrid_reference(dev,case,disc,scale):
    start=time.perf_counter();bias=case['bias'];temp=case['temperature'];lo,hi=bounds(dev,bias,temp)
    poles=sorted([z for z in disc['poles'] if -2*z.imag<1e-5],key=lambda z:z.real)
    windows=[(max(lo,z.real-64*(-2*z.imag)),min(hi,z.real+64*(-2*z.imag)),z) for z in poles]
    if any(a[1]>=b[0] for a,b in zip(windows[:-1],windows[1:])):raise RuntimeError('Overlapping precision windows')
    a=MPChain(dev,bias,temp,50);b=MPChain(dev,bias,temp,80);va=0.;vb=0.;eb=0.
    for l,r,z in windows:
        v,_,_=integrate_window(a,l,r,z.real,-2*z.imag,(10,20),1e-9);va+=v
        v,e,_=integrate_window(b,l,r,z.real,-2*z.imag,(16,32),1e-10);vb+=v;eb+=e
    edges=np.unique(np.r_[np.linspace(lo,hi,int(np.ceil((hi-lo)/.004))+1),[x for l,r,_ in windows for x in (l,r)]])
    pref=current_prefactor(dev.mass);calls=0
    def f(en):
        nonlocal calls
        calls+=1
        return pref*direct_banded(dev,en,bias)*float(supply(en,bias,temp,dev.mu_eV))
    outside=0.;err=0.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always',IntegrationWarning)
        for l,r in zip(edges[:-1],edges[1:]):
            if any(a<(l+r)/2<b for a,b,_ in windows):continue
            v,e=quad(f,l,r,epsabs=max(1e-10,abs(scale)*1e-11)/len(edges),epsrel=1e-8,limit=200)
            outside+=v;err+=e
    value=outside+vb;delta=abs(va-vb)/abs(value)
    return value,{'MP50_calls':a.calls,'MP80_calls':b.calls,'MP_relative_agreement':delta,
                  'outside_calls':calls,'outside_warnings':len(caught),'MP_windows':len(windows),
                  'MP_reported_error':err+eb,'MP_seconds':time.perf_counter()-start,
                  'accepted':len(windows)>0 and not caught and delta<1e-6 and err+eb<1e-6*abs(value)}


def reference(dev,case,cfg,scale):
    start=time.perf_counter();disc=discover(dev,case['bias'],case['temperature'])
    q,qe,calls,messages=independent_reference(dev,case,cfg,disc,scale)
    row={'double_J_A_m2':q,'double_quad_error':qe,'double_banded_calls':calls,'double_warnings':json.dumps(messages)}
    if messages:
        truth,mp=hybrid_reference(dev,case,disc,q);row.update(mp);row['reference_type']='hybrid_MP_repair'
    else:
        fine=SeededIntegral(dev,case['bias'],case['temperature'],cfg,mode='uniform',cap=.00025).finish(1e-7)
        delta=abs(fine.value/q-1);truth=q;row.update(fine_J_A_m2=fine.value,fine_evaluations=fine.evaluations,
                                                 fine_relative_agreement=delta,accepted=delta<1e-6,
                                                 reference_type='banded_quad_plus_unseeded_fine')
    row.update(reference_J_A_m2=truth,reference_seconds=time.perf_counter()-start)
    return truth,row


def fresh(cfg):
    out=P/'results';rng=np.random.default_rng(cfg['seed']);cases=[]
    for c in cfg['contexts']:
        for i,h in enumerate(rng.normal(.28,cfg['height_sigma_eV'],cfg['draws_per_context'])):
            cases.append({**c,'sample':i,'height_eV':float(h)})
    (out/'cases.json').write_text(json.dumps(cases,indent=2)+'\n');rows=[];refs=[];timings=[];selections={}
    for c in cases:
        dev=make_device(c,c['height_eV']);saved={};times={m:[] for m in METHODS}
        for rep in range(cfg['timing_repetitions']):
            order=METHODS[rep:]+METHODS[:rep]
            for method in order:
                st,d,elapsed=candidate(dev,c,cfg,method)
                times[method].append(elapsed)
                timings.append({'id':c['id'],'sample':c['sample'],'method':method,'repeat':rep,'seconds':elapsed})
                if rep==0:saved[method]=(st,d)
        truth,audit=reference(dev,c,cfg,saved['uniform'][0].value)
        refs.append({'id':c['id'],'sample':c['sample'],**audit})
        for method,(st,disc) in saved.items():
            e=abs(st.value-truth);selection=st.selection if method=='selective' else {}
            rows.append({'id':c['id'],'sample':c['sample'],'method':method,'J_A_m2':st.value,'actual_abs_error_A_m2':e,
                         'relative_error':e/abs(truth),'estimated_error_A_m2':st.error,'covered':bool(e<=st.error),
                         'passed':bool(e<=cfg['rtol']*abs(truth)+cfg['absolute_floor_A_m2']),
                         'reference_accepted':audit['accepted'],'energy_evaluations':st.evaluations,
                         'coarse_evaluations':selection.get('coarse_evaluations',0),
                         'eigen_seeds':selection.get('eigen_seeds',0 if disc is None else len(disc['energies'])),
                         'retained_seeds':selection.get('retained_seeds',0 if disc is None else len(disc['energies'])),
                         'pole_residual_calls':0 if disc is None else disc['residual_calls'],
                         'seconds_median':float(np.median(times[method])),'seconds_min':min(times[method]),'seconds_max':max(times[method]),
                         'setup_seconds':selection.get('setup_seconds',0 if disc is None else disc['setup_seconds'])})
            if method=='selective':selections[f"{c['id']}/{c['sample']}"]=selection
        write_csv(out/'benchmark.csv',rows);write_csv(out/'references.csv',refs);write_csv(out/'timings.csv',timings)
        (out/'selection_trace.json').write_text(json.dumps(selections,indent=2)+'\n')
        print(c['id'],c['sample'],'reference',audit['accepted'],'selective seeds',saved['selective'][0].selection['retained_seeds'],flush=True)


def regression(cfg):
    old=json.loads((ROOT/'legacy_phase2/protocol.json').read_text());cc={**cfg,**old};rows=[];details=[]
    lookup={(r['context'],int(r['item'])):float(r['J_A_m2']) for r in csv.DictReader((ROOT/'legacy_phase2/references.csv').open())}
    tasks=[]
    for split in ['development','holdout','confirmation']:
        for c in old[split]:
            heights=np.random.default_rng(c['seed']).normal(.28,old['height_sigma_eV'],old['samples'])
            tasks.append((c,'cv',[model(c,height=float(h)) for h in heights]))
    for c in old['designs']:tasks.append((c,'decision',[model(c,well=w) for w in c['wells_nm']]))
    out=P/'results_regression'
    for c,task,devices in tasks:
        case={'bias':c['bias_V'],'temperature':c['temperature_K']}
        for method in METHODS:
            tick=time.perf_counter();states=[];coarse=0;seeds=0;retained=0
            for dev in devices:
                st,_,_=candidate(dev,case,cc,method,finish=False);states.append(st)
                if method=='selective':
                    coarse+=st.selection['coarse_evaluations'];seeds+=st.selection['eigen_seeds'];retained+=st.selection['retained_seeds']
            result=goal_cv(states,cc) if task=='cv' else goal_decision(states,cc)
            elapsed=time.perf_counter()-tick
            truth=np.array([lookup[c['id'],i] for i in range(len(states))]);y=np.array([s.value for s in states]);eps=np.array([s.error for s in states]);error=abs(y-truth)
            target_error=abs(coefficient_of_variation(y)/coefficient_of_variation(truth)-1) if task=='cv' else 1-truth[result['chosen']]/max(truth)
            rows.append({'id':c['id'],'task':task,'method':method,'items':len(states),'energy_evaluations':sum(s.evaluations for s in states),
                         'coarse_evaluations':coarse,'eigen_seeds':seeds,'retained_seeds':retained,'seconds':elapsed,
                         'undercovered_items':int(sum(error>eps)),'maximum_underestimate_ratio':float(max(error/np.maximum(eps,1e-30))),
                         'target_error':target_error,'passed':bool(target_error<=cc['cv_relative_budget'] if task=='cv' else target_error<=cc['decision_regret_budget'])})
            details.extend({'id':c['id'],'method':method,'item':i,'J_A_m2':float(y[i]),'actual_error_A_m2':float(error[i]),'estimated_error_A_m2':float(eps[i])} for i in range(len(y)))
        print('legacy',c['id'],flush=True)
    write_csv(out/'benchmark.csv',rows);write_csv(out/'sample_errors.csv',details)


def main():
    cfg=json.loads((P/'protocol.json').read_text());fresh(cfg);regression(cfg)
    files=[P/'protocol.json',P/'run_phase4.py',P/'src/selective.py']+sorted((ROOT/'src').glob('*.py'))+sorted((ROOT/'vendor').glob('*.py'))
    (P/'results/provenance.json').write_text(json.dumps({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2)+'\n')

if __name__=='__main__':main()
