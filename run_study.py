"""Run the frozen 24-case study. References are computed after candidate stops."""
from pathlib import Path
import sys,json,csv,time,hashlib,platform,warnings,argparse
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'vendor')]
from scipy.integrate import quad,IntegrationWarning
import scipy
from resonance import make_device,discover,bounds,seed_edges
from seeded import SeededIntegral
from config import config
from reference import direct_banded
from negf import supply,current_prefactor


def write_csv(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader();wr.writerows(rows)


def run_method(device,case,cfg,method):
    start=time.perf_counter();d=None
    if method!='uniform':d=discover(device,case['bias'],case['temperature'],refine=method=='pole')
    st=SeededIntegral(device,case['bias'],case['temperature'],cfg,d,method).finish(cfg['rtol'])
    elapsed=time.perf_counter()-start
    return st,d,elapsed


def independent_reference(device,case,cfg,discovery,reference_scale):
    """QUADPACK + independent banded Green solve; guidance is SHARED.

    Agreement with two unseeded fine grids is required. Shared guidance is
    disclosed: this is an empirical reference audit, never a no-missed-pole proof.
    """
    lo,hi=bounds(device,case['bias'],case['temperature'])
    poles=discovery['poles']
    edges=seed_edges(lo,hi,.004,[z.real for z in poles],[-2*z.imag for z in poles],[.25,1,4,16,64])
    pref=current_prefactor(device.mass);calls=0
    def fun(en):
        nonlocal calls
        calls+=1
        return pref*direct_banded(device,en,case['bias'])*float(supply(en,case['bias'],case['temperature'],device.mu_eV))
    value=0.;error=0.;messages=[]
    # Absolute allocation is deliberately looser than machine precision but
    # much tighter than the candidate's 1e-3 relative target.
    epsabs=max(1e-8,abs(reference_scale)*1e-11)/len(edges)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always',IntegrationWarning)
        for l,r in zip(edges[:-1],edges[1:]):
            val,err=quad(fun,l,r,epsabs=epsabs,epsrel=1e-9,limit=200)
            value+=val;error+=err
        messages=[str(w.message) for w in caught]
    return value,error,calls,messages


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--protocol',default='protocol.json')
    parser.add_argument('--results',default='results')
    args=parser.parse_args()
    protocol_path=ROOT/args.protocol
    cfg=config(json.loads(protocol_path.read_text()));out=ROOT/args.results;out.mkdir(exist_ok=True)
    rng=np.random.default_rng(cfg['seed']);cases=[]
    for context in cfg['contexts']:
        for i,height in enumerate(rng.normal(.28,cfg['height_sigma_eV'],cfg['draws_per_context'])):
            cases.append({**context,'sample':i,'height_eV':float(height)})
    (out/'cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    rows=[];refs=[];traces={};timings=[]
    for k,case in enumerate(cases):
        device=make_device(case,case['height_eV'],cfg['spacing_nm']);n=len(device.profile(case['bias'])[0])
        methods=['uniform','spectral','pole'];runs={m:[] for m in methods};saved={}
        for rep in range(cfg['timing_repetitions']):
            # Rotate order to reduce a fixed-order timing advantage. All setup
            # work is repeated; no free discovery cache between methods.
            for method in methods[rep:]+methods[:rep]:
                st,discovery,elapsed=run_method(device,case,cfg,method)
                runs[method].append(elapsed)
                timings.append({'id':case['id'],'sample':case['sample'],'method':method,'repeat':rep,'seconds':elapsed,
                                'setup_seconds':0 if discovery is None else discovery['setup_seconds']})
                if rep==0:saved[method]=(st,discovery)
        # Neither true current nor reference error is available to candidate.
        fine=[];fine_failures=[]
        for cap in cfg['reference_panel_eV']:
            st=SeededIntegral(device,case['bias'],case['temperature'],cfg,mode='uniform',cap=cap)
            try:
                st.finish(cfg['reference_rtol'])
                fine_failures.append('')
            except RuntimeError as exc:
                # Keep the incomplete iterate for diagnosis; never label it
                # an accepted reference or abort remaining adversarial cases.
                fine_failures.append(str(exc))
            fine.append(st)
        pole_discovery=saved['pole'][1]
        q,qe,qcalls,messages=independent_reference(device,case,cfg,pole_discovery,fine[-1].value)
        agreement=max(abs(fine[0].value-q),abs(fine[1].value-q))/max(abs(q),1e-30)
        accepted=agreement<=cfg['reference_accept_relative_agreement'] and len(messages)==0 and not any(fine_failures)
        refs.append({'id':case['id'],'sample':case['sample'],'height_eV':case['height_eV'],
                     'reference_A_m2':q,'quad_error_A_m2':qe,'banded_calls':qcalls,
                     'fine1_A_m2':fine[0].value,'fine2_A_m2':fine[1].value,
                     'fine1_evaluations':fine[0].evaluations,'fine2_evaluations':fine[1].evaluations,
                     'maximum_agreement_relative':agreement,'reference_accepted':accepted,
                     'quad_warnings':json.dumps(messages),'fine_grid_failures':json.dumps(fine_failures)})
        key=f"{case['id']}/{case['sample']}";traces[key]=pole_discovery['trace']
        for method,(st,discovery) in saved.items():
            abs_error=abs(st.value-q);rel=abs_error/max(abs(q),1e-30)
            rows.append({'id':case['id'],'sample':case['sample'],'height_eV':case['height_eV'],'sites':n,
                         'method':method,'J_A_m2':st.value,'estimated_error_A_m2':st.error,
                         'absolute_error_A_m2':abs_error,'relative_error':rel,
                         'underestimate_ratio':abs_error/max(st.error,1e-300),
                         'target_pass':rel<=cfg['rtol']+cfg['absolute_floor_A_m2']/max(abs(q),1e-30),
                         'indicator_covers':abs_error<=st.error,'reference_accepted':accepted,
                         'energy_evaluations':st.evaluations,'initial_evaluations':st.initial_evaluations,
                         'pole_residual_calls':0 if discovery is None else discovery['residual_calls'],
                         'pole_recursion_site_visits':0 if discovery is None else n*discovery['residual_calls'],
                         'eigen_seconds':0 if discovery is None else discovery['eigen_seconds'],
                         'setup_seconds':0 if discovery is None else discovery['setup_seconds'],
                         'poles':0 if discovery is None else len(discovery['poles']),
                         'seconds_median':float(np.median(runs[method])),
                         'seconds_min':min(runs[method]),'seconds_max':max(runs[method])})
        print(key,'reference agreement',f'{agreement:.3g}','accepted',accepted,flush=True)
        write_csv(out/'benchmark.csv',rows);write_csv(out/'references.csv',refs)
    write_csv(out/'timings.csv',timings)
    (out/'pole_trace.json').write_text(json.dumps(traces,indent=2)+'\n')
    provenance={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),
                'threads':'runner launched with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1',
                'protocol_sha256':hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
                'files':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for folder in ('src','vendor') for p in sorted((ROOT/folder).glob('*.py'))}}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')

if __name__=='__main__':main()
