"""Post-outcome reference repair only; candidate methods are unchanged."""
from pathlib import Path
import sys,json,csv,time,warnings
import numpy as np
ROOT=Path(__file__).resolve().parent;sys.path[:0]=[str(ROOT/'src'),str(ROOT/'vendor')]
from scipy.integrate import quad,IntegrationWarning
from resonance import make_device,discover,bounds
from highprecision import MPChain,integrate_window
from reference import direct_banded
from negf import current_prefactor,supply
from run_study import write_csv


def main():
    folder=ROOT/'results_stress';cases=json.loads((folder/'cases.json').read_text())
    refs=list(csv.DictReader((folder/'references.csv').open()));rows=[]
    for case,ref in zip(cases,refs):
        if ref['reference_accepted']=='True':continue
        start=time.perf_counter();dev=make_device(case,case['height_eV']);bias=case['bias'];temp=case['temperature']
        disc=discover(dev,bias,temp);lo,hi=bounds(dev,bias,temp)
        # Only narrow poles need MP; broad states remain in the ordinary banded
        # quadrature. Windows must be disjoint to prevent double integration.
        poles=sorted([z for z in disc['poles'] if -2*z.imag<1e-5],key=lambda z:z.real)
        windows=[(max(lo,z.real-64*(-2*z.imag)),min(hi,z.real+64*(-2*z.imag)),z) for z in poles]
        if any(a[1]>=b[0] for a,b in zip(windows[:-1],windows[1:])):
            raise RuntimeError('Overlapping MP windows require a merged-window implementation')
        mpa=MPChain(dev,bias,temp,50);mpb=MPChain(dev,bias,temp,80)
        local_a=0.;local_b=0.;local_err=0.
        for l,r,z in windows:
            a,ea,_=integrate_window(mpa,l,r,z.real,-2*z.imag,(10,20),1e-9)
            b,eb,_=integrate_window(mpb,l,r,z.real,-2*z.imag,(16,32),1e-10)
            local_a+=a;local_b+=b;local_err+=eb
        edges=np.unique(np.r_[np.linspace(lo,hi,int(np.ceil((hi-lo)/.004))+1),[v for l,r,_ in windows for v in (l,r)]])
        calls=0;pref=current_prefactor(dev.mass)
        def fun(en):
            nonlocal calls
            calls+=1
            return pref*direct_banded(dev,en,bias)*float(supply(en,bias,temp,dev.mu_eV))
        outside=0.;outside_err=0.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always',IntegrationWarning)
            for l,r in zip(edges[:-1],edges[1:]):
                mid=(l+r)/2
                if any(a<mid<b for a,b,_ in windows):continue
                a,e=quad(fun,l,r,epsabs=max(1e-10,abs(float(ref['reference_A_m2']))*1e-11)/len(edges),epsrel=1e-8,limit=200)
                outside+=a;outside_err+=e
        a=outside+local_a;b=outside+local_b;disagreement=abs(a-b)/max(abs(b),1e-30)
        accepted=disagreement<1e-6 and len(caught)==0 and (outside_err+local_err)<1e-6*abs(b)
        rows.append({'id':case['id'],'sample':case['sample'],'reference_MP_A_m2':b,'mp_windows':len(windows),
                     'MP50_vs_MP80_relative':disagreement,'reported_error_A_m2':outside_err+local_err,
                     'outside_banded_calls':calls,'MP50_calls':mpa.calls,'MP80_calls':mpb.calls,
                     'outside_warning_count':len(caught),'accepted_posthoc_audit':accepted,
                     'original_double_reference_relative_difference':abs(float(ref['reference_A_m2'])-b)/max(abs(b),1e-30),
                     'seconds':time.perf_counter()-start})
        write_csv(folder/'precision_audit.csv',rows)
        print(case['id'],case['sample'],'windows',len(windows),'MP agree',disagreement,'accepted',accepted,'seconds',rows[-1]['seconds'],flush=True)

if __name__=='__main__':main()
