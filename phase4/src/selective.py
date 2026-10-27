"""Supply-aware selective spectral seeding; an empirical risk heuristic.

Every spectral solve and coarse transport point is counted. No reference values
are consulted. Estimated resonance areas are NOT mathematical upper bounds.
"""
from pathlib import Path
import sys,time
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'vendor')]
from resonance import discover,bounds
from seeded import SeededIntegral
from negf import transport,supply,current_prefactor


def choose_seeds(energies,widths,areas,budget,width_limit):
    """Drop the smallest proxy areas while total dropped proxy <= budget.

    Broad modes rely on ordinary adaptive quadrature. This empirical width
    filter is part of the hypothesis and is explicitly not a safety guarantee.
    """
    energies=np.asarray(energies);widths=np.asarray(widths);areas=np.asarray(areas)
    if energies.shape!=widths.shape or energies.shape!=areas.shape:raise ValueError('Shape mismatch')
    if budget<0 or np.any(areas<0) or not np.all(np.isfinite(areas)):raise ValueError('Invalid area proxy')
    narrow=np.flatnonzero((widths>=0)&(widths<width_limit));selected=set(narrow.tolist());omitted=0.
    for idx in narrow[np.argsort(areas[narrow],kind='stable')]:
        if omitted+areas[idx]<=budget:
            selected.remove(int(idx));omitted+=float(areas[idx])
        else:break
    return np.array(sorted(selected),int),omitted


class SelectiveIntegral(SeededIntegral):
    def __init__(self,device,bias,temp,cfg):
        start=time.perf_counter();lo,hi=bounds(device,bias,temp);control=cfg['selection_rule']
        nodes,weights=leggauss(control['coarse_order'])
        edges=np.linspace(lo,hi,int(np.ceil((hi-lo)/cfg['initial_panel_eV']))+1)
        mid=(edges[:-1]+edges[1:])/2;half=np.diff(edges)/2
        en=(mid[:,None]+half[:,None]*nodes).ravel();vals=[]
        for j in range(0,len(en),cfg['batch_limit']):
            e=en[j:j+cfg['batch_limit']]
            vals.append(transport(device,e,bias)['T']*supply(e,bias,temp,device.mu_eV))
        pref=current_prefactor(device.mass)
        coarse=float(pref*np.sum(half*(np.concatenate(vals).reshape(-1,len(nodes))@weights)))
        discovery=discover(device,bias,temp,refine=False)
        proxy=pref*np.pi/2*discovery['widths']*abs(supply(discovery['energies'],bias,temp,device.mu_eV))
        budget=control['risk_fraction']*(cfg['rtol']*abs(coarse)+cfg['absolute_floor_A_m2'])
        ids,omitted=choose_seeds(discovery['energies'],discovery['widths'],proxy,budget,control['width_threshold_eV'])
        reduced={**discovery,'energies':discovery['energies'][ids]}
        setup=time.perf_counter()-start
        super().__init__(device,bias,temp,cfg,reduced,'spectral')
        self.evaluations+=len(en);self.initial_evaluations+=len(en)
        self.selection={'coarse_J_A_m2':coarse,'coarse_evaluations':len(en),
                        'eigen_seeds':len(discovery['energies']),'retained_seeds':len(ids),
                        'omitted_narrow_area_proxy_A_m2':omitted,'proxy_budget_A_m2':budget,
                        'eigen_seconds':discovery['eigen_seconds'],'setup_seconds':setup,
                        'trace':[{'energy_eV':float(e),'first_order_width_eV':float(w),
                                  'area_proxy_A_m2':float(a),'retained':bool(i in ids)}
                                 for i,(e,w,a) in enumerate(zip(discovery['energies'],discovery['widths'],proxy))]}
