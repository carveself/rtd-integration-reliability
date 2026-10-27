"""Shared paired quadrature; setup is counted in total timings by the runner."""
import numpy as np
from integration import PanelIntegral
from resonance import bounds,seed_edges
from negf import current_prefactor,KB,Q

class SeededIntegral(PanelIntegral):
    def __init__(self,device,bias,temp,cfg,discovery=None,mode='uniform',cap=None):
        self.device=device;self.bias=bias;self.temp=temp;self.cfg=cfg
        self.evaluations=0;self.refinements=0;self.pref=current_prefactor(device.mass)
        vt=KB*temp/Q;lo,hi=bounds(device,bias,temp)
        self.tail=self.pref*vt**2*(np.exp((device.mu_eV-hi)/vt)+np.exp((device.mu_eV-bias-hi)/vt))
        centers=[];scales=[]
        if mode=='spectral':
            centers=list(discovery['energies']);scales=[cfg['fixed_seed_scale_eV']]*len(centers)
        elif mode=='pole':
            centers=[p.real for p in discovery['poles']]
            scales=[-2*p.imag for p in discovery['poles']]
            # Retain EVERY spectral seed as a fallback for broad/rejected roots.
            centers+=list(discovery['energies']);scales += [cfg['fixed_seed_scale_eV']]*len(discovery['energies'])
        elif mode!='uniform':raise ValueError(mode)
        edges=seed_edges(lo,hi,cap or cfg['initial_panel_eV'],centers,scales,cfg['pole_offsets_in_widths'])
        self.l=edges[:-1];self.r=edges[1:];self.v,self.e=self.panels(self.l,self.r)
        self.initial_evaluations=self.evaluations
