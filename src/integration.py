"""Cached global Gauss panel refinement shared by both candidate and baseline.

The paired-rule error is an empirical estimate, NOT a rigorous upper bound.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from models import Device
from negf import transport,supply,current_prefactor,KB,Q

RULES=[leggauss(8),leggauss(16)]

class PanelIntegral:
    def __init__(self,device,bias,temp,cfg):
        self.device=device;self.bias=bias;self.temp=temp;self.cfg=cfg;self.evaluations=0;self.refinements=0
        self.pref=current_prefactor(device.mass)
        vt=KB*temp/Q;lo=max(0.,-bias)+1e-12;hi=max(device.mu_eV,device.mu_eV-bias)+30*vt
        if hi>=min(4*device.hopping,4*device.hopping-bias):raise ValueError('Contact band ceiling')
        self.tail=self.pref*vt**2*(np.exp((device.mu_eV-hi)/vt)+np.exp((device.mu_eV-bias-hi)/vt))
        edges=np.linspace(lo,hi,int(np.ceil((hi-lo)/cfg['initial_panel_eV']))+1)
        self.l=edges[:-1];self.r=edges[1:];self.v,self.e=self.panels(self.l,self.r)
    def values(self,energies):
        flat=energies.ravel();out=[]
        for start in range(0,len(flat),self.cfg['batch_limit']):
            en=flat[start:start+self.cfg['batch_limit']]
            out.append(transport(self.device,en,self.bias)['T']*supply(en,self.bias,self.temp,self.device.mu_eV))
        self.evaluations+=len(flat)
        return np.concatenate(out).reshape(energies.shape)
    def panels(self,l,r):
        mid=(l+r)/2;half=(r-l)/2;ints=[]
        for nodes,weights in RULES:
            en=mid[:,None]+half[:,None]*nodes
            ints.append(self.pref*half*(self.values(en)@weights))
        return ints[1],np.abs(ints[1]-ints[0])
    @property
    def value(self):return float(np.sum(self.v))
    @property
    def error(self):return float(self.cfg['error_safety']*np.sum(self.e)+self.tail)
    def refine(self):
        n=min(self.cfg['refinement_batch_panels'],len(self.l))
        ids=np.argsort(self.e)[-n:];keep=np.ones(len(self.l),bool);keep[ids]=False
        l=self.l[ids];r=self.r[ids];mid=(l+r)/2
        nl=np.r_[l,mid];nr=np.r_[mid,r];v,e=self.panels(nl,nr)
        self.l=np.r_[self.l[keep],nl];self.r=np.r_[self.r[keep],nr]
        self.v=np.r_[self.v[keep],v];self.e=np.r_[self.e[keep],e];self.refinements+=1
    def finish(self,rtol):
        while self.error>rtol*abs(self.value)+1e-6:
            if self.refinements>=2000:raise RuntimeError('Refinement limit')
            self.refine()
        return self
