import unittest,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'vendor')]
from resonance import spectral_seeds,sigma_complex,discover,seed_edges,bounds,make_device
from seeded import SeededIntegral
from integration import PanelIntegral
from models import Device
from negf import self_energy,transport
from reference import direct_banded
from config import config
from highprecision import MPChain,integrate_window
import json
CFG=config(json.loads((ROOT/'protocol.json').read_text()))

class Toy(PanelIntegral):
    """An analytically integrable narrow spike with exact known center/width."""
    def __init__(self,center,gamma,guided):
        self.center=center;self.gamma=gamma;self.cfg=CFG;self.pref=1.;self.tail=0.
        self.evaluations=0;self.refinements=0
        edges=seed_edges(0.,.1,.008,[center] if guided else [],[2*gamma] if guided else [],CFG['pole_offsets_in_widths'])
        self.l=edges[:-1];self.r=edges[1:];self.v,self.e=self.panels(self.l,self.r)
    def values(self,en):
        self.evaluations+=en.size
        return 1e9/(1+((en-self.center)/self.gamma)**2)

class Physics(unittest.TestCase):
    def test_continuation_matches_retarded_real_axis(self):
        t=56.;e=np.array([.001,.01,.1,1.])
        np.testing.assert_allclose(sigma_complex(e+0j,-.07,t),self_energy(e,-.07,t),rtol=1e-12)
    def test_boundary_pencil_independent_residual(self):
        dev=Device(barrier_nm=7.,well_nm=5.4);bias=.07;lo,hi=bounds(dev,bias,150.)
        vals,_,_=spectral_seeds(dev,bias,lo,hi)
        # A true pencil seed makes the REAL open-boundary Hamiltonian singular.
        _,u=dev.profile(bias);t=dev.hopping
        for e in vals:
            matrix=np.diag(2*t+u)+np.diag(np.full(len(u)-1,-t),1)+np.diag(np.full(len(u)-1,-t),-1)
            matrix[0,0]+=self_energy(e,0.,t).real
            matrix[-1,-1]+=self_energy(e,-bias,t).real
            self.assertLess(np.min(abs(np.linalg.eigvalsh(matrix)-e)),1e-9)
    def test_pole_center_and_linewidth_from_separate_solver(self):
        dev=Device(barrier_nm=8.,well_nm=5.6)
        d=discover(dev,0.,150.);self.assertGreaterEqual(len(d['poles']),1)
        z=min(d['poles'],key=lambda z:z.real);half=-z.imag
        t0=direct_banded(dev,z.real,0.)
        self.assertGreater(t0,.99)
        for sign in (-1,1):self.assertAlmostEqual(direct_banded(dev,z.real+sign*half,0.)/t0,.5,delta=.01)
    def test_discrete_solvers_match_narrow_resonance(self):
        dev=Device(barrier_nm=8.,well_nm=5.6);d=discover(dev,.04,100.)
        z=min(d['poles'],key=lambda z:z.real);en=z.real+np.linspace(-3,3,7)*(-z.imag)
        np.testing.assert_allclose(transport(dev,en,.04)['T'],[direct_banded(dev,e,.04) for e in en],rtol=2e-8)
    def test_no_resonance_guidance_free_chain(self):
        dev=Device(height_eV=0.);bias=.05
        # At bias zero uniform chain T=1; avoid evaluating J=0 as a relative target.
        np.testing.assert_allclose(transport(dev,np.array([.01,.1,.3]),0.)['T'],1.,atol=2e-12)
    def test_extremely_narrow_spike_exposes_blind_indicator(self):
        c=.037173;g=1e-10;exact=1e9*g*(np.arctan((.1-c)/g)+np.arctan(c/g))
        blind=Toy(c,g,False).finish(.001)
        self.assertLess(blind.error,.001*abs(blind.value)+1e-6)
        self.assertGreater(abs(blind.value-exact)/exact,.99)
        guided=Toy(c,g,True).finish(.001)
        self.assertLess(abs(guided.value-exact)/exact,.001)
    def test_partition_covers_once_and_excludes_outside(self):
        edges=seed_edges(0.,.2,.008,[.000001,.199999],[1e-7,1e-7],[.25,1,64,1024])
        self.assertEqual(edges[0],0.);self.assertEqual(edges[-1],.2)
        self.assertTrue(np.all(np.diff(edges)>0))
        self.assertAlmostEqual(np.sum(np.diff(edges)),.2)
    def test_reject_spectral_contact_cut(self):
        with self.assertRaises(ValueError):spectral_seeds(Device(),0.,0.,.1)
    def test_layer_geometry_and_asymmetric_offset(self):
        case={'kind':'double','barrier_nm':6.5,'right_barrier_nm':8.,'well_nm':5.,'well_offset_eV':.012}
        dev=make_device(case,.28);_,u=dev.profile(0.)
        self.assertEqual(len(u),275);self.assertEqual(np.count_nonzero(u==.012),50)
    def test_arbitrary_precision_free_chain_conservation(self):
        chain=MPChain(Device(height_eV=0.),0.,100.,50)
        self.assertLess(abs(chain.transmission('.05')-1),chain.ctx.mpf('1e-40'))
    def test_mp_local_quadrature_against_analytic_spike(self):
        chain=MPChain(Device(),0.,100.,50);m=chain.ctx
        c=.04;g=1e-11;left=c-128*g;right=c+128*g
        def integrand(e):
            chain.calls+=1
            return m.mpf('1e9')/(1+((e-m.mpf(c))/m.mpf(g))**2)
        chain.current_integrand=integrand
        val,_,_=integrate_window(chain,left,right,c,2*g,(10,20),1e-10)
        exact=1e9*m.mpf(g)*(m.atan((m.mpf(right)-m.mpf(c))/m.mpf(g))-m.atan((m.mpf(left)-m.mpf(c))/m.mpf(g)))
        self.assertLess(abs(val-float(exact))/float(exact),1e-9)

if __name__=='__main__':unittest.main(verbosity=2)
