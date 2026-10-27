from pathlib import Path
import sys,unittest,json
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'phase4/src'),str(ROOT/'src'),str(ROOT/'vendor')]
from selective import choose_seeds,SelectiveIntegral
from seeded import SeededIntegral
from models import Device
CFG=json.loads((ROOT/'phase4/protocol.json').read_text())

class SelectionTests(unittest.TestCase):
    def test_cumulative_discard_budget(self):
        ids,omitted=choose_seeds(np.arange(4),np.ones(4)*1e-4,np.array([.04,.06,1.,2.]),.11,.008)
        np.testing.assert_array_equal(ids,[2,3]);self.assertAlmostEqual(omitted,.10)
    def test_higher_budget_never_retains_more(self):
        args=(np.arange(5),np.ones(5)*1e-4,np.array([.01,.03,.2,1.,10.]))
        a,_=choose_seeds(*args,.04,.008);b,_=choose_seeds(*args,.25,.008)
        self.assertTrue(set(b)<=set(a))
    def test_no_narrow_modes(self):
        ids,area=choose_seeds(np.array([.1]),np.array([.1]),np.array([1.]),0.,.008)
        self.assertEqual(len(ids),0);self.assertEqual(area,0.)
    def test_invalid_area_rejected(self):
        with self.assertRaises(ValueError):choose_seeds(np.array([0]),np.array([0]),np.array([np.nan]),1.,.008)
    def test_counts_all_discarded_coarse_work(self):
        dev=Device(height_eV=0.);s=SelectiveIntegral(dev,.05,150.,CFG)
        base=SeededIntegral(dev,.05,150.,CFG,{'energies':np.array([])},'spectral')
        self.assertEqual(s.selection['retained_seeds'],0)
        self.assertEqual(s.evaluations,base.evaluations+s.selection['coarse_evaluations'])
    def test_physical_current_agrees_with_strict_reference(self):
        dev=Device(barrier_nm=7.,well_nm=5.4)
        candidate=SelectiveIntegral(dev,.07,150.,CFG).finish(.001)
        reference=SeededIntegral(dev,.07,150.,CFG,mode='uniform',cap=.00025).finish(1e-8)
        self.assertLess(abs(candidate.value/reference.value-1),.001)

if __name__=='__main__':unittest.main(verbosity=2)
