"""Constant-mass layer models compatible with the frozen NEGF kernel."""
from pathlib import Path
import sys
from dataclasses import dataclass
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'vendor'))
from negf import Device, KINETIC

@dataclass(frozen=True)
class Layers:
    layers: tuple
    spacing_nm: float=.1
    mass: float=.067
    mu_eV: float=.08
    spacer_nm: float=4.

    def __post_init__(self):
        for width,edge in self.layers:
            if width<=0 or abs(width/self.spacing_nm-round(width/self.spacing_nm))>1e-8:
                raise ValueError('Layers must contain integer cells')
        if self.layers[0][1]!=0 or self.layers[-1][1]!=0:
            raise ValueError('Lead-like spacers must have zero offset')
    @property
    def hopping(self):return KINETIC/self.mass/self.spacing_nm**2
    def profile(self,bias=0.):
        widths=[w for w,_ in self.layers]
        counts=[round(w/self.spacing_nm) for w in widths]
        x=(np.arange(sum(counts))+.5)*self.spacing_nm
        u=np.repeat([e for _,e in self.layers],counts)
        active=sum(widths)-widths[0]-widths[-1]
        return x,u-bias*np.clip((x-widths[0])/active,0,1)


def model(case,height=.28,well=None,spacing=.1):
    w=case.get('well_nm',5.) if well is None else well
    b=case['barrier_nm']
    if case['kind']=='double':
        return Device(barrier_nm=b,well_nm=w,height_eV=height,spacing_nm=spacing)
    if case['kind']!='triple':raise ValueError('Unknown structure')
    layers=((4.,0.),(b,height),(w,0.),(b,height),(case['second_well_nm'],0.),(b,height),(4.,0.))
    return Layers(layers=layers,spacing_nm=spacing)
