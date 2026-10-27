"""Separate discretized solve and continuum scattering composition references."""
import numpy as np
from scipy.linalg import solve_banded
from negf import KINETIC, self_energy


def direct_banded(device,energy,bias=0.):
    """Same discrete problem, different solver: two full inverse columns."""
    t=device.hopping; _,u=device.profile(bias); n=len(u)
    sl=self_energy(energy,0.,t); sr=self_energy(energy,-bias,t)
    ab=np.zeros((3,n),complex); ab[0,1:]=t; ab[2,:-1]=t
    ab[1]=energy-2*t-u; ab[1,0]-=sl; ab[1,-1]-=sr
    rhs=np.zeros((n,2),complex); rhs[0,0]=1.; rhs[-1,1]=1.
    g=solve_banded((1,1),ab,rhs)
    return float((-2*sl.imag)*(-2*sr.imag)*abs(g[-1,0])**2)


def cascade(a,b):
    """Redheffer star product in amplitude normalization, including evanescence."""
    r1,tr1,tl1,rr1=a; r2,tr2,tl2,rr2=b
    den=1-rr1*r2
    return (r1+tr1*r2*tl1/den,tr1*tr2/den,tl2*tl1/den,rr2+tl2*rr1*tr2/den)


def continuum(energy,layers,mass=.067,left_edge=0.,right_edge=0.):
    """Piecewise constant continuum Schrodinger scattering (constant mass).

    Each layer: (length_nm, potential_eV). Avoid exact internal turning points.
    No discrete-chain self-energy or recursive Green matrix is used.
    """
    if energy <= max(left_edge,right_edge): return 0.
    potentials=[left_edge]+[v for _,v in layers]+[right_edge]
    k=np.sqrt((energy-np.array(potentials)+0j)*mass/KINETIC)
    if np.any(abs(k)<1e-12): raise ValueError('Avoid exact turning-point energy')
    s=(0j,1+0j,1+0j,0j)
    for j in range(len(k)-1):
        den=k[j]+k[j+1]
        interface=((k[j]-k[j+1])/den,2*k[j+1]/den,2*k[j]/den,(k[j+1]-k[j])/den)
        s=cascade(s,interface)
        if j<len(layers):
            p=np.exp(1j*k[j+1]*layers[j][0]); s=cascade(s,(0j,p,p,0j))
    return float((k[-1].real/k[0].real)*abs(s[2])**2)


def rectangular_exact(energy,height,width,mass=.067):
    if energy <= 0: return 0.
    if abs(energy-height)<1e-12:
        return 1/(1+mass*height*width*width/(4*KINETIC))
    if energy<height:
        z=np.sqrt(mass*(height-energy)/KINETIC)*width
        return 1/(1+height**2*np.sinh(z)**2/(4*energy*(height-energy)))
    z=np.sqrt(mass*(energy-height)/KINETIC)*width
    return 1/(1+height**2*np.sin(z)**2/(4*energy*(energy-height)))
