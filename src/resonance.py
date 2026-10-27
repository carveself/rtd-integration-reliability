"""Empirical spectral/pole guidance, not a certified resonance enumeration.

The transport integrand is always evaluated with the frozen real-axis kernel.
Complex poles serve only as mesh hints; no Lorentzian current approximation.
"""
import time
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import root
from models import Layers
from negf import KB, Q


def make_device(case, height, spacing=.1):
    b=case['barrier_nm']; w=case['well_nm']; off=case.get('well_offset_eV',0.)
    layers=[(4.,0.),(b,height),(w,off)]
    if case['kind']=='triple':
        layers += [(b,height),(case['second_well_nm'],0.)]
    layers += [(case.get('right_barrier_nm',b),height),(4.,0.)]
    return Layers(tuple(layers),spacing_nm=spacing)


def bounds(device,bias,temp):
    vt=KB*temp/Q
    return max(0.,-bias)+1e-12,max(device.mu_eV,device.mu_eV-bias)+30*vt


def spectral_seeds(device,bias,low,high):
    """Exact real-part boundary pencil inside BOTH contact bands.

    Re Sigma=(E-edge-2t)/2. Thus A v = E B v, B_end=1/2.
    Convert to a symmetric tridiagonal ordinary eigenproblem. Widths are
    FIRST-ORDER endpoint-coupling estimates, not rigorous resonance widths.
    """
    t=device.hopping; _,u=device.profile(bias); n=len(u)
    if low<=max(0.,-bias) or high>=min(4*t,4*t-bias):
        raise ValueError('Spectral pencil requires two open contacts')
    d=2*t+u.copy(); d[0]-=t; d[-1]-=t-bias/2
    b=np.ones(n); b[[0,-1]]=.5
    val,vec=eigh_tridiagonal(d/b,-t/np.sqrt(b[:-1]*b[1:]),
                            select='v',select_range=(low,high),check_finite=True)
    v=vec/np.sqrt(b[:,None])
    gam_l=np.sqrt(4*t*t-(val-2*t)**2)
    gam_r=np.sqrt(4*t*t-(val+bias-2*t)**2)
    widths=gam_l*v[0]**2+gam_r*v[-1]**2
    return val,widths,np.argmax(abs(v),axis=0)


def sigma_complex(z,edge,t):
    """Outgoing analytic continuation across interior of the contact band.

    This is the resonance sheet (Im z<0), not the physical resolvent below
    its cut. Principal sqrt is continuous in our restricted lower strip.
    """
    x=z-edge-2*t
    return .5*x-.5j*np.sqrt(4*t*t-x*x)


def pole_residual(z,diagonal,t,bias,pivot):
    # Schur complement at an eigenvector antinode. Zeros are poles when its
    # subblocks are nonsingular; a residual test alone does not prove coverage.
    gl=1/(z-diagonal[0]-sigma_complex(z,0.,t))
    for j in range(1,pivot):gl=1/(z-diagonal[j]-t*t*gl)
    gr=1/(z-diagonal[-1]-sigma_complex(z,-bias,t))
    for j in range(len(diagonal)-2,pivot,-1):gr=1/(z-diagonal[j]-t*t*gr)
    if pivot==0:return z-diagonal[0]-sigma_complex(z,0.,t)-t*t*gr
    if pivot==len(diagonal)-1:return z-diagonal[-1]-sigma_complex(z,-bias,t)-t*t*gl
    return z-diagonal[pivot]-t*t*(gl+gr)


def discover(device,bias,temp,refine=True):
    start=time.perf_counter();lo,hi=bounds(device,bias,temp)
    tick=time.perf_counter();energies,widths,pivots=spectral_seeds(device,bias,lo,hi)
    eigen_seconds=time.perf_counter()-tick
    _,u=device.profile(bias);d=2*device.hopping+u
    trace=[];poles=[];calls=0
    for en,width,pivot in zip(energies,widths,pivots):
        row={'seed_eV':float(en),'first_order_width_eV':float(width),'pivot':int(pivot)}
        if not refine or width>.03:
            row.update(status='broad_seed' if refine else 'spectral_only',residual_calls=0)
            trace.append(row);continue
        ncall=0
        def fun(x):
            nonlocal ncall
            ncall+=1
            z=complex(x[0],x[1]);f=pole_residual(z,d,device.hopping,bias,int(pivot))
            return [f.real,f.imag]
        sol=root(fun,[en,-max(width/2,1e-10)],tol=1e-10)
        z=complex(*sol.x);err=abs(complex(*fun(sol.x)));calls+=ncall
        accepted=(lo<z.real<hi and 1e-12<-2*z.imag<.05 and err<1e-7
                  and abs(z.real-en)<.03)
        # sol.success can be false at roundoff despite a small physical residual.
        row.update(status='accepted' if accepted else 'rejected',pole_real_eV=z.real,
                   pole_width_eV=-2*z.imag,residual_eV=err,residual_calls=ncall,
                   solver_success=bool(sol.success))
        if accepted and not any(abs(z-q)<max(1e-10,-z.imag*.01) for q in poles):poles.append(z)
        trace.append(row)
    return {'energies':energies,'widths':widths,'poles':poles,'trace':trace,
            'residual_calls':calls,'eigen_seconds':eigen_seconds,
            'setup_seconds':time.perf_counter()-start}


def seed_edges(lo,hi,cap,centers,scales,offsets):
    edges=list(np.linspace(lo,hi,int(np.ceil((hi-lo)/cap))+1))
    for center,scale in zip(centers,scales):
        edges.append(float(center))
        for factor in offsets:
            edges.extend([center-scale*factor,center+scale*factor])
    return np.unique([float(e) for e in edges if lo<=e<=hi])
