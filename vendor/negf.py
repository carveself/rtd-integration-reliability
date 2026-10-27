"""Coherent open-boundary transport; eV, nm internally, SI current density."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import hbar, m_e, elementary_charge as Q, Boltzmann as KB, h
from numpy.polynomial.legendre import leggauss

KINETIC = hbar**2 / (2*m_e*Q) * 1e18  # eV nm^2

@dataclass(frozen=True)
class Device:
    spacing_nm: float = 0.1
    spacer_nm: float = 4.0
    barrier_nm: float = 3.0
    well_nm: float = 5.0
    height_eV: float = 0.28
    mass: float = 0.067
    mu_eV: float = 0.08

    def __post_init__(self):
        if not (0.02 <= self.spacing_nm <= .4 and self.mass > 0 and self.height_eV >= 0):
            raise ValueError('Invalid mesh, mass or barrier height')
        for length in (self.spacer_nm, self.barrier_nm, self.well_nm):
            if length <= 0 or abs(length/self.spacing_nm-round(length/self.spacing_nm)) > 1e-8:
                raise ValueError('Every layer must contain an integer number of cells')

    @property
    def layers(self):
        return [(self.spacer_nm,0.), (self.barrier_nm,self.height_eV),
                (self.well_nm,0.), (self.barrier_nm,self.height_eV), (self.spacer_nm,0.)]

    @property
    def hopping(self):
        return KINETIC / self.mass / self.spacing_nm**2

    def profile(self, bias=0.):
        """Cell centers; voltage drops linearly between first/last barrier faces."""
        counts=[round(d/self.spacing_nm) for d,_ in self.layers]
        x=(np.arange(sum(counts))+.5)*self.spacing_nm
        offsets=np.repeat([u for _,u in self.layers],counts)
        active=2*self.barrier_nm+self.well_nm
        ramp=np.clip((x-self.spacer_nm)/active,0.,1.)
        return x, offsets-bias*ramp


def self_energy(energy, edge, hopping):
    """Exact retarded real-axis self-energy for a semi-infinite matched chain.

    Stable small quadratic root outside the band. No artificial imaginary eta.
    """
    x=np.asarray(energy,dtype=float)-edge-2*hopping
    inside=np.abs(x) < 2*hopping
    root=np.sqrt(np.maximum(x*x-4*hopping*hopping,0.))
    outside=2*hopping*hopping/np.where(inside,1.,x+np.where(x>=0,1.,-1.)*root)
    return np.where(inside, .5*x-.5j*np.sqrt(np.maximum(4*hopping*hopping-x*x,0.)), outside+0j)


def transport(device, energy, bias=0., shift=0., local=False):
    """Vectorized recursive Green function; return T/R and optional injectivities.

    shift moves all energies/edges and tests gauge invariance. E is the *total*
    longitudinal energy in the same reference as both lead band edges.
    """
    e=np.atleast_1d(np.asarray(energy,float))
    if e.ndim != 1 or not np.all(np.isfinite(e)) or not np.isfinite(bias+shift):
        raise ValueError('Finite one-dimensional energies and finite bias/shift are required')
    _,u=device.profile(bias); u=u+shift
    t=device.hopping
    sl=self_energy(e,shift,t); sr=self_energy(e,shift-bias,t)
    gl=-2*sl.imag; gr=-2*sr.imag
    a=e[:,None]-(2*t+u)[None,:]+0j
    a[:,0]-=sl; a[:,-1]-=sr
    # Band-edge singular points are outside our quadrature (interior Gauss nodes).
    if np.any((gl==0)&(gr==0)):
        raise ValueError('Choose energies inside at least one contact band; avoid exact band edges')
    right=np.empty_like(a)
    right[:,-1]=1/a[:,-1]
    for j in range(len(u)-2,-1,-1):
        right[:,j]=1/(a[:,j]-t*t*right[:,j+1])
    c0=right[:,0]
    end=c0.copy()
    cols_l=np.empty_like(a) if local else None
    if local: cols_l[:,0]=end
    for j in range(1,len(u)):
        end=-t*right[:,j]*end
        if local: cols_l[:,j]=end
    transmission=gl*gr*np.abs(end)**2
    reflection=np.abs(1-1j*gl*c0)**2
    result={'T':transmission,'R':reflection,'Gamma_L':gl,'Gamma_R':gr}
    if local:
        left=np.empty_like(a); left[:,0]=1/a[:,0]
        for j in range(1,len(u)):
            left[:,j]=1/(a[:,j]-t*t*left[:,j-1])
        cols_r=np.empty_like(a); cols_r[:,-1]=left[:,-1]
        for j in range(len(u)-2,-1,-1):
            cols_r[:,j]=-t*left[:,j]*cols_r[:,j+1]
        # Contact-resolved spectral functions, per site, per spin, in 1/eV.
        al=gl[:,None]*abs(cols_l)**2; ar=gr[:,None]*abs(cols_r)**2
        diag=1/(a-np.column_stack([np.zeros(len(e)),t*t*left[:,:-1]])
                  -np.column_stack([t*t*right[:,1:],np.zeros(len(e))]))
        result.update(A_L=al,A_R=ar,LDOS=(al+ar)/(2*np.pi*device.spacing_nm),
                      spectral_identity_error=np.max(abs(al+ar+2*diag.imag)),
                      bond_left=2*t*gl[:,None]*np.imag(cols_l[:,:-1].conj()*cols_l[:,1:]),
                      bond_right=2*t*gr[:,None]*np.imag(cols_r[:,:-1].conj()*cols_r[:,1:]))
    return result


def supply(energy,bias,temperature,mu=.08):
    """Integrated 2D transverse Fermi supply in eV; spin included in J prefactor."""
    if temperature <= 0: raise ValueError('Temperature must be positive')
    vt=KB*temperature/Q
    e=np.asarray(energy)
    return vt*(np.logaddexp(0.,(mu-e)/vt)-np.logaddexp(0.,(mu-bias-e)/vt))


def current_prefactor(mass):
    # J = q*m/(2*pi^2*hbar^3) integral dE_J T * kBT_J log(...).
    # Both integral factors here use eV, hence q^3; includes spin 2, one valley.
    return Q**3*(mass*m_e)/(2*np.pi**2*hbar**3)


def integrate_current(device,bias,temperature=300.,rtol=2e-7,panel_eV=.004,tail_kT=30.):
    """Adaptive paired 8/16 Gauss panels; finite tail bound is returned explicitly.

    The initial panel cap is essential: nested agreement alone can miss a narrow
    peak. Convergence experiments halve it and compare a separate quadrature.
    """
    if not (temperature>0 and rtol>0 and panel_eV>0 and tail_kT>0):
        raise ValueError('Positive temperature and quadrature controls are required')
    vt=KB*temperature/Q
    low=max(0.,-bias)+1e-12
    high=max(device.mu_eV,device.mu_eV-bias)+tail_kT*vt
    if high>=min(4*device.hopping,4*device.hopping-bias):
        raise ValueError('Energy cutoff reaches the discretized contact band ceiling')
    edges=np.linspace(low,high,int(np.ceil((high-low)/panel_eV))+1)
    l=edges[:-1]; r=edges[1:]
    x8,w8=leggauss(8); x16,w16=leggauss(16)
    total=0.; error=0.; evaluations=0
    for depth in range(24):
        mid=(l+r)/2; half=(r-l)/2
        values=[]
        for x,w in ((x8,w8),(x16,w16)):
            en=mid[:,None]+half[:,None]*x
            trans=transport(device,en.ravel(),bias)['T'].reshape(en.shape)
            y=trans*supply(en,bias,temperature,device.mu_eV)
            values.append(half*(y@w)); evaluations+=en.size
        err=np.abs(values[1]-values[0])
        tol=rtol*np.abs(values[1])+1e-13*(r-l)/(high-low)
        if error+np.sum(err) <= rtol*abs(total+np.sum(values[1]))+1e-13:
            total+=np.sum(values[1]); error+=np.sum(err); break
        accept=err<=tol
        total+=np.sum(values[1][accept]); error+=np.sum(err[accept])
        if np.all(accept): break
        ll=l[~accept]; rr=r[~accept]; mm=(ll+rr)/2
        l=np.concatenate([ll,mm]); r=np.concatenate([mm,rr])
    else: raise RuntimeError('Energy quadrature failed to converge')
    # Since T<=1, log(1+exp(x))<=exp(x), each reservoir tail <=VT^2 exp((mu-Emax)/VT).
    tail=vt*vt*(np.exp((device.mu_eV-high)/vt)+np.exp((device.mu_eV-bias-high)/vt))
    pref=current_prefactor(device.mass)
    return {'J_A_m2':pref*total,'estimated_error_A_m2':pref*error,
            'tail_bound_A_m2':pref*tail,'evaluations':evaluations,'Emax_eV':high}
