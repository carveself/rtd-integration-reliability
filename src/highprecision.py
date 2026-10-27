"""Independent arbitrary-precision continuant determinant for local peak audits.

Constant-mass chain only. High precision resolves floating-point cancellation;
this does not certify discovery of all resonances or model accuracy.
"""
import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from negf import KB,Q,current_prefactor

class MPChain:
    def __init__(self,device,bias,temp,digits=50):
        self.ctx=mp.mp.clone();self.ctx.dps=digits;m=self.ctx
        self.t=m.mpf(float(device.hopping));self.t2=self.t**2
        self.diag=[m.mpf(float(x)) for x in (2*device.hopping+device.profile(bias)[1])]
        self.bias=m.mpf(float(bias));self.vt=m.mpf(float(KB*temp/Q));self.mu=m.mpf(float(device.mu_eV))
        self.pref=m.mpf(float(current_prefactor(device.mass)));self.calls=0
        self.numerator=self.t**(2*(len(self.diag)-1))
    def transmission(self,energy):
        self.calls+=1;m=self.ctx;e=m.mpf(energy);t=self.t
        def sigma(edge):
            x=e-edge-2*t
            return (x-m.j*m.sqrt(4*t*t-x*x))/2
        sl=sigma(m.mpf(0));sr=sigma(-self.bias)
        prev=m.mpc(1);cur=e-self.diag[0]-sl
        for j in range(1,len(self.diag)):
            a=e-self.diag[j]-(sr if j==len(self.diag)-1 else 0)
            prev,cur=cur,a*cur-self.t2*prev
        return (-2*sl.imag)*(-2*sr.imag)*self.numerator/abs(cur)**2
    def current_integrand(self,energy):
        m=self.ctx;e=m.mpf(energy);v=self.vt
        supply=v*(m.log1p(m.exp((self.mu-e)/v))-m.log1p(m.exp((self.mu-self.bias-e)/v)))
        return self.pref*self.transmission(e)*supply


def integrate_window(chain,left,right,center,width,orders=(10,20),rtol=1e-9):
    """All energy arithmetic uses MP; Gaussian nodes/weights have double accuracy.

    Geometric panel edges are formed in MP so a linewidth smaller than a
    float energy ulp would not collapse the local mapping (not studied here).
    """
    m=chain.ctx;c=m.mpf(float(center));w=m.mpf(float(width));l=m.mpf(float(left));r=m.mpf(float(right))
    edges=[l,r,c]
    for f in (.25,1,4,16,64):
        edges += [c-w*f,c+w*f]
    edges=sorted(set(x for x in edges if l<=x<=r))
    rules=[leggauss(n) for n in orders];calls0=chain.calls
    def panel(a,b):
        mid=(a+b)/2;half=(b-a)/2;ints=[]
        for nodes,weights in rules:
            ints.append(half*sum(m.mpf(float(weight))*chain.current_integrand(mid+half*m.mpf(float(node)))
                                 for node,weight in zip(nodes,weights)))
        return ints[1],abs(ints[1]-ints[0])
    stack=[(a,b,0) for a,b in zip(edges[:-1],edges[1:])];value=m.mpf(0);error=m.mpf(0)
    while stack:
        a,b,depth=stack.pop();v,e=panel(a,b)
        if e<=rtol*abs(v)+m.mpf('1e-15'):
            value+=v;error+=e
        elif depth<20:
            mid=(a+b)/2;stack.extend([(a,mid,depth+1),(mid,b,depth+1)])
        else:raise RuntimeError('MP local refinement limit')
    return float(value),float(error),chain.calls-calls0
