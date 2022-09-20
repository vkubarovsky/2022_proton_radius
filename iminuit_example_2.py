import numpy as np
from scipy.stats import norm
from iminuit import Minuit, cost
truth=100., 200., 0.3, 0.1, 0.7, 0.2
def scaled_cdf(xe,n1,n2,mu1,sigma1,mu2,sigma2):
       return n1*norm.cdf(xe,mu1,sigma1) + n2*norm.cdf(xe,mu2,sigma2)
xe=np.linspace(0,1)
m = np.diff(scaled_cdf(xe, *truth))
n=np.random.default_rng(1).poisson(m)  # generate random distribution

c=cost.ExtendedBinnedNLL(n, xe, scaled_cdf)
m=Minuit(c, *truth)
m.interactive()
