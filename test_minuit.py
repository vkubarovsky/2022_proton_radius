#!/usr/local/bin/python3
import numpy as np
from matplotlib import pyplot as plt
import math
import os
from scipy.optimize import curve_fit

hc=0.1973    #   GeV fm
pi=math.pi
alpha=1/137
e2=4*pi*alpha

Mp=0.938
Mjpsi=3.097

def tmin(E):
    m1=0
    m2=Mp
    m3=Mjpsi 
    m4=Mp
    s=Mp**2+2*Mp*E
    ssqrt=math.sqrt(s)
    E1cm=(s+m1**2-m2**2)/2/ssqrt
    E3cm=(s+m3**2-m4**2)/2/ssqrt
    P1cm=math.sqrt(E1cm**2-m1**2)
    P3cm=math.sqrt(E3cm**2-m3**2)
    return ((m1**2-m3**2-m2**2+m4**2)/2/ssqrt)**2-(P1cm-P3cm)**2
def tmax(E):
    m1=0
    m2=Mp
    m3=Mjpsi 
    m4=Mp
    s=Mp**2+2*Mp*E
    ssqrt=math.sqrt(s)
    E1cm=(s+m1**2-m2**2)/2/ssqrt
    E3cm=(s+m3**2-m4**2)/2/ssqrt
    P1cm=math.sqrt(E1cm**2-m1**2)
    P3cm=math.sqrt(E3cm**2-m3**2)
    return ((m1**2-m3**2-m2**2+m4**2)/2/ssqrt)**2-(P1cm+P3cm)**2

def int(E,c2,rm):
    pi=3.1415
    si=Mp**2+2*E*Mp
    Ecmi=(si-Mp**2)/(2*si**0.5)
    ms=math.sqrt(12.)*0.1973/rm
    G2i=Mp**2*ms**8/3.*(1/(-tmin(E)+ms**2)**3-1/(-tmax(E)+ms**2)**3)
    b=9
    return 1/(64*pi*si)/Ecmi**2*e2*(2/3*c2)**2*(16*pi**2*Mp/b)**2*G2i 

def int0(E,const0,ms):
    return const0*ms**8/3.*(1/(-tmin(E)+ms**2)**3-1/(-tmax(E)+ms**2)**3)

def exp_i (x,p0,p1):
    return p0*p1*np.exp(-p1*x)

def exp_t (x,p0,p1):
    return p0*p1*np.exp(-p1*x)/(math.exp(p1*tmin(Ebeam))-math.exp(p1*tmax(Ebeam)))

def dipole (x,p0,p1):
     Norm=p1**8/3*(1/(-tmin(Ebeam)+p1**2)**3-1/(-tmax(Ebeam)+p1**2)**3)
     return p0/(1+x/p1**2)**4/Norm


#    CLAS12 ds/dt

t1_CLAS   = np.array([0.301879, 0.703882,  1.19263,   1.90433])
ds1_CLAS  = np.array([0.237853, 0.125103,  0.101574,  0.025458])
sds1_CLAS = np.array([0.0428135,0.0262715, 0.0241746, 0.00641551])

t2_CLAS   = np.array([0.259479,  0.671546,  1.346,      1.78621])
ds2_CLAS  = np.array([0.397352,  0.210328,  0.0532372, 0.0903125])
sds2_CLAS = np.array([0.0953645, 0.0504788, 0.0175683, 0.0234813])

EminD_CLAS=np.array([ 9.7, 10.2])
EmaxD_CLAS=np.array([10.2, 10.6])
ED_CLAS=(EminD_CLAS+EmaxD_CLAS)/2

xsD_CLAS =np.array([0.2744, 0.3912])
sxsD_CLAS=np.array([0.04346,0.07381])

pexp_CLAS= np.array([1.283, 1.405])
spexp_CLAS=np.array([0.2953,0.3118])
#                xs     slope

rm_CLAS  = (6*pexp_CLAS)**0.5 * hc
srm_CLAS = rm_CLAS/2 * (spexp_CLAS/pexp_CLAS)

tdip1=t1_CLAS-tmin(ED_CLAS[0])
tdip2=t2_CLAS-tmin(ED_CLAS[1])

tdata=np.linspace(0,5.5,100)

#  My own fit of the ds/d(t-tmin) CLAS12 data by exp functin
popt3, pcov3 = curve_fit(exp_i, t1_CLAS, ds1_CLAS,sigma=sds1_CLAS,absolute_sigma=True)
popt4, pcov4 = curve_fit(exp_i, t2_CLAS, ds2_CLAS,sigma=sds2_CLAS,absolute_sigma=True)

err3=np.sqrt(np.diag(pcov3))
err4=np.sqrt(np.diag(pcov4))
#   Use t not t-tmin
Ebeam=ED_CLAS[0]
popt1, pcov1 = curve_fit(exp_t, tdip1, ds1_CLAS,sigma=sds1_CLAS,absolute_sigma=True)
Ebeam=ED_CLAS[1]
popt2, pcov2 = curve_fit(exp_t, tdip2, ds2_CLAS,sigma=sds2_CLAS,absolute_sigma=True)
err1=np.sqrt(np.diag(pcov1))
err2=np.sqrt(np.diag(pcov2))

#   My own fit of the ds/dt CLAS12 data by dipole formula

Ebeam=ED_CLAS[0]
dip1, sdip1 = curve_fit(dipole, tdip1, ds1_CLAS,sigma=sds1_CLAS,absolute_sigma=True)
Ebeam=ED_CLAS[1]
dip2, sdip2 = curve_fit(dipole, tdip2, ds2_CLAS,sigma=sds2_CLAS,absolute_sigma=True)
errd1=np.sqrt(np.diag(sdip1))
errd2=np.sqrt(np.diag(sdip2))

#   ms fit parameter to calculate proton radius later
ms_d=ED_CLAS*0; sms_d=ED_CLAS*0
ms_d[0]=dip1[1]
ms_d[1]=dip2[1]
sms_d[0] = errd1[1]
sms_d[1] = errd2[1]

#   Proton Radius
rm_d=(12)**0.5/ms_d*hc
srm_d=rm_d*sms_d/ms_d

#   Joseph fit by dipole formula
ms_J=np.array([1.195,1.271])
sms_J=np.array([0.2986,0.2811])
rm_J=(12)**0.5/ms_J*hc
srm_J=rm_J*sms_J/ms_J


plt.figure (num=10)
plt.title  (r'CLAS12   $\frac{d\sigma}{dt}$  9.95 $GeV^2$')
plt.xlabel (r'$-t, GeV^2$')
plt.ylabel (r'$\frac{d\sigma}{dt},nb/GeV^2$')

plt.errorbar(tdip1,ds1_CLAS,yerr=sds1_CLAS,fmt='o', color='Black',label='Exp')
#plt.errorbar(tdip1,ds1_CLAS,yerr=sds1_CLAS,fmt='o',  color='Blue',label='Dip[ole')
plt.plot(tdata, exp_t(tdata, *popt1), 'r-', label=\
         'Exp    : p0='+f'{popt1[0]:.4f}' + '+/-' + f'{err1[0]:.4f}'+'  '+\
         'p1='+f'{popt1[1]:.4f}' + '+/-' + f'{err1[1]:.4f}')

plt.plot(tdata, dipole(tdata, *dip1), 'b-', label=\
         'Dipole: p0='+f'{dip1[0]:.4f}' + '+/-' + f'{errd1[0]:.4f}'+'  '+\
         'p1='+f'{dip1[1]:.4f}' + '+/-' + f'{errd1[1]:.4f}')

plt.xlabel('t')
plt.ylabel(r'$d\sigma/dt$')
plt.ylim(0,0.3)
plt.grid()
plt.legend()
pdf_file='dsdt_CLAS1_minuit.pdf'
plt.savefig(pdf_file)
#os.system('open -a Preview '+pdf_file)

###

plt.figure (num=20)
plt.title  (r'CLAS12   $\frac{d\sigma}{dt}$  $10.4 GeV^2$')
plt.xlabel (r'$-t, GeV^2$')
plt.ylabel (r'$\frac{d\sigma}{dt},nb/GeV^2$')

plt.errorbar(tdip2,ds2_CLAS,yerr=sds2_CLAS,fmt='o', color='Black',label='Exp')
#plt.errorbar(tdip2,ds2_CLAS,yerr=sds2_CLAS,fmt='o', color='Blue',label=r'Dipole')

plt.plot(tdata, exp_t(tdata, *popt4), 'r-', label=\
         'Exp    : p0='+f'{popt2[0]:.4f}' + '+/-' + f'{err2[0]:.4f}'+'  '+\
         'p1='+f'{popt2[1]:.4f}' + '+/-' + f'{err2[1]:.4f}')

plt.plot(tdata, dipole(tdata, *dip2), 'b-', label=\
         'Dipole: p0='+f'{dip2[0]:.4f}' + '+/-' + f'{errd2[0]:.4f}'+'  '+\
         'p1='+f'{dip2[1]:.4f}' + '+/-' + f'{errd2[1]:.4f}')

plt.xlabel('t-tmin')
plt.ylabel(r'$d\sigma/dt$')
plt.ylim(0,0.5)
plt.grid()
plt.legend()
pdf_file='dsdt_CLAS2_minuita.pdf'
#
plt.savefig(pdf_file)
#os.system('open -a Preview '+pdf_file)


#
######  MINUIT
#
# everything in iminuit is done through the Minuit object, so we import it
from iminuit import Minuit
def dipmin (x,p0,p1):
     Norm=p1**8/3*(1/(-tmin(Ebeam)+p1**2)**3-1/(-tmax(Ebeam)+p1**2)**3)
     return p0/(1+x/p1**2)**4/Norm

def cost(p0,p1):
    ym = dipmin(tdip1, p0, p1)
    res = (ds1_CLAS - ym)**2 / sds1_CLAS*2 # studentized residuals
    return np.sum(res ** 2)
#
# least-squares score function = sum of data residuals squared
Ebeam=9.95
p0=1. ; p1=1.
def LSQ(x,p0, p1):
    return np.sum((ds1_CLAS - dipole(x, p0, p1)) ** 2 / sds1_CLAS ** 2)

# set start values via keywords for p0 and p1
#m = Minuit(dipole, p0=1., p1=1., error_p0=0.1, error_p1=0.1, errordef=1)
m = Minuit(cost, p0=1., p1=1.)
m.params


