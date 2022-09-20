#!/usr/local/bin/python3
import numpy as np
from matplotlib import pyplot as plt
import math
import os
from scipy.optimize import curve_fit

#
#  D. Kharzeev "The mas radius of the proton", 
#  Hall-C
#  Gluex
#

hc=0.1973    #   GeV fm
pi=math.pi
alpha=1/137
e2=4*pi*alpha

Mp=0.938
Mjpsi=3.097


def G(ms,t):
    return Mp/(1-t/ms**2)**2

def m_s(rm):
    return math.sqrt(12.)*hc/rm

def dsdt(E,c2,rm,t):
    s=Mp**2+2*E*Mp
    Ecm=(s-Mp**2)/(2*math.sqrt(s))
    ms=m_s(rm)
    G2=G(ms,t)**2
    b=9
    return 1/(64*pi*s)/Ecm**2*e2*(2/3*c2)**2*(16*pi**2*Mp/b)**2*G2

def dsdt0(const0,ms,t):
    return const0/(1-t/ms**2)**4


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
 

    
#     HALL-C 

lab=np.array([ ' 9.10', ' 9.25', ' 9.40', ' 9.55', ' 9.70', ' 9.85','10.00','10.15','10.30','10.45']) #,10.72])
Emini=np.array([ 9.10, 9.25, 9.40, 9.55, 9.70, 9.85,10.00,10.15,10.30,10.45]) #,10.72])
Emaxi=np.array([ 9.25, 9.40, 9.55, 9.70, 9.85,10.00,10.15,10.30,10.45,10.60]) #,10.72])
msi  =np.array([ 1.69, 2.30, 1.59, 1.52, 1.18, 1.26, 1.43, 1.31, 1.33, 1.32]) #, 1.24])
smsi =np.array([ 0.45, 0.47, 0.22, 0.17, 0.16, 0.17, 0.16, 0.14, 0.19, 0.20])
c2i  =np.array([34.10,27.37,46.34,53.57,75.39,70.06,63.75,72.74,70.29,75.32]) #, 32.0])
sc2i =np.array([10.60, 4.30, 7.35, 7.63,16.88,14.65, 9.24,10.20,12.07,13.59])
rmi  =np.array([0.404,0.297,0.430,0.450,0.579,0.544,0.477,0.520,0.515,0.517]) #, 0.55])
srmi =np.array([0.107,0.061,0.059,0.052,0.079,0.075,0.053,0.057,0.074,0.079])
ds   =rmi*0
xstot=rmi*0
Ei   =rmi*0


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
#
#     PLOTS of my fit of the ds/dt CLAS12 data  
#
tdata=np.linspace(0,5.5,100)

plt.figure (num=10,dpi=120)
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
pdf_file='dsdt_CLAS1.pdf'
plt.savefig(pdf_file)
os.system('open -a Preview '+pdf_file)
###
plt.figure (num=20,dpi=120)
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
pdf_file='dsdt_CLAS2.pdf'
#
plt.savefig(pdf_file)
os.system('open -a Preview '+pdf_file)
         

#
#   Hall-C ds/dt
#
plt.figure (num=0,dpi=120)
plt.title  (r'Hall C   $\frac{d\sigma}{dt}$')
plt.xlabel (r'$-t, GeV^2$')
plt.ylabel (r'$\frac{d\sigma}{dt},nb/GeV^2$')
plt.grid()
#plt.yscale("log")

j = np.arange(10)
for i in j:
    E  = (Emini[i]+Emaxi[i])/2
    Ei[i] = E
    c2 = c2i[i]
    rm = rmi[i]
    x0=np.linspace(-tmin(E),5.5,1000)
    y0=dsdt(E,c2,rm,-x0)
    xstot[i]=int(E,c2,rm)
    plt.plot(x0,y0,label=lab[i])
#plt.plot(x0,y0,label="10.72,Kharzeev")
#plt.ylim(bottom=0.01,top=20)
#plt.legend(["blue", "green","red"], loc ="upper right")
#plt.legend(prop={'size': 10})
plt.yscale("log")
plt.legend()
pdf_file='dsdt_hallc.pdf'
plt.savefig(pdf_file)
#os.system('open -a Preview '+pdf_file)

#
# Gluex ds/dt
#
plt.figure (num=1,dpi=120)
plt.title  (r'Gluex $\frac{d\sigma}{dt}$')
plt.xlabel (r'$-t, GeV^2$')
plt.ylabel (r'$\frac{d\sigma}{dt},nb/GeV^2$')
plt.yscale("log")
plt.grid()
#plt.xlim(0.0,10)
#plt.ylim(0.001,3)
EiG   = np.array([ 8.930, 9.860,10.820])
CiG   = np.array([ 3.121, 2.303, 4.184])
sCiG  = np.array([ 2.230, 0.400, 0.541])
msiG  = np.array([ 1.089, 1.453, 1.314])
smsiG = np.array([ 0.172, 0.074, 0.049])
rmiG  = math.sqrt(12.)/msiG*hc
srmiG = rmiG*smsiG/msiG
lab   = np.array(['8.93 GeV2', '9.86 GeV2', '10.82 GeV2'])
xstotG=EiG*0.
j = np.arange(3)
for i in j:
    E    = EiG[i]
    const= CiG[i]
    msG  = msiG[i]
    x0=np.linspace(-tmin(E),-tmax(E),1000)
    y0=dsdt0(const,msG,-x0)
    xstotG[i]=int0(E,const,msG)
    plt.plot(x0,y0, label=lab[i])
plt.legend()
pdf_file='dsdt_Gluex.pdf'
plt.savefig(pdf_file)
#os.system('open -a Preview '+pdf_file)

#
#  ms Fit parameter
#
plt.figure (num=3,dpi=120)
plt.title  ('ms Fit parameter')
plt.xlabel (r'$E_\gamma, GeV$')
plt.ylabel ('ms Fit Parameter')
plt.grid()
#plt.ylim(bottom=0,top=10)
plt.errorbar(EiG,msiG,yerr=smsiG, marker='o', color='Red', label='Gluex')
plt.errorbar(Ei, msi, yerr=smsi,  marker='o', color='Blue',label='Hall-C')
#plt.errorbar(ED_CLAS, ms_d, yerr=sms_d,  marker='o', color='Black',label='Me')
plt.errorbar(ED_CLAS, ms_J, yerr=sms_J,  marker='o', color='Black',label='CLAS12')
plt.legend()
pdf_file = 'ms.pdf'
plt.savefig(pdf_file)
#os.system('open -a Preview '+pdf_file)

#
#  Proton Radius r=sqrt(12)/ms*hc or sqrt(6)/<t>*hc
#
plt.figure (num=4,dpi=120)
plt.title  ('Proton Radius, fm')
plt.xlabel (r'$E_\gamma, GeV$')
plt.ylabel ('Proton Radius, fm')
plt.grid()
#plt.ylim(bottom=0,top=10)
plt.errorbar(EiG,rmiG,yerr=srmiG,marker='o',           color='Red' , label='Gluex')
plt.errorbar(Ei, rmi, yerr=srmi, marker='o',           color='Blue',   label='Hall-C')
plt.errorbar(ED_CLAS, rm_CLAS, yerr=srm_CLAS, fmt='o', color='fuchsia',label='CLAS12 Exponential')
#plt.errorbar(ED_CLAS+0.05, rm_d, yerr=srm_d, fmt='o',  color='Black',    label='CLAS12 Dipole')
plt.errorbar(ED_CLAS+0.05, rm_J, yerr=srm_J, fmt='o',  color='aqua',    label='CLAS12 Dipole')
plt.legend(loc='upper left')
pdf_file = 'rm.pdf'
os.system('rm '+pdf_file)
plt.savefig(pdf_file)
#os.system('open -a Preview '+pdf_file)

#
#####     Proton radius Weighted average
#

j=np.arange(10)
sum1=0; sum2=0
for i in j:
    sum1=rmi[i]/srmi[i]**2+sum1
    sum2=1/srmi[i]**2+sum2
RM_C=sum1/sum2
SRM_C=1/sum2**0.5

j=np.arange(3)
sum3=0; sum4=0
for i in j:
    sum3=rmiG[i]/srmiG[i]**2+sum3
    sum4=1/srmiG[i]**2+sum4
RM_G=sum3/sum4
SRM_G=1/sum4**0.5

sum5=sum1+sum3
sum6=sum2+sum4
RM=sum5/sum6
SRM=1/sum6**0.5

rme_CLAS=( rm_CLAS[0]/srm_CLAS[0]**2 + rm_CLAS[1]/srm_CLAS[1]**2 )/(1/srm_CLAS[0]**2+1/srm_CLAS[1]**2)
srme_CLAS=(1/(1/srm_CLAS[0]**2 + 1/srm_CLAS[1]**2))**0.5

rmd_CLAS=( rm_d[0]/srm_d[0]**2 + rm_d[1]/srm_d[1]**2 )/(1/srm_d[0]**2+1/srm_d[1]**2)
srmd_CLAS=(1/(1/srm_d[0]**2 + 1/srm_d[1]**2))**0.5

rmd_J=( rm_J[0]/srm_J[0]**2 + rm_J[1]/srm_J[1]**2 )/(1/srm_J[0]**2+1/srm_J[1]**2)
srmd_J=(1/(1/srm_J[0]**2 + 1/srm_J[1]**2))**0.5
                                                    

###################################################################


#  Total xs

#plt.ylim(bottom=0,top=10)

#    Gluex
sxstotG=xstotG*((sCiG/CiG)**2+0*(smsiG/msiG))**0.5

#    Gluex 18 points total xs
#   E       xs      ds       Emin Emax
#  0  8.29  0.04411 0.02869 8.2 8.38
#  1  8.47  0.13752 0.02549 8.38 8.56
#  2  8.65  0.25215 0.02842 8.56 8.74
#  3  8.83  0.32990 0.01463 8.74 8.92
#  4  9.01  0.22654 0.05680 8.92 9.1
#  5  9.19  0.18879 0.02169 9.1 9.28
#  6  9.37  0.51052 0.02281 9.28 9.46
#  7  9.55  0.78648 0.06868 9.46 9.64
#  8  9.73  0.51886 0.02176 9.64 9.82
#  9  9.91  0.70360 0.11635 9.82 10
# 10 10.09  0.79588 0.07384 10 10.18
# 11 10.27  0.81328 0.04666 10.18 10.36
# 12 10.45  1.37679 0.06611 10.36 10.54
# 13 10.63  0.96900 0.11616 10.54 10.72
# 14 10.81  1.14295 0.05098 10.72 10.9
# 15 10.99  1.09218 0.03088 10.9 11.08
# 16 11.17  1.56346 0.12898 11.08 11.26
# 17 11.35  1.84710 0.02760 11.26 11.44

EgGluex=np.array([ 8.29,   8.47,   8.65,   8.83,   9.01,   9.19,   9.37,   9.55,   9.73,   9.91,  10.09,  10.27,  10.45,  10.63,  10.81,  10.99,  11.17,  11.35])
xsGluex=np.array([ 0.04411,0.13752,0.25215,0.32990,0.22654,0.18879,0.51052,0.78648,0.51886,0.70360,0.79588,0.81328,1.37679,0.96900,1.14295,1.09218,1.56346,1.84710])
sxsGluex=np.array([0.02869,0.02549,0.02842,0.01463,0.05680,0.02169,0.02281,0.06868,0.02176,0.11635,0.07384,0.04666,0.06611,0.11616,0.05098,0.03088,0.12898,0.02760])

#    Hall-C
sxstot=xstot*((2*sc2i/c2i)**2+0*(smsi/msi))**0.5


#    CLAS12 total cross section

EgCLAS = np.array([8.714,   9.128,    9.526,  9.929, 10.359])
xsCLAS = np.array([0.0581, 0.136346, 0.1984, 0.2604, 0.4384])
sxsCLAS= np.array([0.0159,  0.0548,   0.0411, 0.0334, 0.0748])

plt.figure (num=2,dpi=120)
plt.title  (r'$\sigma_{total}$')
plt.xlabel (r'$E_\gamma, GeV$')
plt.ylabel (r'$\sigma_{total},nb$')
plt.grid()


plt.errorbar(EgGluex,xsGluex, yerr=sxsGluex, marker='o', color='lime',  label='Gluex preliminary')
plt.errorbar(EiG,xstotG, yerr=sxstotG,       marker='o', color='Red',   label=r'Gluex $\frac{d\sigma}{dt}$ integrated')
plt.errorbar(Ei,  xstot, yerr=sxstot,        marker='o', color='Blue',  label=r'Hall-C $\frac{d\sigma}{dt}$ integrated')
plt.errorbar(EgCLAS, xsCLAS, yerr=sxsCLAS,   marker='o', color='Black', label='CLAS12')
#plt.errorbar(ED_CLAS, xsD_CLAS,sxsD_CLAS,marker='o',     color='Grey',  label='CLAS12 integrated')
plt.legend()
pdf_file = 'xs_total.pdf'
plt.savefig(pdf_file)
#os.system('open -a Preview '+pdf_file)

####   xs Total Log scale

plt.figure (num=5,dpi=120)
plt.title  (r'$\sigma_{total}$')
plt.xlabel (r'$E_\gamma, GeV$')
plt.ylabel (r'$\sigma_{total},nb$')
plt.grid()

plt.yscale("log")
plt.ylim(0.03,2.5)

plt.errorbar(EgGluex,xsGluex, yerr=sxsGluex, marker='o',color='lime',label='Gluex preliminary')
plt.errorbar(EiG,xstotG, yerr=sxstotG, marker='o', color='Red',  label=r'Gluex $\frac{d\sigma}{dt}$ integrated')
plt.errorbar(Ei,  xstot, yerr=sxstot, marker='o', color='Blue', label=r'Hall-C $\frac{d\sigma}{dt}$ integrated')
plt.errorbar(EgCLAS, xsCLAS, yerr=sxsCLAS, marker='o', color='Black', label='CLAS12')
#plt.errorbar(ED_CLAS, xsD_CLAS,sxsD_CLAS,marker='o',color='Grey',label='CLAS12 integrated')

plt.legend()

pdf_file = 'xs_total_log_1.pdf'
plt.savefig(pdf_file)
#os.system('open -a Preview '+pdf_file)


# Compare fits



    

    

