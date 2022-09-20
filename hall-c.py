#!/usr/local/bin/python3
import numpy as np
from matplotlib import pyplot as plt
import math
import os

#
#  H. W. Koch and J. W. Motz, Rev. Mod. Phys.31, 920 (1959).
#


def G(ms,t):
    return Mp/(1-t/ms**2)**2
def dsdt(E,c2,rm,t):
    pi=3.1415
    s=Mp**2+2*E*Mp
    Ecm=(s-Mp**2)/(2*math.sqrt(s))
    rmgev2=rm/0.1973
    ms=math.sqrt(12.)/rmgev2
    G2=G(ms,t)**2
    b=9
    return 1/(64*pi*s)/Ecm**2*(2/3*c2)**2*(16*pi**2*Mp/b)**2*G2 

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
    rmgev2=rm/0.1973
    ms=math.sqrt(12.)/rmgev2
    G2i=Mp**2*ms**8/3.*(1/(-tmin(E)+ms**2)**3-1/(-tmax(E)+ms**2)**3)
    b=9
    return 1/(64*pi*si)/Ecmi**2*(2/3*c2)**2*(16*pi**2*Mp/b)**2*G2i 

def int1(E,c2,rm):
    ms=math.sqrt(12.)/rm*0.1973
    return dsdt(E,c2,rm,tmin(E))/(G(ms,tmin(E))**2)*Mp**2*ms**8/3.*(1/(-tmin(E)+ms**2)**3-1/(-tmax(E)+ms**2)**3)

def fill_y(E,c2,rm,xx,yy):
    xx = np.linspace(-tmin(E),-tmax(E),1000)
    yy=xx*0
    for i in range(np.size(xx)):
        s=Mp**2+2*E*Mp
        Ecm=(s-Mp**2)/(2*math.sqrt(s))
        rmgev2=rm/0.1973
        ms=math.sqrt(12.)/rmgev2
        yy[i] = dsdt(E,c2,rm,-x[i])
        print(i,xx[i],yy[i])

Mp=0.938
Mjpsi=3.097

E=(10.72+10.72)/2; c2=32.00; rm=0.55
#x0=np.linspace(-tmin(E),-tmax(E),1000); y0=dsdt(E,c2,rm,-x0)
x0=np.linspace(-tmin(E),3.5,1000); y0=dsdt(E,c2,rm,-x0)
pi=3.1415
s=Mp**2+2*E*Mp
Ecm=(s-Mp**2)/(2*math.sqrt(s))
rmgev2=rm/0.1973
ms=math.sqrt(12.)/rmgev2
G2=G(ms,tmin(E))**2
b=9
xs=1/(64*pi*s)/Ecm**2*(2/3*c2)**2*(16*pi**2*Mp/b)**2*G2 


E=(9.10+9.25)/2; c2=34.00; rm=0.404
x1=np.linspace(-tmin(E),-tmax(E),1000); y1=dsdt(E,c2,rm,-x1)
int(E,c2,rm)

E=(9.70+9.85)/2; c2=75.39; rm=0.579
x2=np.linspace(-tmin(E),-tmax(E),1000); y2=dsdt(E,c2,rm,-x2)
int(E,c2,rm)

E=(10.45+10.60)/2; c2=75.32; rm=0.74
x3=np.linspace(-tmin(E),-tmax(E),1000); y3=dsdt(E,c2,rm,-x3)
int(E,c2,rm)

#E=10.72; c2=33.23; rm=0.55 # fm

plt.figure (num=0,dpi=120)
plt.title  (r'$\frac{d\sigma}{dt}$')
plt.xlabel (r'$-t, GeV^2$')
plt.ylabel (r'$\frac{d\sigma}{dt},nb/GeV^2$')
plt.grid()
#plt.yscale("log")

plt.plot(x0,y0,label="10.72,Kharzeev")
plt.plot(x1,y1,label="9.1-9.25")
plt.plot(x2,y2,label="9.7-9.85")
plt.plot(x3,y3,label="10.45-10.60")
plt.ylim(bottom=0.001)
plt.legend(["blue", "green","red"], loc ="upper right")
plt.legend(prop={'size': 10})

#plt.show()
#os.system("rm -f FF.pdf")
plt.savefig('FF.pdf')
os.system("open -a preview FF.pdf")



    

