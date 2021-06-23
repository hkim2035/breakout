import sys
import os.path
import math
from zipfile import ZIP_DEFLATED
import numpy as np
from numpy.lib.npyio import zipfile_factory

def rock_strength_parameter(fric, coh):
    Co = 2*coh*(math.sqrt(fric**2+1)+fric)  # uniaxial compressive strength
    q = (math.sqrt(fric**2+1)+fric)**2
    C1 = (1.+0.6*fric)*Co                   # biaixal plane strength (Wiebols and Cook, 1968)

    Xa = (q-1.)/(math.sqrt(3.)*(2.+q))
    Xk = Co*math.sqrt(3.)/(2.+q)
    T = Co/12.                              # tensile strength

    return Co, T

def pstress(r, zeta, S1, S2, S3, So):
    
    global SHmax, SHmin, Sv, nu, Pm, Po

    theta = 2.*zeta*3.141593/180.
    cos1 = math.cos(theta)
    sin1 = math.sin(theta)

    So = 0.5*(SHmax + SHmin)*(1+r**2)
    So = So - 0.5*(SHmax-SHmin)*(1+3*r**4)*cos1
    So = So - r**2*(Pm-Po)

sys.argv

if os.path.isfile(sys.argv[1]) == False :
    print("Input file error.")
    exit

f = open(sys.argv[1],"r")
input = f.readlines()

if len(input) != 16:
    print("Input file format error.")
    exit

fric = float(input[8])
coh = float(input[9])
nu = float(input[10])
Sv = float(input[11])
SHmax = float(input[12])
SHmin = float(input[13])
Po = float(input[14])
Pm = float(input[15])

drad = 0.01
const = 3.141593/180.
radius = 5. #cm
xl = 7.     #cm
yl = 7.     #cm

zeta = 90 - np.linspace(0,90,num=91)
angl = zeta 

rad22 = []
for x in zeta:
    if x >= 45:
        tmp = yl/math.sin(x*3.141593/180.)
    else:
        tmp = xl/math.cos(x*3.141593/180.)
    rad22.append(tmp)

Nrad2 = [(yy/radius - 1.)/drad for yy in rad22]
Nrad3 = 6./drad
r = [radius/yy for yy in rad22]

# calculate stress field
def pstress(r, zeta, Sl, S2, S3, So):

    global SHmax, SHmin, Sv, nu, Pm, Po

    theta = 2.*zeta*3.141593/180.
    cos1 = math.cos(theta)
    sin1 = math.sin(theta)

    So = 0.5*(SHmax+SHmin) * (1.+r**2)
    So = So - 0.5*(SHmax-SHmin)*(1.+3*r**4)*cos1
    So = So - r**2*(Pm-Po)

    Sr = 0.5*(SHmax+SHmin)*(1.-r**2)
    Sr = Sr + 0.5*(SHmax-SHmin)*cos1*(1.-4*r**2+3*r**4)
    Sr = Sr + r**2*(Pm-Po)
    Sro = -0.5*(SHmax-SHmin)*sin1*(1.+2*r**2-3*r**4)
    Sz = Sv - 2*nu*(SHmax-SHmin)*cos1*r**2
    Sp1 = 0.5*(So+Sr)
    Sp1 = Sp1 + 0.5*math.sqrt((So-Sr)**2+4*Sro**2)
    Sp2 = 0.5*(So+Sr)
    Sp2 = Sp2 - 0.5*math.sqrt((So-Sr)**2+4*Sro**2)

    if Sp2 >= Sz: 
        S1 = Sp1
        S2 = Sp2
        S3 = Sz

    if (Sp1 >= Sz) & (Sz > Sp2):
        S1= Sp1
        S2 = Sz
        S3 = Sp2

    if Sz >= Sp1:
        S1 = Sz
        S2 = Sp1
        S3 = Sp2

# calculate rock stress


