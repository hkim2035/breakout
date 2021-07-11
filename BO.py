import sys
import os.path
import math
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def rock_strength_parameter(fric, coh):
    Co = 2*coh*(math.sqrt(fric**2+1)+fric)  # uniaxial compressive strength
    q = (math.sqrt(fric**2+1)+fric)**2
    C1 = (1.+0.6*fric)*Co                   # biaixal plane strength (Wiebols and Cook, 1968)

    Xa = (q-1.)/(math.sqrt(3.)*(2.+q))
    Xk = Co*math.sqrt(3.)/(2.+q)
    T = Co/12.                              # tensile strength

    return Co, T

# calculate principal stress field around wellbore
def pstress(radius, r, zeta):

    # r = a/r. r is the radial coordinate system and a is the wellbore radius.

    global SHmax, SHmin, Sv, nu, Pm, Po

    eSHmax = SHmax - Po
    eSHmin = SHmin - Po
    eSv = Sv - Po

    theta = 2.*math.radians(zeta)
    cos1 = math.cos(theta)
    sin1 = math.sin(theta)

    So = 0.5*(eSHmax+eSHmin) * (1.+radius**2/r**2)
    So = So - 0.5*(eSHmax-eSHmin)*(1.+3*radius**4/r**4)*cos1
    So = So - radius**2/r**2*(Pm-Po)

    Sr = 0.5*(eSHmax+eSHmin)*(1.-radius**2/r**2)
    Sr = Sr + 0.5*(eSHmax-eSHmin)*cos1*(1.-4*radius**2/r**2+3*radius**4/r**4)
    Sr = Sr + radius**2/r**2*(Pm-Po)

    Sz = eSv - 2*nu*(eSHmax-eSHmin)*cos1*radius**2/r**2

    Tau = -0.5*(eSHmax-eSHmin)*sin1*(1.+2*radius**2/r**2-3*radius**4/r**4)

    Sp1 = 0.5*(So+Sr)
    Sp1 = Sp1 + 0.5*math.sqrt((So-Sr)**2+4*Tau**2)
    Sp2 = 0.5*(So+Sr)
    Sp2 = Sp2 - 0.5*math.sqrt((So-Sr)**2+4*Tau**2)

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

    return [zeta, r, S1, S2, S3]


if __name__ == '__main__': 

    if os.path.isfile(sys.argv[1]) == False :
        print("Input file error.")
        exit

    f = open(sys.argv[1],"r")
    input = f.readlines()

    if len(input) != 18:
        print("Input file format error.")
        exit

    fric = float(input[9])
    coh = float(input[10])
    nu = float(input[11])
    radius = float(input[12])
    Sv = float(input[13])
    SHmax = float(input[14])
    SHmin = float(input[15])
    Po = float(input[16])
    Pm = float(input[17])

    plotr = 1.5 # wellbore radius a *

    # calculate stress field
    azimuths = np.radians(np.linspace(0, 91, 450))
    zeniths = np.arange(radius*5000, radius*plotr*5000, 1)/500
    
    stress = []
    for itheta in azimuths:
        for ir in zeniths:
            result = pstress(radius, ir, itheta)
            stress.append(result)
    
    ptheta = np.array(stress).T[0]
    pr = np.array(stress).T[1]
    pS1 = np.array(stress).T[2]
    pS2 = np.array(stress).T[3]
    pS3 = np.array(stress).T[4]
    
    bo_low = int(np.min(pS3))
    bo_upper = int(np.max(pS1))
    bo_mid = (bo_low + bo_upper)/2
    x1 = [pr[i]*math.cos(ptheta[i]) for i in range(0,len(pr))]
    y1 = [pr[i]*math.sin(ptheta[i]) for i in range(0,len(pr))]

    list = [x for x in zip(x1, y1, pS1, pS2, pS3)]

    fig = make_subplots(rows=1, cols=3)

    fig.append_trace(go.Scatter(x=x1, y=y1, mode='markers',
        marker=dict(size=1, color=pS1, showscale=True)),
        row=1, col=1)

    fig.append_trace(go.Scatter(x=x1, y=y1, mode='markers',
        marker=dict(size=1, color=pS2, showscale=False)),
        row=1, col=2)

    fig.append_trace(go.Scatter(x=x1, y=y1, mode='markers',
        marker=dict(size=1, color=pS3,  showscale=False)),
        row=1, col=3)

    fig.update_layout(
    autosize=False,
    width=1600,
    height=500,
    margin=dict(
        l=50,
        r=50,
        b=50,
        t=50,
        pad=4
    ),
    paper_bgcolor="White",
)

    fig.show()

    #cntr = plotly.go(plt.pcolormesh([x1, y1], pS1)
    #cntr = plt.tricontourf(x1, y1, pS1, levels=15, cmap="RdBu_r")
    #plt.colorbar(cntr)
    
    #fig, ax = plt.subplots()

    #plt.gca().set_aspect('equal', adjustable='box')
    #plt.draw()
    #plt.show()

