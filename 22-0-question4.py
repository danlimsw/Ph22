from math import *
from RootFinder import *
from numpy import *
from pylab import *
def radvelplot(phi,phase,N):
    # Plots the radial velocity as a function of fraction of orbit period.
    # phi = Angle of line-of-sight to the Earth
    # phase = Initial phase of the binary (radians)
    # N = Number of steps
    T = 27906.98161
    ecc = 0.617139
    a = 2.34186*299792.458 #km/s
    zetalist = []
    N = abs(int(N))
    deltat = 2.0*T/N #two cycles modelled
    for i in range(N):
        zeta = newtonsimple('(%f/(2*pi))*(x-%f*sin(x))-%f' %(T,ecc,deltat*i+phase*T/(2*pi)),'(%f/(2*pi))*(1-%f*cos(x))'%(T,ecc),0.2,0.00000001)
        zetalist = zetalist + [zeta]
    xlist = [a*(cos(zetalist[i])-ecc) for i in range(N)]
    ylist = [a*sqrt(1-ecc**2)*sin(zetalist[i]) for i in range(N)]
    vx = [(xlist[i+1]-xlist[i])/deltat for i in range(N-1)]
    vy = [(ylist[i+1]-ylist[i])/deltat for i in range(N-1)]
    tlist = [(2.0/N)*i for i in range(N-1)]
    radvel = [vx[i]*cos(phi)+vy[i]*sin(phi) for i in range(N-1)]
    plot(tlist,radvel)
    xlabel('Relative Phase')
    ylabel('Radilal Velocity, km/s')
    show()

