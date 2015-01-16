from math import *
from RootFinder import *
from numpy import *
from pylab import *
T = 27906.98161 #seconds
ecc = 0.617139
a = 2.34186 #lightseconds
zetalist = []
N = 1000 # number of steps
for i in range(N):
    zeta = newtonsimple('(%f/(2*pi))*(x-%f*sin(x))-%f' %(T,ecc,100.0*i),'(%f/(2*pi))*(1-%f*cos(x))'%(T,ecc),0.2,0.00000001)
    zetalist = zetalist + [zeta]
    #print i
#tlist = [100.0*i for i in range(N)]
xlist = [a*(cos(zetalist[i])-ecc) for i in range(N)]
ylist = [a*sqrt(1-ecc**2)*sin(zetalist[i]) for i in range(N)]
plot(xlist,ylist)
xlabel('x-position, lightseconds')
ylabel('y-position,lightseconds')
show()

