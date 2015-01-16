"""
This function numerically solves the
simple harmonic oscilator equation
using k/m=1. The input is the start step
size h, and the methods used are the explicit
Euler method and midpoint method. 
h = Starting Step size (time)
The number of subintervals is chosen
so that three cycles are plotted.
The evaluation is performed for step-sizes
h, h/2, h/4, h/8 and h/16. The error
(given by x_analytic - x for each step)
is integrated over the three cycles to
give the global error. The error at each
step-size is then plotted for each method
on the same axis. 
"""
from math import *
from numpy import *
from pylab import *
import sys
def eulerrun(h):
    #main loop for Explicit Euler Method
    #initial conditions x=1(max), v=0
    x = 1
    v = 0
    xlist = [1]
    vlist = [0]
    h = min(h,0.1)
    N = int((3*2*pi)/h)
    det = 1 + h**2
    for n in range(N):
        x1 = x + h*v
        v1 = v - h*x
        xlist = xlist + [x1]
        vlist = vlist + [v1]
        x = x1
        v = v1
    t = [h*n for n in range(N+1)]
    # analytic solution is just cos(t)
    anx = [cos(t[i]) for i in range(N+1)]
    globalerror = abs(sum([anx[i]-xlist[i] for i in range(N+1)]))
    return globalerror
def midpointrun(h):
    #main loop for Midpoint Method
    #initial conditions x=1(max), v=0
    x = 1
    v = 0
    xlist = [1]
    vlist = [0]
    h = min(h,0.1)
    N = int((3*2*pi)/h)
    for n in range(N):
        xmid = x + (h/2)*v
        vmid = v - (h/2)*x
        x1 = x + h*vmid
        v1 = v - h*xmid
        xlist = xlist + [x1]
        vlist = vlist + [v1]
        x = x1
        v = v1
    t = [h*n for n in range(N+1)]
    # analytic solution is just cos(t)
    anx = [cos(t[i]) for i in range(N+1)]
    globalerror = abs(sum([anx[i]-xlist[i] for i in range(N+1)]))
    return globalerror
if len(sys.argv)==2:
    [name,h] = sys.argv
    h=abs(float(h))
    listeulerh = [eulerrun(h)]
    listmidpointh = [midpointrun(h)]
    for i in range(4):
        h = h/2
        listeulerh = listeulerh + [eulerrun(h)]
        listmidpointh = listmidpointh + [midpointrun(h)]
    eulerplot = plot([h,h/2,h/4,h/8,h/16],listeulerh,'bo',label='Euler')
    plot([h,h/2,h/4,h/8,h/16],listeulerh,'b-')
    midpointplot = plot([h,h/2,h/4,h/8,h/16],listmidpointh,'ro',label='Midpoint')
    plot([h,h/2,h/4,h/8,h/16],listmidpointh,'r-')
    xticks([h,h/2,h/4,h/8,h/16], ['h','h/2','h/4','h/8','h/16'])
    xlabel('Step Size h')
    ylabel('Global error')
    legend()
    show()
else:
    print 'Incorrect number of arguments'
