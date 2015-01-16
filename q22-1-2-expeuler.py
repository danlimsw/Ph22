"""
This function numerically solves the
simple harmonic oscilator equation
using k/m=1. The input is the step
size h, and the methods used are the explicit
Euler method and midpoint method. 
h = Starting Step size (time)
The number of subintervals is chosen
so that three cycles are plotted. The error at each
step is then plotted for each method
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
    error = [anx[i]-xlist[i] for i in range(N+1)]
    return error
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
    error = [anx[i]-xlist[i] for i in range(N+1)]
    return error
if len(sys.argv)==2:
    [name,h] = sys.argv
    h=abs(float(h))
    N = int((3*2*pi)/h)
    t = [h*n for n in range(N+1)]
    listeuler = eulerrun(h)
    listmidpoint = midpointrun(h)
    # eulerplot = plot(t,listeuler,'b',label='Euler')
    midpointplot = plot(t,listmidpoint,'r',label='Midpoint')
    xlabel('Time/s')
    ylabel('Error')
    legend()
    show()
else:
    print 'Incorrect number of arguments'
