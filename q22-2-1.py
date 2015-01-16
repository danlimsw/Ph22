"""
Problem Set 22-2 Three Body Problem
Author: Soon Wei Daniel Lim
Last Updated: 11 Jan 2015
"""
from rungekutta12 import *
from math import *
from numpy import *
from pylab import *
import csv
# physical variables
def runqn1(h,steps,alpha):
    """
    Numerically simulates the motion of an asteroid
    in the Jupiter-Sun corrotating frame.
    Input:
    h = step size (seconds)
    steps = number of steps
    alpha = initial angle with respect to the Jupiter-Sun line (rad)
    Output: Plot of asteroid position (y against x),
    CSV file saved as 'RungeKuttaOutput.csv' with all data.
    """
    G = 6.6742*10**-11 # m^3 kg^-1 s^-2
    M2 = 1.989*10**30  # kg, Sun
    M1 = 1.899*10**27 # kg, Jupiter
    R = 778.3*10**9 # m, Semimajor axis of Jupiter's orbit
    T = 3.743*10**8 # s, Period of Jupiter's orbit
    omega = 1.679*10**-8 # angular velocity of J-S system
    func = ['x3','x4','-6.6742*10**-11*1.899*10**27*(x1-1.989*10**30*778.3*10**9/(1.899*10**27+1.989*10**30))/(((x1-1.989*10**30*778.3*10**9/(1.899*10**27+1.989*10**30))**2+x2**2)**(3.0/2))-6.6742*10**-11*1.989*10**30*(x1+1.899*10**27*778.3*10**9/(1.899*10**27+1.989*10**30))/(((x1+1.899*10**27*778.3*10**9/(1.899*10**27+1.989*10**30))**2+x2**2)**(3.0/2))+2.0*1.679*10**-8*x4+x1*(1.679*10**-8)**2','-6.6742*10**-11*1.899*10**27*x2/(((x1-1.989*10**30*778.3*10**9/(1.899*10**27+1.989*10**30))**2+x2**2)**(3.0/2))-6.6742*10**-11*1.989*10**30*x2/(((x1+1.899*10**27*778.3*10**9/(1.899*10**27+1.989*10**30))**2+x2**2)**(3.0/2))-2.0*1.679*10**-8*x3+x2*(1.679*10**-8)**2']
    init = [cos(alpha)*R*(M2-M1)/(M2+M1),R*sin(alpha),0,0]
    t0 = 0 # initial time
    h = float(h)
    steps = abs(int(steps))
    rungeODE(func,init,t0,h,steps)
def runqn2(h,steps):
    """
    Numerically simulations the motion of three equal masses
    of mass 1 interacting under gravity with G=1, placed on
    the vertices of an equilateral triangle of side length =1.
    The initial conditions are chosen to ensure uniform circular
    motion.
    Input:
    h = step size (arbitrary unit)
    steps - number of steps
    Output: Plot of mass positions, output csv file 'RungeKuttaOutput.csv'
    """
    func = ['x7','x8','x9','x10','x11','x12','-(x1-x3)/(((x1-x3)**2+(x2-x4)**2)**(3.0/2))-(x1-x5)/(((x1-x5)**2+(x2-x6)**2)**(3.0/2))','-(x2-x4)/(((x1-x3)**2+(x2-x4)**2)**(3.0/2))-(x2-x6)/(((x1-x5)**2+(x2-x6)**2)**(3.0/2))','-(x3-x1)/(((x3-x1)**2+(x4-x2)**2)**(3.0/2))-(x3-x5)/(((x3-x5)**2+(x4-x6)**2)**(3.0/2))','-(x4-x2)/(((x3-x1)**2+(x4-x2)**2)**(3.0/2))-(x4-x6)/(((x3-x5)**2+(x4-x6)**2)**(3.0/2))','-(x5-x1)/(((x5-x1)**2+(x6-x2)**2)**(3.0/2))-(x5-x3)/(((x5-x3)**2+(x6-x4)**2)**(3.0/2))','-(x6-x2)/(((x5-x1)**2+(x6-x2)**2)**(3.0/2))-(x6-x4)/(((x5-x3)**2+(x6-x4)**2)**(3.0/2))']
    init = [0,0,-0.5,-sqrt(3)/2,0.5,-sqrt(3)/2,-1,0,0.5,-sqrt(3)/2,0.5,sqrt(3)/2]
    rungeODE(func,init,0,h,steps)
    # extract data from CSV output file
    filename = open('RungeKuttaOutput.csv','rb')
    reader = csv.reader(filename,delimiter=',')
    x1=[]
    x2=[]
    x3=[]
    x4=[]
    x5=[]
    x6=[]
    for row in reader:
        x1 = x1 + [row[1]]
        x2 = x2 + [row[2]]
        x3 = x3 + [row[3]]
        x4 = x4 + [row[4]]
        x5 = x5 + [row[5]]
        x6 = x6 + [row[6]]
    x1out=[float(x1[i]) for i in range(1,len(x1))]
    x2out=[float(x2[i]) for i in range(1,len(x1))]
    x3out=[float(x3[i]) for i in range(1,len(x1))]
    x4out=[float(x4[i]) for i in range(1,len(x1))]
    x5out=[float(x5[i]) for i in range(1,len(x1))]
    x6out=[float(x6[i]) for i in range(1,len(x1))]
    plot(x1out,x2out,'r',x3out,x4out,'g',x5out,x6out,'b')
    xlabel('X')
    ylabel('Y')
    show()
def runqn3(h,steps):
    """
    Numerically simulates the choreographic orbit with
    the given initial conditions in the worksheet.
    Input:
    h = step size (arbitrary units)
    steps = number of steps
    Output: Plot of mass positions, output csv file 'RungeKuttaOutput.csv'
    """
    func = ['x7','x8','x9','x10','x11','x12','-(x1-x3)/(((x1-x3)**2+(x2-x4)**2)**(3.0/2))-(x1-x5)/(((x1-x5)**2+(x2-x6)**2)**(3.0/2))','-(x2-x4)/(((x1-x3)**2+(x2-x4)**2)**(3.0/2))-(x2-x6)/(((x1-x5)**2+(x2-x6)**2)**(3.0/2))','-(x3-x1)/(((x3-x1)**2+(x4-x2)**2)**(3.0/2))-(x3-x5)/(((x3-x5)**2+(x4-x6)**2)**(3.0/2))','-(x4-x2)/(((x3-x1)**2+(x4-x2)**2)**(3.0/2))-(x4-x6)/(((x3-x5)**2+(x4-x6)**2)**(3.0/2))','-(x5-x1)/(((x5-x1)**2+(x6-x2)**2)**(3.0/2))-(x5-x3)/(((x5-x3)**2+(x6-x4)**2)**(3.0/2))','-(x6-x2)/(((x5-x1)**2+(x6-x2)**2)**(3.0/2))-(x6-x4)/(((x5-x3)**2+(x6-x4)**2)**(3.0/2))']
    init = [0,0,0.97000436,-0.24308753,-0.97000436,0.24308753,-0.93240737,-0.86473146,0.466203685,0.43236573,0.466203685,0.43236573]
    rungeODE(func,init,0,h,steps)
    # extract data from CSV output file
    filename = open('RungeKuttaOutput.csv','rb')
    reader = csv.reader(filename,delimiter=',')
    x1=[]
    x2=[]
    x3=[]
    x4=[]
    x5=[]
    x6=[]
    for row in reader:
        x1 = x1 + [row[1]]
        x2 = x2 + [row[2]]
        x3 = x3 + [row[3]]
        x4 = x4 + [row[4]]
        x5 = x5 + [row[5]]
        x6 = x6 + [row[6]]
    x1out=[float(x1[i]) for i in range(1,len(x1))]
    x2out=[float(x2[i]) for i in range(1,len(x1))]
    x3out=[float(x3[i]) for i in range(1,len(x1))]
    x4out=[float(x4[i]) for i in range(1,len(x1))]
    x5out=[float(x5[i]) for i in range(1,len(x1))]
    x6out=[float(x6[i]) for i in range(1,len(x1))]
    plot(x1out,x2out,'r',x3out,x4out,'g',x5out,x6out,'b')
    show()
