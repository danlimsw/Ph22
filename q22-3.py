from math import *
from numpy import *
from pylab import *
import csv
from random import *
import matplotlib.animation as animation
"""
Ph22.3 N-body Simulations
Author: Soon Wei Daniel Lim
Last Updated: 12 Jan 2015
"""
def nstep(data,dt,a):
    """
    Performs an N-body symplectic Euler integrator update
    to the data input with timestep dt.
    Assumptions: G=1, m=1
    Input:
    data = Nested list of depth 2. Each element of the list
    is a 4-member list of a particle's [x,y,vx,vy]
    dt = Time step (arbitrary units)
    a = Force Softening constant
    Output:
    new = updated data set
    """
    # convert all data to floats
    old = [[float(data[i][j]) for j in range(4)] for i in range(len(data))]
    dt = float(dt)
    N = len(data)
    a = float(a)
    change = [] # stores all change data over one timestep.
    for i in range(N):
        therest = range(N)
        therest.pop(i) # remove yourself from iterable list
        dx3 = sum([dt*(old[j][0]-old[i][0])/(((old[j][0]-old[i][0])**2+(old[j][1]-old[i][1])**2+a**2)**1.5) for j in therest])
        dx1 = (old[i][2]+dx3)*dt
        dx4 = sum([dt*(old[j][1]-old[i][1])/(((old[j][0]-old[i][0])**2+(old[j][1]-old[i][1])**2+a**2)**1.5) for j in therest])
        dx2 = (old[i][3]+dx4)*dt
        change = change + [[dx1,dx2,dx3,dx4]]
    new = [[old[i][j]+change[i][j] for j in range(4)] for i in range(N)]
    return new
def runqn1(N,dt,interv,a):
    """
    Generates N particles randomly scattered around the unit square.
    Initial velocity is set to zero for all particles.
    Animates the resultant gravitational interaction.
    Input:
    N = number of particles
    dt = Timestep
    Interval = Milliseconds between steps
    a = Force softening constant
    Output:
    Animation of the N particles.
    """
    N = int(N)
    dt = float(dt)
    interv = int(interv)
    data = [[random(),random(),0,0] for i in range(N)] # positions are random in the unit square
    fig,ax = subplots()
    xdata = [data[i][0] for i in range(N)]
    ydata = [data[i][1] for i in range(N)]
    line, = ax.plot(xdata,ydata,'.')
    ax.set_ylim(-1,1)
    ax.set_xlim(-1,1)
    def datagen():
        data = datagen.data
        dt = datagen.dt
        a = datagen.a
        yield data,dt,a
    datagen.data = data
    datagen.dt = dt
    datagen.a = a
    def run(inp):
        data,dt,a = inp
        new = nstep(data,dt,a)
        N = len(data)
        xdata = [new[i][0] for i in range(N)]
        ydata = [new[i][1] for i in range(N)]
        line.set_xdata(xdata)
        line.set_ydata(ydata)
        datagen.data = new
        return line,
    ani = animation.FuncAnimation(fig,run,datagen,interval=interv,blit=False)
    show()
def runqn2(N,dt,interv,a,v):
    """
    Does exactly the same thing as runqn1,
    but places the particles in a unit sphere
    and sets the initial velocity as v (random directions)
    Input:
    N = number of particles
    dt = Timestep
    Interval = Milliseconds between steps
    a = Force softening constant
    v = Initial velocity
    Output:
    Animation of the N particles.
    """
    N = int(N)
    dt = float(dt)
    interv = int(interv)
    v = float(v)
    data = []
    for i in range(N):
        theta1 = random()*2*pi #fixes initial position
        theta2 = random()*2*pi #fixes initial velocity direction
        r = random() # fixes initial position (radial)
        data = data +[[r*cos(theta1),r*sin(theta1),v*cos(theta2),v*sin(theta2)]] # positions are random in the unit circle with initial velocity v
    fig,ax = subplots()
    xdata = [data[i][0] for i in range(N)]
    ydata = [data[i][1] for i in range(N)]
    line, = ax.plot(xdata,ydata,'.')
    ax.set_ylim(-1,1)
    ax.set_xlim(-1,1)
    def datagen():
        data = datagen.data
        dt = datagen.dt
        a = datagen.a
        yield data,dt,a
    datagen.data = data
    datagen.dt = dt
    datagen.a = a
    def run(inp):
        data,dt,a = inp
        new = nstep(data,dt,a)
        N = len(data)
        xdata = [new[i][0] for i in range(N)]
        ydata = [new[i][1] for i in range(N)]
        line.set_xdata(xdata)
        line.set_ydata(ydata)
        datagen.data = new
        return line,
    ani = animation.FuncAnimation(fig,run,datagen,interval=interv,blit=False)
    show()
