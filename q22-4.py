from math import *
from numpy import *
from pylab import *
import csv
from random import *
from time import *
import matplotlib.animation as animation
"""
Ph22.4 N-body Simulations with Approximations
Author: Soon Wei Daniel Lim
Last Updated: 14 Jan 2015
"""
class Quad(object):
    def __init__(self,bound=[[0,1],[1,1],[0,0],[1,0]],theta=0.1):
        """
        Bound = [TL,TR,BL,BR]
        """
        self.bound = bound # coordinates of corners
        self.length = abs(self.bound[0][1]-self.bound[2][1]) # length of side
        self.position = [None,None] # weighted average of contained particles
        self.children = [None,None,None,None]
        self.mass = 0.0
        self.theta = theta
    def split4(self):
        if self.children[0] != None: return # already split
        point1 = self.bound[0]
        point2 = [(self.bound[1][0]+self.bound[0][0])/2.0,self.bound[0][1]]
        point3 = self.bound[1]
        point4 = [self.bound[0][0],(self.bound[0][1]+self.bound[2][1])/2.0]
        point5 = [(self.bound[1][0]+self.bound[0][0])/2.0,(self.bound[0][1]+self.bound[2][1])/2.0]
        point6 = [self.bound[1][0],(self.bound[0][1]+self.bound[2][1])/2.0]
        point7 = self.bound[2]
        point8 = [(self.bound[1][0]+self.bound[0][0])/2.0,self.bound[2][1]]
        point9 = self.bound[3]
        self.children = [Quad([point1,point2,point4,point5],self.theta),Quad([point2,point3,point5,point6],self.theta),Quad([point4,point5,point7,point8],self.theta),Quad([point5,point6,point8,point9],self.theta)]
        # print 'Split Success'
    def insertin(self,position):
        if self.position[0] == None: # quad is currently empty, so put one particle
            # print 'Found empty quad'
            self.position = position
            self.mass = 1.0
        elif self.children[0] == None: # quad has exactly one particle
            # print 'Quad has one particle already'
            if self.position == position: return 0 # no repeated particles
            self.split4()
            self.mass = 2.0
            self.choosequadrant(self.position) # fit old particle
            self.choosequadrant(position) # fit new particle
            self.position = [(self.position[0]+position[0])/2.0,(self.position[1]+position[1])/2.0] # average of old pos and added pos
        else: # quad has already been split
            # print 'Entered split quad'
            self.mass = self.mass + 1.0
            self.choosequadrant(position)
            self.position = [(self.mass*self.position[0]+position[0])/(1.0+self.mass),(self.mass*self.position[1]+position[1])/(1.0+self.mass)] # weighted average
    def choosequadrant(self,position):
        if self.children[0] == None: return # not split yet
        if position[0] <= (self.bound[1][0]+self.bound[0][0])/2.0:
            if position[1] <= (self.bound[0][1]+self.bound[2][1])/2.0:
                self.children[2].insertin(position)
                # print 'Inserted in quadrant 3'
            else:
                self.children[0].insertin(position)
                # print 'Inserted in quadrant 1'
        else:
            if position[1] <= (self.bound[0][1]+self.bound[2][1])/2.0:
                self.children[3].insertin(position)
                # print 'Inserted in quadrant 4'
            else:
                self.children[1].insertin(position)
                # print 'Inserted in quadrant 2'
##    def mass(self):
##        if self.position[1] == None and self.children[1] == None: # no position nor children = empty
##            return 0
##        if self.position[1] != None and self.children[1] == None: # only one particle
##            return 1
##        else:
##            return sum([self.children[i].mass for i in range(4)])
    def acc(self,position,a=0.1,m=1):
        theta = self.theta
        if self.mass == 0.0:
            return [0.0,0.0]
        elif self.mass == 1.0:
            if self.isin(position): # particle is in the current quad
                return [0.0,0.0]
            else: # some other particle is in the quad
                return [m**2*(self.position[0]-position[0])/(d(self.position,position)*(d(self.position,position)**2+a**2)),m**2*(self.position[1]-position[1])/(d(self.position,position)*(d(self.position,position)**2+a**2))]
        else: # there is more than 1 particle in current quad (definitely have children)
            if self.isin(position): # particle is in the current quad, check children
                allacc = [self.children[i].acc(position,theta,a) for i in range(4)]
                return [sum([allacc[i][0] for i in range(4)]),sum([allacc[i][1] for i in range(4)])]
            else: # particle is outside. Now compare to theta.
                if self.length/d(self.position,position) < theta: # far enough away
                    return [m**2*self.mass*(self.position[0]-position[0])/(d(self.position,position)*(d(self.position,position)**2+a**2)),m**2*self.mass*(self.position[1]-position[1])/(d(self.position,position)*(d(self.position,position)**2+a**2))]
                else: # still too close, check children
                    allacc = [self.children[i].acc(position,theta,a) for i in range(4)]
                    return [sum([allacc[i][0] for i in range(4)]),sum([allacc[i][1] for i in range(4)])]
    def isin(self,position):
        if (position[0] <= self.bound[1][0]) and (position[0] > self.bound[0][0]) and (position[1] <= self.bound[0][1]) and (position[1] > self.bound[2][1]):
            return True
        else: return False
def d(a,b):
    return sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)
def nstep(data,dt,a=0.1,theta=0.1,m=1):
    """
    Perform the N-body simulation using the Barnes-Hut Algorithm
    Input:
    data = Nested list of depth 2, each element being [x,y,vx,vy]
    dt = Time step
    a = Force Softening Constant
    theta = Barnes-Hut opening angle
    Output:
    new = Updated data set
    """
    # start = time()
    N = len(data)
    a = float(a)
    dt = float(dt)
    change = [] # stores all change data over one timestep
    # generate Quad tree
    xmin = min([data[i][0] for i in range(N)])
    xmax = max([data[i][0] for i in range(N)])
    ymin = min([data[i][1] for i in range(N)])
    ymax = max([data[i][1] for i in range(N)])
    xlen = xmax - xmin
    ylen = ymax - ymin
    side = 1.1*max(xlen,ylen) # ensures all particles are inside, not on boundary
    xmin = xmin - 0.05*side
    ymin = ymin - 0.05*side
    root = Quad([[xmin,ymin+side],[xmin+side,ymin+side],[xmin,ymin],[xmin+side,ymin]],theta)
    # populate quad tree
    # end1 = time()
    # print 'Made tree', end1-start
    for i in range(N):
        root.insertin([data[i][0],data[i][1]])
    #end2 = time()
    # print 'Inserted elements', end2-end1
    # calculate changes over one timestep
    for i in range(N):
        acc = root.acc([data[i][0],data[i][1]],a,m)
        dvx = acc[0]*dt
        dvy = acc[1]*dt
        dx = (data[i][2]+dvx)*dt # symplectic Euler
        dy = (data[i][3]+dvy)*dt # symplectic Euler
        change = change + [[dx,dy,dvx,dvy]]
    # end3 = time()
    # print 'Calculated changes',end3-end2
    new = [[data[i][j]+change[i][j] for j in range(4)] for i in range(N)]
    # end4 = time()
    # print 'Included changes',end4-end3
    return new
def runqn1(N,dt,interv=1,a=0.1,theta=0.1,m=1):
    """
    Generates N particles randomly scattered around the unit square.
    Initial velocity is set to zero for all particles.
    Animates the resultant gravitational interaction via Barnes-Hut
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
    ax.set_ylim(0,1)
    ax.set_xlim(0,1)
    def datagen():
        data = datagen.data
        dt = datagen.dt
        a = datagen.a
        theta = datagen.theta
        m = datagen.m
        yield data,dt,a,theta,m
    datagen.data = data
    datagen.dt = dt
    datagen.a = a
    datagen.theta = theta
    datagen.m = m
    def run(inp):
        data,dt,a,theta,m = inp
        new = nstep(data,dt,a,theta,m)
        N = len(data)
        xdata = [new[i][0] for i in range(N)]
        ydata = [new[i][1] for i in range(N)]
        line.set_xdata(xdata)
        line.set_ydata(ydata)
        datagen.data = new
        return line,
    ani = animation.FuncAnimation(fig,run,datagen,interval=interv,blit=False)
    show()
def runqn2(N,dt,interv=1,a=0.1,v=0,theta=0.1,m=1):
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
        theta = datagen.theta
        m = datagen.m
        yield data,dt,a,theta,m
    datagen.data = data
    datagen.dt = dt
    datagen.a = a
    datagen.theta = theta
    datagen.m = m
    def run(inp):
        data,dt,a,theta,m = inp
        new = nstep(data,dt,a,theta,m)
        N = len(data)
        xdata = [new[i][0] for i in range(N)]
        ydata = [new[i][1] for i in range(N)]
        line.set_xdata(xdata)
        line.set_ydata(ydata)
        datagen.data = new
        return line,
    ani = animation.FuncAnimation(fig,run,datagen,interval=interv,blit=False)
    show()
def runtime(N,dt,theta,samples):
    """
    Measures the average runtime for the
    slow step (Barnes-Hut routine) using
    N particles scattered on the unit circle.
    Input:
    N = Number of particles
    dt = Time-step
    theta = Opening angle
    samples = Number of samples to take runtime for
    Output:
    Mean = Average runtime for all samples
    stderr = Standard Error of the mean
    """
    N = int(N)
    dt = float(dt)
    theta = float(theta)
    samples = int(samples)
    data = []
    v = 0
    for i in range(N):
        theta1 = random()*2*pi #fixes initial position
        theta2 = random()*2*pi #fixes initial velocity direction
        r = random() # fixes initial position (radial)
        data = data +[[r*cos(theta1),r*sin(theta1),v*cos(theta2),v*sin(theta2)]] # positions are random in the unit circle with initial velocity v
    start = time()
    output = []
    for i in range(samples):
        data = nstep(data,dt,0.1,theta,1)
        end = time()
        output = output + [end-start]
        start = time()
    mean = average(output)
    stderr = std(output)/sqrt(samples)
    return mean,stderr
def comparetime(N,dt,thetalist,samples):
    """
    Performs runtime(N,dt,theta,samples) for
    each theta in thetalist, then plots the values
    on graph.
    Input:
    N = Number of particles
    dt = Time-step
    thetalist = List of Opening angles to test
    samples = Number of samples for each run
    Output:
    Plot of Average Runtime against Theta
    """
    N = int(N)
    dt = float(dt)
    samples = int(samples)
    output = []
    outputerr = []
    for theta in thetalist:
        run = runtime(N,dt,theta,samples)
        output = output + [run[0]]
        outputerr = outputerr + [run[1]]
    errorbar(thetalist,output,outputerr)
    xlabel('Theta')
    ylabel('Average Runtime/s')
    show()
def energy(data):
    """
    Calculates the total kinetic and potential
    energy of a nested list of elements
    [x,y,vx,vy]
    Input:
    data = Nested list of [x,y,vx,vy]
    Output:
    energy = Total energy
    """
    N = len(data)
    KE = 0.5*sum([data[i][2]**2+data[i][3]**2 for i in range(N)])
    PE = sum(sum([[-1.0/(d(data[i][0:2],data[j][0:2])) for i in range(j)] for j in range(N)]))
    return KE+PE

def runenergy(N,dt,a,theta,steps):
    """
    Performs the Barnes-Hut algorithm and
    calculates the total energy at each step.
    Plots the total energy as a function of time.
    Input: N = number of particles
    dt = Time step
    a = Force softening constant
    theta = opening angle
    steps = number of timesteps
    Output:
    Plot of total energy against time.
    """
    seed(138)
    N = int(N)
    dt = float(dt)
    data = []
    steps = int(steps)
    v = 0
    for i in range(N):
        theta1 = random()*2*pi #fixes initial position
        theta2 = random()*2*pi #fixes initial velocity direction
        r = random() # fixes initial position (radial)
        data = data +[[r*cos(theta1),r*sin(theta1),v*cos(theta2),v*sin(theta2)]]
    output = [energy(data)]
    for i in range(steps-1):
        data = nstep(data,dt,a,theta,1)
        output = output + [energy(data)]
    tlist = [dt*i for i in range(steps)]
    plot(tlist,output)
    xlabel('Time')
    ylabel('Energy')
    ylim([-900,-700])
    show()
def runenergy2(N,dt,a,thetalist,steps):
    """
    Performs the Barnes-Hut Algorithm and
    calculates the total energy at each step.
    Repeats the calculation for each theta
    in thetalist (using the same data), then plots
    the fractional deviation
    in energy at the end of each simulation
    against theta.
    Input:
    N = number of particles
    dt = Timestep
    a = Force softening constant
    thetalist = list of theta values to consider
    steps = number of timesteps
    Output:
    Plot of fractional deviation of energy
    (final energy - initial energy)/initial energy
    against theta
    """
    seed(138)
    N = int(N)
    dt = float(dt)
    data = []
    steps = int(steps)
    v = 0
    for i in range(N):
        theta1 = random()*2*pi #fixes initial position
        theta2 = random()*2*pi #fixes initial velocity direction
        r = random() # fixes initial position (radial)
        data = data +[[r*cos(theta1),r*sin(theta1),v*cos(theta2),v*sin(theta2)]]
    startdata = data
    initial = energy(data)
    output = []
    for j in range(len(thetalist)):
        theta = thetalist[j]
        data = startdata
        for i in range(steps-1):
            data = nstep(data,dt,a,theta,1)
        output = output + [(energy(data)-initial)/initial]
    plot(thetalist,output)
    xlabel('Theta')
    ylabel('Fractional Energy Deviation')
    show()
