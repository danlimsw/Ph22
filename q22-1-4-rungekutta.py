from math import *
from numpy import *
from pylab import *
from MyVector import *
"""
Runge-Kutta ODE Solver (maximum 10 variables!)
Input:
func = List of N strings, each representing the first derivative
    of the ith variable in terms of the variables x1...xN
init = List of N floats, initial conditions for x1 ... xN
t0 = Initial time (float)
h = Step size (seconds)
steps = Number of steps to evaluate
Output:
CSV file of t,x1,...,xN values for each step
"""
def rungekutta(vecold,t,h,func):
    # algorithm routine
    """
    Runge-Kutta step function
    vecold = N-vector of variable values [x1,x2,...,xN] at time t 
    t = initial time (float number)
    h = Time step size
    func = List of N strings that for the function of derivatives
    dxi/dt from independent variable t in terms of variables x1, x2 ... xN.
    Output = vecnew, N-list of parameter values at t+h.
    """
    if len(vecold) != len(func):
        return 'Length Error: Variable vector length: %d, Derivative function \
length: %d' %(len(vecold),len(func))
    N = len(vecold)
    t = float(t)
    h = float(h)
    vecold = [float(vecold[i]) for i in range(N)]
    prek1 = funceval(func,vecold,t)
    k1 = [h*prek1[i] for i in range(N)]
    prek2 = funceval(func,[vecold[i]+k1[i]/2.0 for i in range(N)],t+h/2.0)
    k2 = [h*prek2[i] for i in range(N)]
    prek3 = funceval(func,[vecold[i]+k2[i]/2.0 for i in range(N)],t+h/2.0)
    k3 = [h*prek3[i] for i in range(N)]
    prek4 = funceval(func,[vecold[i]+k3[i] for i in range(N)],t+h)
    k4 = [h*prek4[i] for i in range(N)]
    vecnew = [vecold[i]+k1[i]/6.0+k2[i]/3.0+k3[i]/3.0+k4[i]/6.0 for i in range(N)]
    return vecnew
def funceval(func,values,t):
    """
    Evaluates the vector-valued function func in terms
    of the variables x1, x2, ... , xN at values determined by
    the N-list values and at time t.
    """
    if len(func) != len(values):
        return 'Length Error: Function vector length: %d, \
Values vector length: %d' %(len(func),len(values))
    if len(values) > 10:
        return 'Too many variables'
    N = len(func)
    t = float(t)
    values = [float(values[i]) for i in range(N)]
    values = values + [0 for i in range(10-N)]
    values = values + [t]
    outlist = []
    for i in range(N):
        function = lambda x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,t:eval(func[i])
        out = function(*values)
        outlist = outlist + [out]
    return outlist
def rungeODE(func,init,t0,h,steps):
    # driver routine
    steps = abs(int(steps))
    t = float(t0)
    h = float(h)
    if len(func) != len(init):
        return 'Length error: Function vector length: %d, \
Initial condition value length: %d' %(len(func),len(init))
    if len(func) > 10:
        return 'Too many variables'
    N = len(func)
    init = [float(init[i]) for i in range(N)]
    savefile = open('RungeKuttaoutput.csv','w')
    savefile.write('t,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10\n')
    init = init + [0 for i in range(10-N)]
    func = func + ['0' for i in range(10-N)]
    # highly inefficient print method while I find a way to unpack in Python 2
    savefile.write('%f,'%t)
    for i in range(9):
        savefile.write('%f,'%init[i])
    savefile.write('%f\n'%init[9])
    xlist = [init[0]]
    ylist = [init[1]]
    for i in range(steps):
        output = rungekutta(init,t,h,func)
        t = t + h
        xlist = xlist + [output[0]]
        ylist = ylist + [output[1]]
        savefile.write('%f,'%t)
        for i in range(9):
            savefile.write('%f,'%init[i])
        savefile.write('%f\n'%init[9])
        init = output
    savefile.close()
    tlist = [float(t0)+i*h for i in range(steps+1)]
    subplot(311)
    plot(tlist,xlist)
    xlabel('Time/s')
    ylabel('x1')
    subplot(312)
    plot(tlist,ylist)
    xlabel('Time/s')
    ylabel('x2')
    subplot(313)
    plot(xlist,ylist)
    xlabel('x1')
    ylabel('x2')
    show()
