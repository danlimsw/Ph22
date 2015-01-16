"""
Root Finder
Author: Soon Wei Daniel Lim
Updated: 19 Dec 2014
"""
from math import *
def bisection(f,a,b,p):
    # Bisection Root Finder
    # f = Function of x, string
    # a = Guess, lower bracket
    # b = Guess, upper bracket
    # p = Precision (desired value of b-a)
    # Output = [a0,b0], bracket containing root
    a = float(a)
    b = float(b)
    fx = lambda x: eval(f)
    p = abs(float(p))
    # Error checks
    print 'f(a)=%f, f(b)=%f' %(fx(a), fx(b))
    if b <= a:
        return 'Error: choose proper bracket ordering.'
    if fx(a) == 0:
        if fx(b) == 0:
            return '%f and %f are roots' % (a,b)
        else:
            return '%f is a root.' % a
    elif fx(b) == 0:
        return '%f is a root.' % b
    elif fx(b)*fx(a) > 0:
        return 'Error: guesses have same sign'
    N = 0
    plist = []
    while b-a >= p:
        if fx((a+b)/2)*fx(a)>0: #midpoint same sign as f(a)
            a = (a+b)/2
            N = N + 1
            plist = plist + [b-a]
            #print '%d iterations, precision %e' %(N,b-a)
        elif fx((a+b)/2)*fx(b)>0: #midpoint same sign as f(b)
            b = (a+b)/2
            N = N + 1
            plist = plist + [b-a]
            #print '%d iterations, precision %e' %(N,b-a)
        else: #midpoint is a zero
            print type((a+b)/2.0)
            m = (a+b)/2.0
            return '%f is a root.' % m
    print plist
    return 'The root is in the interior of [%f,%f] with %d iterations' %(a,b,N)
def newton(f,df,a,t):
    # Newton-Raphson Method
    # f = Function of x, string
    # df = Derivative of f, string
    # a = Initial guess
    # t = Tolerance
    savefile = open('Newton-Rapshson-Output.csv','a')
    a = float(a)
    a0 = a
    t = abs(float(t))
    fx = lambda x: eval(f)
    dfx = lambda x: eval(df)
    if fx(a) == 0:
        return '%f is a root' % a
    N = 0
    if dfx(a) == 0:
        return 'Error: local extrema nearby'
    p0 = abs(float(fx(a))/float(dfx(a))) #initial precision (to normalize later)
    savefile.write('%s,%s,%f,%d,%f,%e\n'%(f,df,a0,0,a,1))
    while abs(fx(a)) >= t:
        a = a - fx(a)/dfx(a)
        if dfx(a) == 0:
            savefile.close()
            return 'Error: local extrema nearby'
        N = N + 1
        p = abs(fx(a)/dfx(a))/p0
        savefile.write('%s,%s,%f,%d,%f,%e\n'%(f,df,a0,N,a,p))
    p = abs(fx(a)/dfx(a))
    savefile.close()
    return 'The root is near %f with precision %e after %d iterations' %(a,p,N)
def secant(f,a,b,t):
    # Secant method
    # f = Function of x, string
    # a = Initial guess 1
    # b = Initial guess 2
    # t = Tolerance
    a = float(a)
    b = float(b)
    a0 = a
    b0 = b
    t = abs(float(t))
    t0 = t
    fx = lambda x: eval(f)
    if fx(a) == 0:
        if fx(b) == 0:
            return '%f and %f are roots' %(a,b)
        else:
            return '%f is a root.' %a
    if fx(b) == 0:
        return '%f is a root.' %b
    if (a == b) | (fx(a)==fx(b)):
        return 'Error: Equal guess or value at guess'
    df = (b - a)/(fx(b) - fx(a))
    p0 = abs(df*fx(b)) #initial precision for normalization
    p = p0
    N = 1
    savefile = open('Secant-Output.csv','a')
    c = b - fx(b)*df
    savefile.write('%s,%f,%f,%e,%d,%f,%e\n' %(f,a0,b0,t0,N,c,1))
    while (abs(fx(c)) >= t) & (b != c) & (fx(b) != fx(c)):
        N = N + 1
        a = b
        b = c
        df = (b - a)/(fx(b) - fx(a))
        c = b - fx(b)*df
        p = abs(df*fx(b))
        savefile.write('%s,%f,%f,%e,%d,%f,%e\n' %(f,a0,b0,t0,N,c,p/p0))
    savefile.close()
    return 'The root is near %f with precision %e after %d iterations' %(c,p,N)
def newtonsimple(f,df,a,t):
    # Newton-Raphson Method
    # f = Function of x, string
    # df = Derivative of f, string
    # a = Initial guess
    # t = Tolerance
    a = float(a)
    t = abs(float(t))
    fx = lambda x: eval(f)
    dfx = lambda x: eval(df)
    if fx(a) == 0:
        return a # a is a root
    N = 0
    if dfx(a) == 0:
        return 'Error: local extrema nearby'
    while abs(fx(a)) >= t:
        a = a - fx(a)/dfx(a)
        if dfx(a) == 0:
            return 'Error: local extrema nearby'
        N = N + 1
    return a
    
    
    
    
