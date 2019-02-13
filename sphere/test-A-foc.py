# test-A-foc.py

from math import sin, pi
from numpy import zeros

def test_f(x):
    if x >= 2.0 and x <= 3.0:
        return sin(2.0*pi*x)
    else:
        return 0.0

def update_in_time(f, a, L, ncells, dt, nsteps):
    dx = float(L)/float(ncells)
    xs = zeros(ncells)
    us = zeros(ncells)
    fs = zeros(ncells+1)

    # Initialise positions and starting values
    for i in range(ncells):
        xs[i] = (i+0.5)*dx
        us[i] = f(xs[i])

    # Now begin timestepping
    for n in range(nsteps):
        # 1. Compute fluxes.
        # 1a. At all INTERNAL interfaces
        for i in range(1,ncells):
            fs[i] = (a/2.0)*(us[i-1] + us[i])
        # 1b. At the BOUNDARY interfaces
        fs[0] = 0.0; fs[-1] = a*us[-1]

        # 2. Compute new u value
        for i in range(ncells):
            us[i] = us[i] - dt*(1.0/dx)*(fs[i+1] - fs[i])
        
    return xs, us

from pylab import plot, show
dt = 0.01
nsteps = 50
L = 5.0
a = 1.0
ncells = 100
dx = L/ncells

# generate analytical solution
xa = zeros(ncells)
ua = zeros(ncells)
T = dt*nsteps
for i in range(ncells):
    xa[i] = (i+0.5)*dx
    ua[i] = test_f(xa[i] - a*T)

xn, un = update_in_time(test_f, a, L, ncells, dt, nsteps)

plot(xa, ua, "-", xn, un, "o")
show()



