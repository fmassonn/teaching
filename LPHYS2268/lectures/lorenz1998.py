# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 11:49:56 2019

@author: massonnetf

Reference:
Lorenz, Edward N., and Kerry A. Emanuel. “Optimal Sites for Supplementary Weather Observations: Simulation with a Small Model.” Journal of the Atmospheric Sciences 55, no. 3 (February 1, 1998): 399–414. https://doi.org/10.1175/1520-0469(1998)055<0399:OSFSWO>2.0.CO;2.

"""
import numpy as np
import matplotlib.pyplot as plt

# Fourth-order Runge-Kutta
# From https://people.sc.fsu.edu/~jburkardt/py_src/rk4/rk4.py
#*****************************************************************************80
#
## RK4_TEST_F evaluates the right hand side of a particular ODE.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    18 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real T, the current time.
#
#    Input, real U, the current solution value.
#
#    Output, real VALUE, the value of the derivative, dU/dT.

def rk4vec ( t0, m, u0, dt, f ):

#*****************************************************************************80
#
## RK4VEC takes one Runge-Kutta step for a vector ODE.
#
#  Discussion:
#
#    Thanks  to Dante Bolatti for correcting the final function call to:
#      call f ( t1, m, u3, f3 )
#    18 August 2016.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    18 August 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, real T0, the current time.
#
#    Input, integer M, the spatial dimension.
#
#    Input, real U0(M), the solution estimate at the current time.
#
#    Input, real DT, the time step.
#
#    Input, function uprime = F ( t, m, u  ) 
#    which evaluates the derivative UPRIME(1:M) given the time T and
#    solution vector U(1:M).
#
#    Output, real U(M), the fourth-order Runge-Kutta solution 
#    estimate at time T0+DT.
#
#  Get four sample values of the derivative.
#
  f0 = f ( t0, m, u0 )

  t1 = t0 + dt / 2.0
  u1 = np.zeros ( m )
  u1[0:m] = u0[0:m] + dt * f0[0:m] / 2.0
  f1 = f ( t1, m, u1 )

  t2 = t0 + dt / 2.0
  u2 = np.zeros ( m )
  u2[0:m] = u0[0:m] + dt * f1[0:m] / 2.0
  f2 = f ( t2, m, u2 )

  t3 = t0 + dt
  u3 = np.zeros ( m )
  u3[0:m] = u0[0:m] + dt * f2[0:m]
  f3 = f ( t3, m, u3 )
#
#  Combine them to estimate the solution U at time T1.
#
  u = np.zeros ( m )
  u[0:m] = u0[0:m] + ( dt / 6.0 ) * ( \
            f0[0:m] \
    + 2.0 * f1[0:m] \
    + 2.0 * f2[0:m] \
    +       f3[0:m] )

  return u

# State vector dimension
M = 40

# Time and length of integration
# Following LE98 paper, one unit of time is 5 days

t0 = 0.0
# Time step, in units "5 days"
dt = 0.01 # This is 6 hr
# Integration time
N = int(((2.3 ) / (5 * dt)))
N = int(10 / (5 * dt))

t = np.array([t0 + n * dt for n in range(N)])

# Forcing 
F = 8.0
# Initial state: Lorenz and Emanuel initial conditions
X0 = np.full(M, F)
X0[20 - 1] = F + 0.008   

# Solution
X = np.full((M, N), np.nan)
X[:, 0] = X0

# Right-hand side of ODE
# X is a column np array of size M
def f(t, n, X):
    if len(X.shape) != 1:
        sys.exit("X not column")
    M = len(X)
    
    value  = np.full(M, np.nan)
    for jM in np.arange(M):
        jMp1 = (jM + 1   ) % M
        jMm2 = (jM + (-2)) % M
        jMm1 = (jM + (-1)) % M

        value[jM] = (X[jMp1] - X[jMm2]) * X[jMm1] - X[jM] + F
    
    return value

for jN in np.arange(1, N):
    # Runge-Kutta 4th order method
    X[:, jN] = rk4vec(t[jN - 1], M, X[:, jN - 1], dt, f)

        
# Reproduction of LE98 Fig. 1
#fig = plt.figure(figsize = (6, 6), dpi = 300)
#plt.grid()
## Max range
#rangemax = np.max(np.max(X, axis = 1) - np.min(X, axis = 1))
#for jN in range(N):
#    plt.plot(range(M), X[:, jN] -  rangemax * jN, "-")
#plt.yticks([F - rangemax * jN for jN in np.arange(0, N)], ["0", "", "", "", "1", "", "", "", "2"])
#plt.xticks(np.arange(4, J + 1, 5), [str(j) for j in np.arange(5, J + 2, 5)])
#plt.xlabel("site")
#plt.ylabel("time (days)")
#plt.savefig("LE98.png")

fig = plt.figure(figsize = (8, 3), dpi = 100)
# Plot solutions
[plt.plot(t, X[jM, :], lw = 0.1, color = "grey") for jM in range(J)]

# Plot one station
plt.plot(t, X[0, :], lw = 2, color = "orange")
# Plot ensemble mean
mean = np.mean(X, axis = 0)
plt.plot(t, mean , lw = 1, color = "k")
# Plot ensemble std
std = np.std(X, axis = 0)
plt.fill_between(t, mean - std, mean + std, color = [0.0, 0.0, 0.0], alpha = 0.1)

plt.xlabel("Days")

fig = plt.figure(figsize = (8, 3), dpi = 100)

plt.plot(X[0, :], X[1, :], lw = 2, color = "orange")
