# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 14:00:55 2020

@author: massonnetf
"""
import numpy as np
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