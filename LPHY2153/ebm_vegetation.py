# Climate-vegetation interactions with a conceptual model
#
# Ref:    Rombouts and Ghil, 2015: doi:10.5194/npg-22-275-2015
# Author: F. Massonnet
# Date:   August 2018

# Imports, closing figures.
import numpy as np
import matplotlib.pyplot as plt
plt.close("all")

# Model parameters
# ----------------
C_T       = 500.0   # Heat capacity of the system [W yr K^-1 m^-2]
Q_0       = 342.5   # Incoming solar radiation    [W m^-2]
gamma     = 0.1    # Death rate of vegetation    [yr^-1]
dt        = 1.0     # Model time step             [years]
tf        = 5000.0  # Integration length          [years]
p         = 0.3     # Fraction of land            [-]
alpha_max = 0.85    # Maximum albedo for ocean 
                    # (ice-covered)               [-]
alpha_min = 0.25    # Minimum albedo for ocean 
                    # (ice-free)                  [-]
T_alphal  = 263.0   # Threshold temperature below
                    # which ocean is ice-covered  [K]
T_alphau  = 300.0   # Threshold temperature above 
                    # which ocean is ice-free     [K]
alpha_v   = 0.1     # Albedo of vegetation        [-]
alpha_g   = 0.4     # Albedo of ground            [-]
B_0       = 200.0   # Intercept linearization of 
                    # outgoing longwave radiation [W m^-2]
B_1       = 2.5     # Slope linearization of out-
                    # going longwave radition     [W K^-1 m^-2]
T_opt     = 283.0   # Optimal growth temperature  [K]
k         = 0.004   # Parameter for width of 
                    # growth curve                [yr^-1 K^-2]
eps = 1e-12         # Threshold below which vegetation is set
                    # to 0.0 [-]

# Initial conditions
# ------------------
T_0 = 297.0  # [K]
A_0 = 0.8    # [-]

# Initialization
# --------------
time = np.arange(0, tf, dt)
nt = len(time)

T = np.empty(nt)
T[:] = np.nan

A = np.empty(nt)
A[:] = np.nan

T[0] = T_0
A[0] = A_0

# External forcing
# ----------------
# A ramp-up / ramp-down shape function is taken between
# yearb = 1850 and yearf = 2350 (mid-point: 2100)
# Set "rcp" to 0.0 for no forcing.

rcp = 0.0                 # [Wm^-2]
yearb, yearf = 1850, 2350

forcing = np.zeros(nt)  # Time series of forcing [W m^-2]

# Linearly increasing part of the forcing
tmp1 = (time - yearb) * rcp / ((yearf - yearb) / 2.0)
tmp1[tmp1 < 0.0] = 0.0

tmp2 = rcp - rcp / ((yearf - yearb) / 2.0) * (time - ((yearf + yearb) / 2.0))
tmp2[tmp2 < 0.0] = 0.0

forcing = np.minimum(tmp1, tmp2)

# Plotting the forcing
plt.figure(figsize = (4, 3))
plt.plot(time, forcing, lw = 3)
plt.ylabel("W/m$^2$")
plt.grid()
plt.tight_layout()
plt.savefig("./forc.png", dpi = 300)

# Specific model functions
# ------------------------

# Ocean albedo
def alpha_o(T):
    
    albedo = (T <= T_alphal) * alpha_max + (T > T_alphau) *  \
             alpha_min + (T <= T_alphau) * (T > T_alphal) *  \
             (alpha_max + (alpha_min - alpha_max) / (T_alphau - T_alphal) * \
             (T - T_alphal))  
            
    return albedo

# Surface albedo
def alpha(T, A):
  
    albedo = (1.0 - p) * alpha_o(T) + \
             p         * (alpha_v * A + alpha_g * (1.0 - A))

    return albedo

# Outgoing longwave radiation
def R_o(T):
    
    outgoing = B_0 + B_1 * (T - T_opt)
    return outgoing

# Coefficient to logistic growth term
def beta(T):
    
    tmp = 1.0 - k * (T - T_opt) ** 2
    
    return (tmp > 0) * tmp

# Jacobian
# --------
def J(T, A):
    if T <= T_alphal or T > T_alphau:
        dalbdoT = 0.0
    else:
        dalbdoT = (alpha_min - alpha_max) / (T_alphau - T_alphal)
    
    if beta(T) >= 0:
        dbetadT = - 2.0 * k * (T - T_opt)
    else:
        dbetadT = 0.0
    
    a = - 1.0 / C_T *(Q_0 * (1 - p) * dalbdoT + B_1)
    
    b = - 1.0 / C_T * Q_0 * p * (alpha_v - alpha_g)
    
    c = dbetadT * A * (1 - A)
    
    d = beta(T) * (1 - 2 * A) - gamma
    
    jacobian = np.array([[a, b], [c, d]])
    
    return jacobian

# Integration
# -----------
# Following Euler explicit scheme:
# dX / dt = f(X) --> X[t+dt] = X[t] + dt * f(X[t])

for jt in np.arange(0, nt - 1):
    
    # 1. Evaluation of the non linear right-hand-side term
        
    f_T = 1.0 / C_T * ((1.0 - alpha(T[jt], A[jt])) * Q_0 - R_o(T[jt]) + forcing[jt])
    f_A = beta(T[jt]) * A[jt] * (1.0 - A[jt]) - gamma * A[jt] 
    
    # 2. Update
    dT = f_T * dt
    dA = f_A * dt
    
    T[jt + 1] = T[jt] + dT
    A[jt + 1] = A[jt] + dA
    
    # 3. Regularize
    if A[jt + 1] < eps:
        A[jt + 1] = 0.0
        
        
# Visualization of the results
# ----------------------------
plt.figure(figsize = (6, 4))
plt.subplot(2, 1, 1)
plt.plot(time, T, lw = 2.0, color = "orange")
plt.ylabel("K")
plt.title("Temperature (T$_\infty$ = " + str(np.round(T[-1], 2)) + " K)")
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(time, A, lw = 2.0, color = "darkgreen")
plt.ylabel("-")
plt.title("Vegetation (A$_\infty$ = " + str(np.round(A[-1], 2)) + ")")
plt.grid()
plt.tight_layout()
plt.savefig("./T_A.png", dpi = 300)


# Plot vector field and trajectory in phase space
# -----------------------------------------------
TT = np.linspace(200, 350, 100)
AA = np.linspace(-0.1, 1.1, 20)
TT, AA = np.meshgrid(TT, AA)

# Vector field (with unit norm)
dTT = 1.0 / C_T * ((1.0 - alpha(TT, AA)) * Q_0 - R_o(TT))
dAA = beta(TT) * AA * (1.0 - AA) - gamma * AA
norm = np.sqrt(dTT ** 2 + dAA ** 2)
dTT = dTT / norm
dAA = dAA / norm

plt.figure(figsize = (12, 6))
plt.quiver(TT, AA, dTT, dAA, color = [0.5, 0.5, 0.5], scale = 50.0, angles = "xy")

# Plot the nullclines, i.e. the ensemble of (T, A) points
# that are fixed points for the T and A equations
# Note that the time-varying forcing term is ignored here

x = np.linspace(200, 350, 10000)
plt.plot(x, [np.max((0, 1 - gamma / beta(xx))) for xx in x], color = "blue")
plt.plot(x, [1.0 / (alpha_v - alpha_g) * (1.0 / p * (1.0 - R_o(xx) /  
         Q_0 - (1.0 - p) * alpha_o(xx)) - alpha_g) for xx in x], color = "red")
plt.scatter(T[0], A[0], color = "black")
plt.plot(T, A, color = "black", lw = 1)
plt.xlabel("Temperature T [K]")
plt.ylabel("Vegetation cover A [-]")
plt.xlim(295, 310)
plt.ylim(0.0, 1)
plt.grid()
plt.savefig("./phase.png", dpi = 300)