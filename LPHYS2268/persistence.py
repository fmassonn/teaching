# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 13:45:17 2020

@author: massonnetf
"""
import numpy as np

# Seed
np.random.seed(1)

# Time axis
nt = 1000
t = np.arange(0, nt)


# AR1 model
# X_{t+1} = mu + alpha * X_t  + std * eps
mu = 0
alpha = 0.8 
std = 1.0
# Generation of data
x = np.full(nt, np.nan)

for jt in np.arange(nt):
    if jt == 0:
        x[jt] = std * np.random.randn()
    else:
        x[jt] = mu + alpha * x[jt - 1] + std * np.random.randn()
        

plt.figure(dpi = 300, figsize = (6, 3))        
plt.grid()
plt.plot(t, x, color = [0.5, 0.5, 0.5], lw = 3)
plt.title(r"AR1 process, $\alpha$ = " + str(alpha))
plt.xlim(0.0,200)
plt.xlabel("Time")
plt.tight_layout()
plt.savefig("./ar1.png")

for lag in np.array([0, 1, 3, 5, 10, 30]):
    plt.figure(figsize = (3, 3), dpi = 300)
    rho = np.corrcoef(x[:nt - lag], x[lag:])[0, 1]
    plt.scatter(x[:nt - lag], x[lag:], 5, color = [0.5, 0.5, 0.5])
    plt.title("Lag = " + str(lag) + "; r = " + str(np.round(rho, 2)))
    plt.xlim(-8, 8)
    plt.ylim(-8, 8)
    plt.grid()
    plt.savefig("./fig_lag" + str(lag).zfill(3) + ".png")


# Auto-correlation function
fig = plt.figure(figsize = (4, 3), dpi = 300)
for lag in np.arange(0, 50):
    if lag == 0:
        label = "Sample autocorrelation"
    else:
        label = None
    rho = np.corrcoef(x[:nt - lag], x[lag:])[0, 1]
    plt.scatter(lag, rho, 50, marker = "x", color = [0.0, 0.0, 1.0],
                label = label)

plt.plot(np.arange(0, 50), alpha ** np.arange(0, 50), color = [0.0, 1.0, 0.5], 
         linestyle = "-", label = "Theoretical value")

plt.xlim(-1.0, 20.0)
plt.plot((0.0, 100.0), (1 / np.exp(1), 1/ np.exp(1)), "r--", label = "1/e")
plt.xlabel(r"$\tau$")
plt.legend()
plt.gca().set_axisbelow(True)
plt.grid()
plt.tight_layout()
plt.savefig("./acf.png")
