# -*- coding: utf-8 -*-
"""
Created on Mon May 13 22:59:45 2019

@author: massonnetf
"""
import numpy as np
import matplotlib.pyplot as plt
# Carnot cycle

gamma = 1.4 * 2 # For design purposes
R = 287.0

fig = plt.figure(figsize = (5, 5))
p_A = 300 * 100.0
p_B = 150 * 100.0
T_h = 600.0
T_c = 200.0

V = np.arange(0.001, 100.0, 0.01)
# Plot 4 points
# A
V_A = R * T_h / p_A
plt.scatter(V_A, p_A, 150, marker = ".", color = "k")
plt.text(V_A, p_A, "  A", ha = "left", fontsize = 20)
plt.xlim(0.0, 40.0)
plt.ylim(0.0, 40000)
plt.xlabel("V", fontsize = 20)
plt.ylabel("P", fontsize = 20)
plt.xticks([],[])
plt.yticks([],[])
#plt.axis("off")
plt.grid()
plt.tight_layout()
plt.savefig("Carnot1.png", dpi = 500)
# Plot isothermal passing at A
plt.plot(V, R * T_h / V, "r--")
plt.savefig("Carnot2.png", dpi = 500)

# A --> B: isothermal expansion. T_B = T_A (pV = constant = R * T_A)
# Point B
V_B = R * T_h / p_B
plt.scatter(V_B, p_B, 150, marker = ".", color = "k")
plt.text(V_B, p_B, " B", ha = "left", fontsize = 20)
plt.savefig("Carnot3.png", dpi = 500)

# B --> C: adiabatic expansion (pV^gamma = constant = p_B * V_B ^ gamma )
# adiabat at B
cst_B = p_B * V_B ** gamma
plt.plot(V, cst_B / V ** gamma, 150, color = "purple")
plt.savefig("Carnot4.png", dpi = 500)


# One must arrive at temperature T_c so that 
# * p_C * V_C = R * T_c
# * p_C * V_C^gamma = p_B * V_B ^gamma
# -> R * T_c * V_C ^(gamma - 1) = p_B * V_B ^gamma
# -> V_C = ...
V_C = (p_B * V_B ** gamma / (R * T_c)) ** (1.0 / (gamma - 1))
p_C = R * T_c / V_C

plt.scatter(V_C, p_C, 150, marker = ".", color = "k")
plt.text(V_C, p_C, " C", ha = "left", fontsize = 20)
plt.savefig("Carnot5.png", dpi = 500)

# Isothermal at C
plt.plot(V, R * T_c / V, "b--")
plt.savefig("Carnot6.png", dpi = 500)


# C --> D: isothermal compression
# Point D is on the adiabat passing at A and on the isotherm passing at C
# p_D * V_D ^gamma = cst_A
# p_D * V_D        = R * T_c
# -> R * T_c * V_D ^ (gamma - 1) = cst_A
# -> V_D = ...
V_D = (cst_A / (R * T_c)) ** (1.0 / (gamma - 1))
p_D = R * T_c / V_D

plt.scatter(V_D, p_D, 150, marker = ".", color = "k")
plt.text(V_D, p_D, "D", va = "top", ha = "left", fontsize = 20)
plt.savefig("Carnot7.png", dpi = 500)

# adiabat at A
cst_A = p_A * V_A ** gamma
plt.plot(V, cst_A / V ** gamma, 150, color = "pink")
plt.savefig("Carnot8.png", dpi = 500)






# Adiabatic
#, 