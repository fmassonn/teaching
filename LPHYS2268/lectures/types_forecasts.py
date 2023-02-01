# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 16:39:24 2020

@author: massonnetf
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

# Illustration of four types of forecasts
def normal(xx, mu = 0, sigma = 1):
    return 1.0 / np.sqrt(2 * np.pi * sigma ** 2) * np.exp(- ((xx - mu) / (np.sqrt(2) * sigma)) ** 2)

    

fig, axs = plt.subplots(2, 1, figsize = (5, 6), dpi = 300)
ax = axs[0]

ax.set_xlim(21.0, 25.0)
ax.set_ylim(0.0, 2.0)
ax.set_xlabel("Physical variable (e.g. temperature in °C)")
ax.set_ylabel("PDF")
ax.grid()
ax.set_axisbelow(True)
ax.set_title("Different forecast PDFs")

xx = np.linspace(-10, 30, 10000)

cc = ax.plot(xx, normal(xx, mu = 23, sigma = 0.5), color = "k", lw = 4, \
        label = "Poor knowledge")
ax.legend()
plt.tight_layout()
plt.savefig("./fig000.png")

aa = ax.plot(xx, normal(xx, mu = 23.5, sigma = 0.3), color = "k", ls = ":", lw = 4, \
        label = "Medium knowledge")
ax.legend()
plt.savefig("./fig001.png")

bb = ax.plot(xx, normal(xx, mu = 24, sigma = 1e-3), color = "k", ls = "--", \
        label =  "Full knowledge")
ax.legend()
plt.savefig("./fig002.png")

aa[0].remove()
bb[0].remove()

#ax.fill_between(xx[xx > -1e9], normal(xx[xx > -1e9], mu = 23, sigma = 0.5), label = "Event 1")  
ax.fill_between(xx[xx > -1e9], np.full(len(xx[xx > -1e9]), 1e9), label = "Event 1") 
ax.fill_between(xx[xx > 22.5], np.full(len(xx[xx > 22.5]), 1e9), label = "Event 2") 
ax.fill_between(xx[xx > 23.0], np.full(len(xx[xx > 23.0]), 1e9), label = "Event 3") 
#ax.fill_between(xx[xx > 22.5], normal(xx[xx > 22.5], mu = 23, sigma = 0.5), label = "Event 2")
##ax.fill_between(xx[xx > 22.5], normal(xx[xx > 22.5], mu = 23, sigma = 0.5), label = "Event 2")
#ax.fill_between(xx[xx > 23.0], normal(xx[xx > 23.0], mu = 23, sigma = 0.5), label = "Event 3")
ax.legend()
plt.savefig("./fig003.png")

aa = ax.plot(xx, normal(xx, mu = 23.5, sigma = 0.3), color = "k", ls = ":", lw = 4, \
        label = "Medium knowledge")
ax.legend()
plt.savefig("./fig004.png")

bb = ax.plot(xx, normal(xx, mu = 24, sigma = 1e-3), color = "k", ls = "--", \
        label =  "Full knowledge")
ax.legend()
plt.savefig("./fig005.png")


ax = axs[1]
ax.grid(axis = "y")
ax.set_axisbelow(True)
p1 = 1 - norm.cdf((23 - 22.5) / 0.5)
p2 = 1 - norm.cdf((23 - 23) / 0.5) - p1
p3 = 1 - (p1 + p2)
aa = ax.bar(1, p1)
bb = ax.bar(2, p2)
cc = ax.bar(3, p3)
ax.set_xticks([1, 2, 3])
ax.set_xlim(-0.5, 4.5)
ax.set_xticklabels(["Event 1\nT<22.5", "Event 2\n22.5<T<23", "Event 3\nT>23.0"])
ax.set_ylabel("Probability")
plt.savefig("./fig006.png")

aa.remove()
bb.remove()
cc.remove()
ax.set_ylabel("")
ax.set_yticklabels([""])
ax.grid(axis = "x")
ax.grid()
ax.text(1, 0.25, "NO", color = "red", va = "center", ha = "center", fontsize = 30, rotation = 30)
ax.text(2, 0.25, "NO", color = "red", va = "center", ha = "center", fontsize = 30, rotation = 30)
ax.text(3, 0.25, "YES", color = "green", va = "center", ha = "center", fontsize = 30, rotation = 30)
plt.savefig("./fig007.png")





