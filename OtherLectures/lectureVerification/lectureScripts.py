#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 20:37:34 2022

@author: massonnetf
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats
from scipy import stats
matplotlib.rcParams['font.family'] = "Arial Narrow"


np.random.seed(2)
# Number of time steps
nt = 3000

t = np.arange(nt)


alpha = 0.0
mu = 19.0
sig = 1.0

o0 = mu
ot = np.full(nt, np.nan)
ot[0] = o0



for jt in np.arange(nt - 1):
    ot[jt + 1] = alpha * (ot[jt] - mu) + sig * np.random.randn() + mu
    

fig, ax = plt.subplots(1, 1, figsize = (8, 3), dpi = 150)

ax.scatter(t, ot, 50, marker = "s", color = "k",  \
           label = "Verification data", zorder = 100)

ax.set_xticks(np.arange(2000, 2030, 1))
ax.set_xlim(2012 - 0.5, 2021 + 0.9)
ax.grid()
ax.set_ylabel("°C")
ax.set_axisbelow(True)
ax.set_title("Sea surface temperature, De Haan, 15th August (synthetic data)")
ax.set_ylim(15.0, 23.0)
#fig.legend(ncol = 4, bbox_to_anchor=(0.87, -0.1), fontsize = 10)
counterFigure = 1
fig.tight_layout()
fig.savefig("./fig" + str(counterFigure).zfill(3) + ".png"); counterFigure += 1


# Defining events

thr2 = 20.0
ax.fill_between(t, thr2 * np.ones(len(t)), 100 * np.ones(len(t)),  \
                color = [0.8, 0.0, 0.2], alpha = 0.5, label = "Warm")
    
    
#fig.legend(ncol = 4, bbox_to_anchor=(0.87, -0.1), fontsize = 10)
fig.savefig("./fig" + str(counterFigure).zfill(3) + ".png"); counterFigure += 1


thr1 = 18.0
ax.fill_between(t, thr1 * np.ones(len(t)), thr2 * np.ones(len(t)),  \
                color = [1.0, 0.5, 0.0], alpha = 0.5, label = "Normal")
#fig.legend(ncol = 4, bbox_to_anchor=(0.87, -0.1), fontsize = 10)
fig.savefig("./fig" + str(counterFigure).zfill(3) + ".png"); counterFigure += 1


ax.fill_between(t, -100 * np.ones(len(t)), thr1 * np.ones(len(t)),  \
                color = [0.0, 0.5, 1.0], alpha = 0.5, label = "Cold")
#fig.legend(ncol = 4, bbox_to_anchor=(0.87, -0.1), fontsize = 10)
fig.savefig("./fig" + str(counterFigure).zfill(3) + ".png"); counterFigure += 1


# Event = in warm
ax.scatter(t[ot > thr2], ot[ot > thr2], 100, marker = "v", \
           color = "green", zorder = 101, label = "Event occurs")
ax.scatter(t[ot < thr2], ot[ot < thr2], 100, marker = "x", \
           color = "red", zorder = 101, label = "Event does not occur")


#fig.legend(ncol = 4, bbox_to_anchor=(0.87, -0.1), fontsize = 10)
fig.savefig("./fig" + str(counterFigure).zfill(3) + ".png"); counterFigure += 1




# Plot the forecast PDFs
xpdf = np.linspace(10, 30, 1000)
leg = "Forecast PDF"

# Our four forecasters

name = ["Alice", "Bob", "Charles", "Damien", "Edwin", "50%-50%"]
nF   = len(name) # Nb. forecasters

muF  = [ot + 0.5 * np.random.randn(nt), \
        mu * np.ones(nt), \
        ot + 0.5 * np.random.randn(nt), \
        ot + 2 * np.random.randn(nt), \
        ot + 0.1 * np.random.randn(nt),
        thr2 * np.ones(nt)]
stF  = [sig, \
        sig, \
        2 * sig, \
        0.5 * sig, \
        0.1 * sig, \
        100 * sig]

# Go through all forecasts
# ------------------------

for jF in np.arange(nF):
        
    myPlots = list()
    myTexts = list()

    for jt in t:
        xpdf = np.linspace(10.0, 30.0, 1000)
        pdf = scipy.stats.norm.pdf(xpdf, muF[jF][jt], stF[jF])
        pdf /= (1.2 * np.max(pdf))
        
        
        # Plot PDF
        if jt >= 2012 and jt <= 2021:
            a = ax.plot(jt + pdf, xpdf, np.ones(len(xpdf)), color = "b")
            
            myPlots.append(a)
  
            # relevant probabilities
            prob = 1 - scipy.stats.norm.cdf(thr2, muF[jF][jt], stF[jF])
            b = ax.text(jt + 0.3, 15.5, str(int(prob * 100)) + " %", \
                        color = "Blue")
            myTexts.append(b)


    fig.savefig("./series" + name[jF] + ".png")

    # delete for following graphs
    [m.pop(0).remove() for m in myPlots]
    [m.remove() for m in myTexts]

    
    # Derive event probabilities
    prob = np.array([1 - scipy.stats.norm.cdf(thr2, muF[jF][jt], stF[jF])\
                      for jt in np.arange(nt)])

    # Attributes        
    obs = list()
    n = list()
    
    REL = list()
    RES = list()
    print("Forecaster " + str(name[jF]))
    print("----------")
    for myBin in np.arange(0.0, 1.0, 0.1):
    
        # Frequency relative to all occurences
        tmp = np.sum(1.0 * (ot >= thr2)[(prob >= myBin) * \
                                        (prob < myBin + 0.1)]) / nt
        # Frequency related to the number occurences in that class
        rel = np.sum(1.0 * (ot >= thr2)[(prob >= myBin) * \
                                        (prob < myBin + 0.1)]) / \
           np.sum(1.0 * (prob >= myBin) * (prob < myBin + 0.1))
           
        obs.append(rel)
    
        ni = np.sum(1.0 * (prob >= myBin) * (prob < myBin + 0.1))
        n.append(ni)
    
    
        tmp2 = np.sum(1.0 * (ot < thr2)[(prob >= myBin) * \
                                        (prob < myBin + 0.1)]) / nt
            

        #print("Probability class: " + str(np.round(myBin, 2)) + \
        #      '--' + str(np.round(myBin + 0.1, 2)))
        
        #print(" p(o = 1|f) = "  + str(np.round(tmp  * 100)) + " %")
        #print(" p(o = 0|f) = "  + str(np.round(tmp2 * 100)) + " %")
        #print(" freq = " + str(np.round(np.sum(1.0 * (prob >= myBin) * \
        #                              (prob < myBin + 0.1)) / nt, 2)))
            
        #print(" p(f): " + str(100 * np.round(np.sum((prob >= myBin) *\
        #                                           (prob < myBin + 0.1)) \
        #                                     / nt, 2)) + " %")
            

    BS = np.mean((prob - (ot > thr2)) ** 2)
    
   

    print("BS : " + str(np.round(BS, 2)))

    # Graphical representation
    # ------------------------
    figR, axR = plt.subplots(figsize = (3, 3), dpi = 300)
    axR.scatter(np.arange(0.05, 1.0, 0.1), obs, n, marker = ".", \
                color = "blue", alpha = 0.8)

    axR.plot((0, 1.0), (0.0, 1.0), color = [0.3, 0.3, 0.3], linestyle = "--")
    axR.set_xlabel("Event forecast probabilities")
    axR.set_ylabel("Observed event frequencies")
    axR.grid()
    axR.set_axisbelow(True)
    figR.tight_layout()
    figR.savefig("./reliability" + str(name[jF]) + ".png", dpi = 300)
    plt.close(figR)
    

    # Graph of discrimination
    # -----------------------
    xpdf = np.linspace(0, 1, 1000)
    figD, axD = plt.subplots(1, 1, figsize = (3, 2), dpi = 300)
    if name[jF] == "Bob": # Clim forecast --> kernel does not work
        axD.plot((prob[0], prob[0]), (0, 100), color = [0.0, 1.0, 0.0], lw = 4 )
        axD.plot((prob[0], prob[0]), (0, 100), color = "Red", lw = 2 )
        
    else:
        # event happened
        kernel = stats.gaussian_kde(prob[ot > thr2])
        pdf = kernel(xpdf).T
        scipy.stats.norm.pdf(xpdf, pdf )
    
        
        axD.plot(xpdf, pdf, lw = 4, color = [0.0, 1.0, 0.0])
    
        # event did not happen
        kernel = stats.gaussian_kde(prob[ot < thr2])
        pdf = kernel(xpdf).T
        scipy.stats.norm.pdf(xpdf, pdf )
        axD.plot(xpdf, pdf, lw = 4, color = "Red")
        axD.set_xlim(0.0, 1.0)
        
    axD.set_xlabel("Forecast probability")
    figD.tight_layout()
    figD.savefig("discrim" +name[jF] + ".png", dpi = 300)
    plt.close(figD)


    # Graph of sharpness
    # ------------------
    figS, axS = plt.subplots(1, 1, figsize = (3, 2), dpi = 300)
    theseBins = np.arange(0.0, 1.1, 0.1)
    H, bins = np.histogram(prob, bins = theseBins)
    axS.bar(bins[1:] + 0.05, H / np.sum(H) * 100, width = 0.09)
    axS.set_title(name[jF])
    axS.set_xlabel("Forecast probabilities")
    axS.set_ylabel("Relative frequencies (%)")
    figS.tight_layout()
    figS.savefig("./sharp" + name[jF] + ".png", dpi = 300)
    

# Other graphs
fig, ax = plt.subplots(1, 1, figsize = (3, 3), dpi = 300)
x = np.linspace(0, 1, 1000)
ax.plot(x, x * (1 - x))
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)
ax.set_ylabel("Brier score")
fig.tight_layout()
fig.savefig("./brier01.png", dpi = 300)

ax.plot(x, 2 * x * (1 - x), color = "orange")
fig.tight_layout()
fig.savefig("./brier02.png", dpi = 300)

ax.plot(x,0.25 * np.ones(len(x)), color = "green")
fig.tight_layout()
fig.savefig("./brier03.png", dpi = 300)

ax.plot(x,0.33 * np.ones(len(x)), color = "purple")
fig.tight_layout()
fig.savefig("./brier04.png", dpi = 300)


# CRPS
fig, ax = plt.subplots(1, 1, dpi = 300, figsize = (4 ,3))
xpdf = np.linspace(10, 30, 1000)
jt = 2019
# Alice
jF = 0
pdf = scipy.stats.norm.pdf(xpdf, muF[0][jt], stF[jF])
cdf = scipy.stats.norm.cdf(xpdf, muF[0][jt], stF[jF])
a = ax.plot(xpdf, pdf, color = "b")

ax.plot((ot[jt], ot[jt]), (0, 1), "k--")
ax.set_xlim(15.0, 25.0)
ax.set_ylim(0.0, 1.0)
ax.set_xlabel("SST ($^\circ$ C)")
fig.tight_layout()
fig.savefig("./crps01.png", dpi = 300)
a.pop(0).remove()

xtmp = xpdf[xpdf > ot[jt]]
pdftmp = pdf[xpdf > ot[jt]]
aaa = ax.fill_between(xtmp, pdftmp)
fig.savefig("./talagrand01.png", dpi = 300)
aaa.remove()

ax.plot(xpdf, cdf, color = "b", linestyle = "--")
fig.savefig("./crps02.png", dpi = 300)

xtmp = xpdf[xpdf < ot[jt]]
cdftmp = cdf[xpdf < ot[jt]]

ax.fill_between(xtmp, cdftmp, color = "blue", alpha = 0.5)
xtmp = xpdf[xpdf > ot[jt]]
cdftmp = cdf[xpdf > ot[jt]]

ax.fill_between(xtmp, cdftmp, 1.0, color = "blue", alpha = 0.5)
fig.savefig("./crps03.png", dpi = 300)



# Talagrands
for jF in np.arange(nF):
    quantiles = [1 - scipy.stats.norm.cdf(ot[jt], muF[jF][jt], stF[jF]) for jt in np.arange(nt)]
    fig, ax = plt.subplots(1, 1, dpi = 300, figsize = (4 ,3))
    theseBins = np.arange(0.0, 1.1, 0.1)
    H, bins = np.histogram(quantiles, bins = theseBins)
    ax.bar(bins[1:] + 0.05, H / np.sum(H) * 100, width = 0.09)
    ax.set_ylabel("Relative frequency (%)")
    ax.set_xlabel("Quantile defined by observation")
    ax.set_title(name[jF])
    fig.tight_layout()
    fig.savefig("./talagrand" + name[jF] + ".png")
    