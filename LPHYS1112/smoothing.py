# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 08:15:54 2019

@author: massonnetf
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

def gsmooth(t_in, x_in, tau, t_out = None, method = None, check = False):
    """
    t_in      is a sorted N by 1 numpy array of independent coordinates
    x_in      is a N by 1 numpy array of data
    tau       is the decorrelation time scale for the smoothing window
              (in units of the time coordinate). When data is tau away from
              the center of the window, it is given a weight 1/e
    t_out  is a M by 1 numpy array of independent coordinates for smoothing. If
           no value is specified, a default number of 1000 points is taken
    method    is the method used to propagate data outside the window
    check     is a boolean to plot smoothing outside window
    
    returns
    
    gaussian smoothed version of x at t_out
    """
    
    # Convert to np array if necessary
    if type(t_in) is list:
        t_in = np.array(t_in)
    if type(x_in) is list:
        x_in = np.array(x_in)
    
    is_sorted = np.all(t_in[:-1] < t_in[1:])
    
    if not is_sorted:
        sys.exit("(gsmooth) t is not sorted")
        
    # Remove all NanS
    t_in = t_in[~np.isnan(x_in)]
    x_in = x_in[~np.isnan(x_in)]
   
    # 1. Mirror the data to sort the edges
    if method == "csym":  # Central symmetry

        x =  np.concatenate((np.array(2.0 * x_in[0] - x_in )[-1:0:-1], \
                             x_in                                    , \
                             np.array(2.0 * x_in[-1] - x_in)[-2::-1]     
                           ))
        t =  np.concatenate((np.array(2.0 * t_in[0] - t_in)[-1:0:-1],  \
                             t_in                                   ,
                             np.array(2.0 * t_in[-1] - t_in)[-2::-1]
                           ))
    elif method == None:
        pass
    else:
        sys.exit("(gsmooth) method unknown")
        
    if t_out == None:
        t_out = np.linspace(t[0], t[-1], 3* 1000) # 3 windows of 1000
    
    delta_t  = t_out[:, None] - t
    #weights  = np.exp(- np.abs(delta_t / tau)) # option exponential
    weights  = np.exp(- (delta_t / tau) ** 2)    # option gaussian
    weights /= np.sum(weights, axis = 1, keepdims=True)
    x_out    = np.dot(weights, x)

    start = np.where(t_out > t_in[0])[0][0]
    end   = np.where(t_out > t_in[-1])[0][0]

    if check:
        plt.figure("check")    
        plt.scatter(t_in,  x_in,  50, marker = "*", color = "black")
        plt.scatter(t, x, 10, marker = "*", color = "green")
        plt.plot(t_out, x_out, "red")
        plt.grid()
        plt.savefig("check.png", dpi = 300)
        plt.close("check")
    
    return t_out[start:end], x_out[start:end]
        
def nurbs(t_in, x_in, order = 3, check = False):
    """
    t_in      is a sorted N by 1 numpy array of independent coordinates
    x_in      is a N by 1 numpy array of data
    order     is the order of the fit
    method    is the method used to propagate data outside the window
    check     is a boolean to plot smoothing outside window
    
    returns
    
    NURBS smoothed version of x_in at t_out
    """
    
    # Convert to np array if necessary
    if type(t_in) is list:
        t_in = np.array(t_in)
    if type(x_in) is list:
        x_in = np.array(x_in)
    
    is_sorted = np.all(t_in[:-1] < t_in[1:])
    
    if not is_sorted:
        sys.exit("(nurbs) t is not sorted")
        
    # Remove all NanS
    t_in = t_in[~np.isnan(x_in)]
    x_in = x_in[~np.isnan(x_in)]
   
    # Number of control points
    n = len(t_in)
    
    x1   = t_in
    x2   = x_in
    
    colors = [np.random.random(3) for jn in range(n)]
    
    if check:    
        plt.figure("fig1", figsize = (3, 3))
        plt.scatter(x1, x2, 100, marker = "*", color = colors)
    
    
    # Degree
    p = order
    
    # Basis functions
    # The basis functions are by layers: the first layer are piece-wise constant
    # functions, the second layer are linear functions (one less), etc.
    # The number of layers must be such that there are as many basis functions
    # of order "p" as control points "n"
    #
    # There must be "n"   basis functions of order "p"
    # So            "n+1" ------------------------ "p-1"
    # So            "n+p" ------------------------ "0" (piecewise cst functions)
    #
    # Since each piecewise function is defined over an interval and a set of 
    # "m" intervals require "m + 1" knots, the basis functions must be defined
    # over a set of n+p+1 knots; this set is known as the knot vector
    
    m = n + p + 1 # Number of control knots
    
    uk = np.linspace(0, 1, m)
    #uk = np.sort(np.random.randn((m)))
    
    # Definition of basis functions
    # -----------------------------
        
    # Create support. Notice that uk = set of knots, u = continuous support
    u = np.linspace(uk[0], uk[-1], 100)
    
    # Create degree zero functions
    # N will be a list of list of basis functions
    N = list()
    mylist = [(u >= uk[jm]) * (u < uk[jm + 1]) * 1.0 for jm in range(m -1)]
    N.append(mylist)
    
    # Create higher order functions, by recurrence (Cox-de Boor formula)
    for jp in range(1, p + 1):
        mylist = [(u - uk[jm]) / (uk[jm + jp] - uk[jm]) * N[jp - 1][jm] + \
        (uk[jm + jp + 1] - u) / (uk[jm + jp + 1] - uk[jm + 1]) * \
        N[jp - 1][jm + 1] for jm in range(m - jp - 1)]
    
        N.append(mylist)
    
    if check:
        # Plot basis functions
        plt.figure("basis", figsize = (2* (p + 1), 4 * (p + 1) ))
        for jp in range(p + 1):
            if jp == p:
                cols = colors
            else:
                cols = [[0.3, 0.3, 0.3] for jm in range(m - jp - 1)]
            plt.subplot(p + 1, 1, jp + 1)
            for jm in range(m - jp - 1):
                plt.plot(u, N[jp][jm], color = cols[jm], linewidth = 4, \
                         label = "N(" + str(jm) + "," + str(jp) + ")")
            
            plt.grid()
            plt.legend()
            plt.title("Order " + str(jp))
        
            # Plot all knots
            for jm in range(m):
                plt.scatter(uk[jm], 0, 100, color = "black", zorder = 1000)
                
        plt.tight_layout()
        plt.savefig("./basis.png", dpi = 300)
        plt.close("basis")
    
    # Fitting
    # -------
    
    # A linear combination of all n basis functions of order p is taken,
    # but the contribution of each function is weighted by the corresponding
    # control point coordinate
    
    Cx1 = np.sum([[N[p][jn] * x1[jn]] for jn in range(n)], axis = 0) / np.sum(N[p], axis = 0)
    #print(np.sum([[N[p][jn] * x1[jn]] for jn in range(n)], axis = 0))
    #print(np.sum(N[p], axis = 0))
    Cx2 = np.sum([[N[p][jn] * x2[jn]] for jn in range(n)], axis = 0) / np.sum(N[p], axis = 0)
    
    if check:
        plt.figure("fig1")
        plt.plot(Cx1[0], Cx2[0], color = "blue", linewidth = 2)
        plt.grid()
        plt.gca().set_axisbelow(True)
        plt.tight_layout()
        plt.savefig("./fig.png", dpi = 300)
        plt.close("fig1")
        
    # Not first and last items (nans)
    i1 = np.where(u >= uk[1])[0][0]
    i2 = np.where(u >= uk[-2])[0][0]
    
    print(i1)
    return Cx1[0][i1:i2], Cx2[0][i1:i2] 


def test():
    #np.random.seed(2)
    nt = 6
    plt.figure("test", figsize = (6, 3))
    t = np.sort(np.random.randn((nt)))
    x = t + 0*np.random.randn((nt))
    plt.scatter(t, x, 50, marker = "*")
    t_out, x_out = gsmooth(t, x, 3, method = "csym")
    t_out, x_out = nurbs(t, x)
    plt.plot(t_out, x_out, color = "red")
    #plt.ylim(-3, 3)
    plt.savefig("./test.png", dpi = 300)
    #plt.close("test")


    

