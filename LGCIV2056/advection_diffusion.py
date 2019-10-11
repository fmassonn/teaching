import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation


# Advection - diffusion with mean flow
# ------------------------------------

# Time and space coordinates
eps = 1e-9
dt = 300 # time step, s
t = np.arange(dt, 12 * 3600.0, dt)
nt = len(t)
x = np.linspace(- 6, 6, 1000)

K = 1e-5#1e-5 # m2/s, typical value for gas

# Mean flow: x, y, z components of velocity
U= 0 * 1e-4 # m/s

# Mesh: space-time
tt, xx = np.meshgrid(t, x)

# Solution  
C = 1.0 / (2.0 * np.sqrt(K * np.pi * tt)) * np.exp(- (xx - U * tt) ** 2 / (4.0 * K * tt))



fig, ax = plt.subplots()

#text = ax.text(0,0,0)
#scatter = ax.scatter(np.random.randint(0,10,5), np.random.randint(0,20,5),marker='+')


def update(jt):
    plt.cla()
    plt.ylim(0.0, 2.0)
    plt.xlim(-5.0, 5.0)
    plt.xlabel("x [m]")
    plt.ylabel("C [-]")
    plt.grid()
    ax.plot(xx[:, 0], C[:, jt], lw = 3, color = "k", label = str(t[jt]))
    plt.title("t = " + str(t[jt]) + " s")
    #text.set_position((iternum, iternum))
    #text.set_text(str(iternum))
    #scatter.set_offsets(np.random.randint(0,10,(5,2)))

ani = animation.FuncAnimation(fig, update, frames = nt, interval=500, blit=False,
                              repeat_delay=2000)
writer = PillowWriter(fps=10)
ani.save("./fig.gif", writer = writer)