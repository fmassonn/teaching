import numpy as np
import matplotlib.pyplot as plt

# Planck's constant
h = 6.62607015e-34  # J*s
# Boltzmann's constant
k = 1.380649e-23  # J/K
# Speed of light
c = 299792458  # m/s

# Temperature of the first black body in Kelvin (Sun)
T1 = 5760
# Temperature of the second black body in Kelvin (Earth)
T2 = 288

# Wavelength range in meters
wavelengths = np.logspace(-7, -4, 1000)

# Irradiance spectra for the two black bodies
irradiance1 = (2*h*c**2/wavelengths**5) * (1/(np.exp(h*c/(wavelengths*k*T1))-1))
irradiance2 = (2*h*c**2/wavelengths**5) * (1/(np.exp(h*c/(wavelengths*k*T2))-1))

#irradiance1 /=np.max(irradiance1)
#irradiance2 /=np.max(irradiance2)

# Plot the irradiance spectra on the same figure
fig, ax = plt.subplots(2, 1)
ax[0].plot(wavelengths*1e6, irradiance1, label='T = {} K'.format(T1), color = "orange")

ax[1].plot(wavelengths*1e6, irradiance2, label='T = {} K'.format(T2), color = "red")
for a in ax.flatten():
    a.set_xscale('log')
    a.set_xlabel("Longueur d'onde (micromètres)")
    a.set_ylabel('Irradiance (W/m$^2$/nm/sr)')
    a.legend()

fig.tight_layout()

plt.savefig("./fig.png", dpi = 300)
