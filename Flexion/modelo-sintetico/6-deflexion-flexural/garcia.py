import numpy as np
from scipy import fftpack
import matplotlib.pyplot as plt
from garcia_functions import *
  
"""
A partir de una topografia y una grilla de rigidez D obtenemos la deflexion
que genera dicha carga topografica a traves del metodo de Garcia (2014).
"""  
  
  
def ploteo(x, y, z, contours=150, title=""):
    plt.axes(aspect='equal')
    plt.contourf(x, y, z, contours)
    plt.colorbar()
    plt.title(title)
    plt.show()


## LEO ARCHIVO DATOS
datos = np.loadtxt("../grillas/datos.txt")
topo = np.loadtxt("../grillas/topografia.npy")
Te = np.loadtxt("../grillas/espesor-elastico-grillado.npy")
x1, x2, y1, y2, nx, ny, rho_t, rho_c, rho_m, t = datos[:]
E = 1e11
nu = 0.25


## GRILLA DE TRABAJO
x = np.linspace(x1, x2, nx)
y = np.linspace(y1, y2, ny)
x, y = np.meshgrid(x, y)









Te[np.isnan(Te)] = 0

print Te
Te_0 = 12000
w0 = deflection_calculation(x, y, topo, rho_t, rho_c, rho_m, Te_0, E, nu)
w = garcia_iteration(x, y, w0, Te, Te_0, rho_t, rho_c, rho_m, E, nu, rms_dif_tolerance=1e-13)

#~ np.savetxt("../grillas/topografia.npy", topo)
np.savetxt("../grillas/deflexion-flexural.npy", w)


ploteo(x, y, Te, title="Te")
ploteo(x, y, w0, title="w0")
ploteo(x, y, w, title="w")

