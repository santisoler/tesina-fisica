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
x1, x2, y1, y2, nx, ny, rho_t, rho_c, rho_m, t = datos[:]
E = 1e11
nu = 0.25


## GRILLA DE TRABAJO
x = np.linspace(x1, x2, nx)
y = np.linspace(y1, y2, ny)
x, y = np.meshgrid(x, y)



ploteo(x, y, topo, title="Topo")



## Flexion variable
#~ Te_0 = 6000.
#~ Te = Te_0*np.ones(np.shape(x))
#~ Te[y>=0] += 1000.
#~ Te[y<0] -= 1000.

Te = 6000 - (y - y.max())*2e-2
np.savetxt("../grillas/espesor-elastico-original.npy", Te)
ploteo(x, y, Te, title="Te")



raw_input("Calcular deflexion?")
Te_0 = 10000
w0 = deflection_calculation(x, y, topo, rho_t, rho_c, rho_m, Te_0, E, nu)
w = garcia_iteration(x, y, w0, Te, Te_0, rho_t, rho_c, rho_m, E, nu, rms_dif_tolerance=1e-13)

#~ np.savetxt("../grillas/topografia.npy", topo)
np.savetxt("../grillas/deflexion-garcia-1.npy", w)


ploteo(x, y, Te, title="Te")
ploteo(x, y, w0, title="w0")
ploteo(x, y, w, title="w")

