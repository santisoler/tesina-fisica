import numpy as np
from scipy import fftpack
import matplotlib.pyplot as plt
from garcia_functions import *
  
def ploteo(x, y, z, contours=150, title=""):
    plt.axes(aspect='equal')
    plt.contourf(x, y, z, contours)
    plt.colorbar()
    plt.title(title)
    plt.show()


## LEO ARCHIVO DATOS
datos = np.loadtxt("grillas/datos.txt")
x1, x2, y1, y2, nx, ny, rho_t, rho_c, rho_m, t = datos[:]
E = 1e11
nu = 0.25


## GRILLA DE TRABAJO
x = np.linspace(x1, x2, nx)
y = np.linspace(y1, y2, ny)
x, y = np.meshgrid(x, y)


## TOPOGRAFIA
h = 3000.
topo = 0*x
topo[x==0] = h

ploteo(x, y, topo, title="Topo")

## Flexion variable
Te_0 = 6000.
Te = Te_0*np.ones(np.shape(x))
Te += (x + y)/np.max(x)*1000.
#~ Te[y>=0] += 1000.
#~ Te[y<0] -= 1000.
ploteo(x, y, Te, title="Te")

w0 = deflection_calculation(x, y, topo, rho_t, rho_c, rho_m, Te_0, E, nu)
w = garcia_iteration(x, y, w0, Te, Te_0, rho_t, rho_c, rho_m, E, nu, rms_dif_tolerance=1e-2)

ploteo(x, y, Te, title="Te")
ploteo(x, y, w0, title="w0")
ploteo(x, y, w, title="w")

np.savetxt("deflexion_garcia.npy", w)
