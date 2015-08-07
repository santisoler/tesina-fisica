import numpy as np
import matplotlib.pyplot as plt
from funciones import *

"""
Vamos a llevar a cabo el siguiente modelo:
    a) Partimos de conocer la anomalia de Bouguer de la zona de estudio.
    b) Tomaremos como punto de partida para la inversion de la raiz, la 
       raiz que se obtiene con el algoritmo de Parker y Oldenburg, haciendo
       uso del filtro pasa bajo de Hamming.
    c) Luego vamos a proceder a corregir esta raiz iterativamente:
        1) Calculamos la anomalia que genera la raiz inicial
        2) Calculamos la diferencia con la anomalia de Bouguer
        3) Calculamos la correccion a la raiz con el metodo de Parker,
           considerando solo la aproximacion a primer orden. Para evitar
           problemas de convergencia voy a aplicar un filtro pasa bajo.
"""


def double_contourf_plot(x, y, z1, z2, title1="1", title2="2"):
    import matplotlib.pyplot as plt
    vmin, vmax = min(np.min(z1), np.min(z2)), max(np.max(z1), np.max(z2))

    ax1 = plt.subplot(121)
    ax1.set_aspect('equal')
    plt.contourf(x, y, z1, 150, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.title(title1)
    
    ax2 = plt.subplot(122)
    ax2.set_aspect('equal')
    plt.contourf(x, y, z2, 150, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.title(title2)
    plt.show()



## LEO ARCHIVOS
topo = np.loadtxt("../grillas/topografia.npy")
bg = np.loadtxt("../grillas/bouguer.npy")
w_original = np.loadtxt("../grillas/deflexion-garcia.npy")
datos = np.loadtxt("../grillas/datos.txt")

x1, x2, y1, y2, nx, ny, rho_t, rho_c, rho_m, t = datos[:]


## GRILLA DE TRABAJO
x = np.linspace(x1, x2, nx)
y = np.linspace(y1, y2, ny)
x, y = np.meshgrid(x, y)


## DEFLEXION INICIAL (PARKER-OLDENBURG)
w_inicial, rms = parker_oldenburg(x, y, bg, t, rho_m-rho_c, rms_min_dif=1e-8, order=15, cut_factor=1)
np.savetxt("../grillas/deflexion-parker-oldenburg.npy", w_inicial)

plt.axes(aspect='equal')
plt.contourf(x, y, w_inicial, 150)
plt.colorbar()
plt.title("Deflexion inicial (Parker-Oldenburg)")
plt.show()


## INVERSION MEJORADA
#~ w_invertida, bg_invertida, rms = iterative_inversion(x, y, bg, w_inicial, t, rho_m-rho_c, rms_min_dif=1e-7, cut_factor=1)
w_invertida, bg_invertida, rms = iterative_inversion(x, y, bg, w_inicial, t, rho_m-rho_c, rms_min_dif=1e-8, park_old_rms_min_dif=1e-8, park_old_order=15, cut_factor=1)
np.savetxt("../grillas/deflexion-invertida.npy", w_invertida)
np.savetxt("../grillas/bouguer-invertida.npy", bg_invertida)

plt.axes(aspect='equal')
plt.contourf(x, y, w_invertida, 150)
plt.colorbar()
plt.title("Deflexion invertida")
plt.show()

double_contourf_plot(x, y, bg, bg_invertida, title1="Bouguer original", title2="Bouguer de la inversion")

