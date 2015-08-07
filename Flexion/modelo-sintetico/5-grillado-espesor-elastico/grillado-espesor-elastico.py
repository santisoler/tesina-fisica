import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

"""
Grillamos los datos de espesor elastico obtenidos con flex-y-lito.
"""


def double_contourf_plot(x, y, z1, z2, title1="1", title2="2"):
    import matplotlib.pyplot as plt
    vmin, vmax = min(np.nanmin(z1), np.nanmin(z2)), max(np.nanmax(z1), np.nanmax(z2))

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

# Lectura de archivos
datos = np.loadtxt("../grillas/datos.txt")
x1, x2, y1, y2, nx, ny, rho_t, rho_c, rho_m, t = datos[:]
Te_datos= np.loadtxt("../grillas/espesor-elastico.txt", skiprows=2, unpack=True)
#~ Te_datos= np.loadtxt("../grillas/espesor-elastico-automatico.txt", skiprows=2, unpack=True)
x_points, y_points = Te_datos[0], Te_datos[1]
Te_scatter = Te_datos[3]


## Grilla de trabajo
x = np.linspace(x1,x2,nx)
y = np.linspace(y1,y2,ny)
x, y = np.meshgrid(x,y)

## Ploteo ubicaciones de puntos
plt.axes(aspect='equal')
plt.plot(x_points, y_points, 'o')
plt.show()

## Te grillado
Te = griddata((x_points, y_points), Te_scatter, (x, y), method='cubic')
np.savetxt("../grillas/espesor-elastico-grillado.npy", Te)

# Te original
Te_original = np.loadtxt("../grillas/espesor-elastico-original.npy")


## Ploteo
double_contourf_plot(x, y, Te, Te_original, title1="Te por flexion", title2="Te original")


extremo = max(abs(np.nanmin(Te-Te_original)), abs(np.nanmax(Te-Te_original)))
plt.axes(aspect='equal')
plt.contourf(x, y, Te-Te_original, 150, cmap='coolwarm', vmin=-extremo, vmax=extremo)
plt.colorbar()
plt.title("Te - Te_original")
plt.show()
