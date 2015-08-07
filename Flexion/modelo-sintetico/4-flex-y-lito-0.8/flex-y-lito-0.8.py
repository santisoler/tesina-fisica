import os
import numpy as np
import matplotlib.pyplot as plt
#~ from lib.double_window_selection import DoubleWindowSelection
from lib.interactive_class import InteractivePlot
from lib.flex_y_lito_functions import *
from scipy import fftpack

"""
En esta variacion del algoritmo, siempre que obtengamos una raiz por flexion
lo haremos considerando toda la grilla de topografia (carga) y el calculo de
rms lo haremos entre ventanas pequenias de las grillas de raiz invertida
y raices flexurales.

En esta version permitimos guardar los datos en un archivo.

"""


## DATOS
outfilename = "../grillas/espesor-elastico.txt"
E = 1e11
nu = 0.25
Te_min = 0
Te_max = 50000
n_Te = 101
padding = 0

# Lectura de archivos
topo = np.loadtxt("../grillas/topografia.npy")
bg = np.loadtxt("../grillas/bouguer-invertida.npy")
w = np.loadtxt("../grillas/deflexion-invertida.npy")
datos = np.loadtxt("../grillas/datos.txt")
x1, x2, y1, y2, nx, ny, rho_t, rho_c, rho_m, t = datos[:]


## Grilla de trabajo
x = np.linspace(x1,x2,nx)
y = np.linspace(y1,y2,ny)
x, y = np.meshgrid(x,y)
print "dx=", abs(x[0][1] - x[0][0])

## Archivo de salida
try:
    os.stat(outfilename)
    Te_done = np.loadtxt("../grillas/espesor-elastico.txt", skiprows=2, unpack=True)
    x_done, y_done = Te_done[0], Te_done[1]
    print "File", outfilename, "opened"    
except:
    outfile = open(outfilename, 'w')
    outfile.write("All quantity in meters\n")
    outfile.write("x_center\ty_center\twin width\tTe\tw rms\n")
    outfile.close()
    x_done, y_done = [], []
    print "File", outfilename, "created"
    
#~ outfile = open(outfilename, 'a')


## ----------
## | PLOTEO |
## ----------
## CONSTRUYO LAYOUT
import matplotlib.gridspec as gridspec
fig = plt.figure()
fig.subplots_adjust(bottom=0.05, top=0.95, left=0.04, right=0.97)
gs = gridspec.GridSpec(3, 3, width_ratios=[1.5,1.5,1])
ax1 = plt.subplot(gs[:,0])
ax2 = plt.subplot(gs[:,1])
ax3 = plt.subplot(gs[0,2])
ax4 = plt.subplot(gs[1,2])
ax5 = plt.subplot(gs[2,2])

ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax4.set_aspect('equal')
ax5.set_aspect('equal')

## Topografia (ax1)
topo_plot = ax1.contourf(x,y,topo,350)
ax1.plot(x_done, y_done, 'o', color='k', ms=2)
ax1.set_xticklabels(['%g' % (0.001 * l) for l in ax1.get_xticks()])
ax1.set_yticklabels(['%g' % (0.001 * l) for l in ax1.get_yticks()])
ax1.set_xlabel("Topografia")
cbar, ax_cbar = cbar_creator(topo, topo_plot, ax1)


## Raiz Invertida (ax2)
w_plot = ax2.contourf(x,y,w,350)
ax2.plot(x_done, y_done, 'o', color='k', ms=2)
ax2.set_xticklabels(['%g' % (0.001 * l) for l in ax2.get_xticks()])
ax2.set_yticklabels(['%g' % (0.001 * l) for l in ax2.get_yticks()])
ax2.set_xlabel("Deflexion Invertida")
cbar, ax_cbar = cbar_creator(w, w_plot, ax2)

## Interactivo
interactivo = InteractivePlot(x, y, topo, w, t, padding,
                              (rho_t,rho_c,rho_m,E,nu),
                              (Te_min,Te_max,n_Te),
                              (ax1,ax2,ax3,ax4,ax5),
                              outfilename
                              )
plt.show()
