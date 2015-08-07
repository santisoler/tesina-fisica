import numpy as np
from scipy import fftpack
import matplotlib.pyplot as plt
from fatiando import utils
#~ from fatiando.constants import G
#~ from fatiando.vis import mpl
from fatiando.mesher import Prism
from fatiando.gravmag import prism

"""
GENERO LA ANOMALIA DE BOUGUER PRODUCIDA POR LA RAIZ QUE OBTUVIMOS DE APLICAR
EL METODO DE GARCIA

"""

def anomaly_calculation(x, y, w, t, rho):
    from fatiando import utils
    from fatiando.mesher import Prism
    from fatiando.gravmag import prism

    ny, nx = np.shape(x)
    prismas = []
    for i in xrange(ny-1):
        for j in xrange(nx-1):
            prisma = Prism(y[i][j], y[i+1][j+1], x[i][j], x[i+1][j+1],
                            t, t + (w[i][j] + w[i+1][j+1])/2.)
            prismas.append(prisma)
    
    gz = -prism.gz(y.ravel(),x.ravel(),np.zeros(nx*ny),prismas,dens=rho)
    gz = utils.mgal2si(np.reshape(gz,(ny,nx)))
    return gz
    

## LEO ARCHIVOS
topo = np.loadtxt("../grillas/topografia.npy")
w_garcia = np.loadtxt("../grillas/deflexion-garcia.npy")
datos = np.loadtxt("../grillas/datos.txt")
x1, x2, y1, y2, nx, ny, rho_t, rho_c, rho_m, t = datos[:]
E = 1e11
nu = 0.25


## GRILLA DE TRABAJO
x = np.linspace(x1, x2, nx)
y = np.linspace(y1, y2, ny)
x, y = np.meshgrid(x, y)


## Grilla de trabajo
x = np.linspace(x1,x2,nx)
y = np.linspace(y1,y2,ny)
x, y = np.meshgrid(x,y)
#~ shape = np.shape(x)
#~ x_dat, y_dat = x.ravel(), y.ravel()

## ----------------------------------------------------------------------------


## CALCULO RAIZ
bg = anomaly_calculation(x, y, w_garcia, t, rho_m-rho_c)


## GUARDO ARCHIVOS
np.savetxt("../grillas/bouguer.npy", bg)


## CONTOURFS
plt.axes(aspect='equal')
plt.contourf(x ,y , topo, 150)
plt.colorbar()
plt.title("Topo")
plt.show()

plt.axes(aspect='equal')
plt.contourf(x, y, w_garcia, 150)
plt.colorbar()
plt.title("Deflexion")
plt.show()

plt.axes(aspect='equal')
plt.contourf(x, y, bg, 150)
plt.colorbar()
plt.title("Anomalia")
plt.show()


