import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack

def cbar_creator(topo,topo_plot,ax):
    import numpy as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)
    ax_cbar = divider.new_vertical(size="5%",pad=0.05)
    fig1 = ax.get_figure()
    fig1.add_axes(ax_cbar)
    
    ticks_min = np.min(topo) + (1000 - np.min(topo)%1000)
    ticks_max = np.max(topo) + (1000 - np.max(topo)%1000)
    ticks = np.linspace(ticks_min, ticks_max, 5)
    cbar = plt.colorbar(topo_plot, orientation='horizontal',cax=ax_cbar, ticks=ticks, ticklocation="top")
    return cbar, ax_cbar                                        
                                    

                                        
def inside_points(x,y,scalars,(x1,x2,y1,y2)):
    dx = abs(x[0][1] - x[0][0])
    dy = abs(y[1][0] - y[0][0])
    assert dx == dy, "dx != dy"
    x0, y0 = np.min(x), np.min(y)
    imin, imax = int(abs(x1-x0)/dx), int(abs(x2-x0)/dx) + 1
    jmin, jmax = int(abs(y1-y0)/dy), int(abs(y2-y0)/dy) + 1
    return x[jmin:jmax, imin:imax], y[jmin:jmax, imin:imax], [s[jmin:jmax, imin:imax] for s in scalars]
            
    
def w_possible_list(x, y, topo, padding, t, rho_t, rho_c, rho_m, E, nu, Te_min, Te_max, n_Te):
    Te_list = np.linspace(Te_min, Te_max, n_Te)
    w_list = [deflection_calculation(x,y,topo,padding,rho_t,rho_c,rho_m,Te_i,E,nu) for Te_i in Te_list]
    return Te_list, w_list
            
            
def rms_calculation(x1,x2):
    ny, nx = np.shape(x1)
    return np.sqrt(np.sum(np.abs(x1-x2)**2)/(nx*ny))


def phi_e(k, Te, rho_c, rho_m, E, nu):
    D = Te**3*E/12./(1-nu**2)
    num = D*k**4
    return 1/(num/(rho_m-rho_c)/9.8 + 1)


def deflection_calculation(x,y,topo,padding,rho_t,rho_c,rho_m,Te,E,nu):
    # Padding
    if padding != 0:
        npad = padding*np.shape(x)[0]
        nad = int(npad)
        topo = np.pad(topo,(npad,npad),'constant',constant_values=0)
    ny, nx = np.shape(topo)
    # Freq
    dx = abs(x[0][1]-x[0][0])
    fx = fftpack.fftfreq(nx,dx)
    fy = fftpack.fftfreq(ny,dx)
    fx, fy = np.meshgrid(fx, fy)
    k = 2*np.pi*np.sqrt(fx**2+fy**2)
    F_w = rho_t/(rho_m - rho_c)*phi_e(k,Te,rho_c,rho_m,E,nu)*fftpack.fft2(topo)
    w = np.real(fftpack.ifft2(F_w))
    if padding != 0:
        w = w[npad:-npad].T[npad:-npad].T
    return w
