import numpy as np
from scipy import fftpack

   
def deflection_calculation(x, y, topo, rho_t, rho_c, rho_m, Te, E, nu, padding=0):
    """
    Calculates the deflection due to a topographic load for a plate of constant
    thickness Te.
    Uses the equation:
    F[w] = rho_t/(rho_m-rho_c)  phi_e(k)  F[topo]
    """
    ny, nx = np.shape(topo)
    dx = abs(x[0][1] - x[0][0])
    dy = abs(y[1][0] - y[0][0])
    if padding != 0:
        ny_pad, nx_pad = ny*padding, nx*padding
        topo = np.pad(topo, (ny_pad,nx_pad), 'constant', constant_values=0)
    else:
        nx_pad, ny_pad = 0, 0 
        
    fx = fftpack.fftfreq(nx + 2*nx_pad, dx)
    fy = fftpack.fftfreq(ny + 2*ny_pad, dy)
    fx, fy = np.meshgrid(fx, fy)
    k = 2*np.pi*np.sqrt(fx**2 + fy**2)
    
    F_w = rho_t/(rho_m - rho_c)*phi_e(k, Te, rho_c, rho_m, E, nu)*fftpack.fft2(topo)
    w = np.real(fftpack.ifft2(F_w))
    
    if padding != 0:
        w = w[ny_pad:-ny_pad, nx_pad:-nx_pad]
    return w


def phi_e(k, Te, rho_c, rho_m, E, nu):
    D = Te**3*E/12./(1.-nu**2)
    num = (D**(1/4.)*k)**4
    return 1/(num/(rho_m-rho_c)/9.8 + 1)



def garcia_iteration(x, y, w0, Te, Te_0, rho_t, rho_c, rho_m, E, nu, rms_dif_tolerance=1e-3):
    ny, nx = np.shape(x)
    dx = abs(x[0][1] - x[0][0])
    dy = abs(y[1][0] - y[0][0])
    assert dx == dy, "dx != dy"
    fx = fftpack.fftfreq(nx, dx)
    fy = fftpack.fftfreq(ny, dy)
    fx, fy = np.meshgrid(fx, fy)
    k = 2*np.pi*np.sqrt(fx**2 + fy**2)
    
    D = E*Te**3/12./(1-nu**2)
    D0 = E*Te_0**3/12./(1-nu**2)
    Dp = D - D0
    phi = phi_e(k, Te_0, rho_c, rho_m, E, nu)*float(9.8*(rho_m-rho_c))**(-1)
    F_w0 = fftpack.fft2(w0)
    Dp_x2 = partial_x2(Dp, x, y)
    Dp_y2 = partial_y2(Dp, x, y)
    Dp_xy = partial_xy(Dp, x, y)
    
    w = w0
    rms = np.inf
    i = 0
    while True:
        w_old = w
        rms_old = rms
        T1 = laplacian(Dp*laplacian(w_old, x, y), x, y)
        T21 = Dp_x2*partial_y2(w_old, x, y)
        T22 = 2*Dp_xy*partial_xy(w_old, x, y)
        T23 = Dp_y2*partial_x2(w_old, x, y)
        fun = T1 - (1-nu)*(T21 - T22 + T23)
        F_w = F_w0 - phi*fftpack.fft2(fun)
        w = np.real(fftpack.ifft2(F_w))
        rms = rms_calculation(w, w_old)
        i += 1
        print "Iteration", i
        print "RMS:", rms, "\n"
        
        if 0 <= rms_old - rms < rms_dif_tolerance:
            print "Iterative Deflection Calculation Finished"
            return w
        elif rms > rms_old:
            print "Increasing RMS value found"
            return w_old
            
        
    return w
    
    
def rms_calculation(x1, x2):
    ny, nx = np.shape(x1)
    rms = np.sqrt(np.sum((x1 - x2)**2)/(nx*ny))
    return rms

#~ def hamming_filter_full(x, y, f, k_cut):
    #~ ny, nx = np.shape(x)
    #~ dx = abs(x[0][1] - x[0][0])
    #~ dy = abs(y[1][0] - y[0][0])
    #~ assert dx == dy, "dx != dy"
    #~ fx = fftpack.fftfreq(nx, dx)
    #~ fy = fftpack.fftfreq(ny, dy)
    #~ fx, fy = np.meshgrid(fx, fy)
    #~ k = 2*np.pi*np.sqrt(fx**2 + fy**2)
    #~ F = fftpack.fft2(f)
    #~ F = hamming_filter(F, k, k_cut)
    #~ return fftpack.ifft2(F)
    #~ 
    
    
## SPECTRAL DERIVATIVES

def hamming_filter(F, k, k_cut):
    mask = (k < k_cut)*(0.5*(1 + np.cos(np.pi*k/k_cut)))
    return F*mask

def laplacian(f, x, y):
    return partial_x2(f, x, y) + partial_y2(f, x, y)    
    
def partial_x2(f, x, y):
    ny, nx = np.shape(x)
    dx = abs(x[0][1] - x[0][0])
    dy = abs(y[1][0] - y[0][0])
    assert dx == dy, "dx != dy"
    
    fx, fy = fftpack.fftfreq(nx, dx), fftpack.fftfreq(ny, dy)
    fx, fy = np.meshgrid(fx, fy)
    kx, ky = 2*np.pi*fx, 2*np.pi*fy
    k = np.sqrt(kx**2 + ky**2)
    
    # Filtrado de funcion
    k_cut = np.pi/dx
    F_f = hamming_filter(fftpack.fft2(f), k, k_cut)
    
    # Resolucion de la derivada
    F_fx2 = -kx**2*F_f
    f_x2 = fftpack.ifft2(F_fx2)
    return np.real(f_x2)
    
def partial_y2(f, x, y):
    ny, nx = np.shape(x)
    dx = abs(x[0][1] - x[0][0])
    dy = abs(y[1][0] - y[0][0])
    assert dx == dy, "dx != dy"
    
    fx, fy = fftpack.fftfreq(nx, dx), fftpack.fftfreq(ny, dy)
    fx, fy = np.meshgrid(fx, fy)
    kx, ky = 2*np.pi*fx, 2*np.pi*fy
    k = np.sqrt(kx**2 + ky**2)
    
    # Filtrado de funcion
    k_cut = np.pi/dx
    F_f = hamming_filter(fftpack.fft2(f), k, k_cut)
    
    # Resolucion de la derivada
    F_fy2 = -ky**2*F_f
    f_y2 = fftpack.ifft2(F_fy2)
    return np.real(f_y2)
    
def partial_xy(f, x, y):
    ny, nx = np.shape(x)
    dx = abs(x[0][1] - x[0][0])
    dy = abs(y[1][0] - y[0][0])
    assert dx == dy, "dx != dy"
    
    fx, fy = fftpack.fftfreq(nx, dx), fftpack.fftfreq(ny, dy)
    fx, fy = np.meshgrid(fx, fy)
    kx, ky = 2*np.pi*fx, 2*np.pi*fy
    k = np.sqrt(kx**2 + ky**2)
    
    # Filtrado de funcion
    k_cut = np.pi/dx
    F_f = hamming_filter(fftpack.fft2(f), k, k_cut)
    
    # Resolucion de la derivada
    F_fxy = -kx*ky*F_f
    f_xy = fftpack.ifft2(F_fxy)
    return np.real(f_xy)
