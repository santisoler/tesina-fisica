import numpy as np
from scipy import fftpack
from scipy.constants import G

### DE USO COMPARTIDO

def hamming_filter(F, k, k_cut):
    mask = (k < k_cut)*(0.5*(1 + np.cos(np.pi*k/k_cut)))
    return F*mask


def rms_window_cut(x1, x2, ratio=0.25):
    ny, nx = np.shape(x1)
    ny_cut, nx_cut = int(ny*ratio), int(nx*ratio)
    x1_cut = x1[ny_cut:-ny_cut, nx_cut:-nx_cut]
    x2_cut = x2[ny_cut:-ny_cut, nx_cut:-nx_cut]
    return rms_calculation(x1_cut, x2_cut)


def rms_calculation(x1, x2):
    ny, nx = np.shape(x1)
    rms = np.sqrt(np.sum((x1 - x2)**2)/(nx*ny))
    return rms



## PRINCIPALES MODULOS

def parker_oldenburg(x, y, bg, t, rho, rms_min_dif=1e3, order=9, cut_factor=1):
    
    def sum_n(k, w, order):
        import math
        sum_n = 0
        for n in range(2,order+1):
            sum_n += (-1)**n*k**(n-1)/math.factorial(n)*fftpack.fft2(w**n)
        return sum_n
    
    def root_calculation(k, w, constant_term, order, k_cut):
        F_w = constant_term + sum_n(k, w, order)
        F_w = hamming_filter(F_w, k, k_cut)
        return np.real(fftpack.ifft2(F_w))
    
    ny, nx = np.shape(x)
    dx = abs(x[0][1] - x[0][0])
    dy = abs(y[1][0] - y[0][0])
    assert dx == dy, "dx != dy"
    fx = fftpack.fftfreq(nx, dx)
    fy = fftpack.fftfreq(ny, dy)
    fx, fy = np.meshgrid(fx, fy)
    k = 2*np.pi*np.sqrt(fx**2 + fy**2)
    k_cut = cut_factor*(np.pi/t)
    
    constant_term = -1/(2*np.pi*G*rho)*np.exp(k*t)*fftpack.fft2(bg)
    
    print "Inversion Parker-Oldenburg"
    w = np.zeros((ny,nx))
    rms = np.inf
    i = 0
    while True:
        i += 1
        print "Iteracion", i
        w_old, rms_old = w, rms
        w = root_calculation(k, w, constant_term, order, k_cut)
        rms = rms_calculation(w_old, w)
        print "RMS:", rms, "\n"
        if 0 <= rms_old - rms < rms_min_dif:
            return w, rms
        elif rms > rms_old:
            print "Increasing RMS value found"
            return w_old, rms_old


def iterative_inversion(x, y, bg, w0, t, rho, rms_min_dif=1e-8, park_old_rms_min_dif=1e-3, park_old_order=15, cut_factor=1):
    w = w0
    print "INITIAL CALCULATION"
    gz = anomaly_calculation(x, y, w, t, rho)
    rms = rms_calculation(bg, gz)
    #~ rms = rms_window_cut(bg, gz)
    print "Initial RMS: ", rms, "\n"
    
    i = 0
    while True:
        i += 1
        print "ITERATION ",i
        w_old, gz_old, rms_old = w, gz, rms
        w = parker_oldenburg_sheet(x, y, bg-gz, w, t, rho, cut_factor=cut_factor, rms_min_dif=park_old_rms_min_dif, order=park_old_order)
        gz = anomaly_calculation(x, y, w, t, rho)
        rms = rms_calculation(bg, gz)
        #~ rms = rms_window_cut(bg, gz)
        print "RMS: ", rms, "\n"
        if 0 < rms_old - rms < rms_min_dif:
            print "Inversion completed"
            return w, gz, rms
        elif rms > rms_old:
            print "Increasing RMS value found"
            return w_old, gz_old, rms_old


def parker_oldenburg_sheet(x, y, delta_g, w0, t, rho, cut_factor=1, order=15, rms_min_dif=1e-3):
        
    def sum_n(k, w, w0, order):
        import math
        sum_n = 0
        for n in range(2,order+1):
            sum_n += (-1)**n*k**(n-1)/math.factorial(n)*fftpack.fft2(w**n - w0**n)
        return sum_n
    
    def root_calculation(k, w, w0, constant_term, order, k_cut):
        F_w_w0 = constant_term - sum_n(k, w, w0, order)
        F_w_w0 = hamming_filter(F_w_w0, k, k_cut)
        w_w0 = np.real(fftpack.ifft2(F_w_w0))
        return w_w0 + w0
        
    ny, nx = np.shape(x)
    dx = abs(x[0][1] - x[0][0])
    dy = abs(y[1][0] - y[0][0])
    assert dx == dy, "dx != dy"
    fx = fftpack.fftfreq(nx, dx)
    fy = fftpack.fftfreq(ny, dy)
    fx, fy = np.meshgrid(fx, fy)
    k = 2*np.pi*np.sqrt(fx**2 + fy**2)
    k_cut = cut_factor*(np.pi/t)
    
    constant_term = -1/(2*np.pi*G*rho)*np.exp(k*t)*fftpack.fft2(delta_g)
    
    w = np.zeros((ny,nx))
    rms = np.inf
    i = 0
    while True:
        i += 1
        #~ print "Iteracion", i
        w_old, rms_old = w, rms
        w = root_calculation(k, w, w0, constant_term, order, k_cut)
        rms = rms_calculation(w_old, w)
        #~ print "RMS:", rms, "\n"
        if 0 <= rms_old - rms < rms_min_dif:
            print "Parker-Oldenburg completed for this iteration"
            print "Parker-Oldenburg RMS:", rms
            return w
        elif rms > rms_old:
            print "Parker Oldenburg has found an increasing RMS value"
            print "Parker-Oldenburg RMS:", rms
            return w_old
    

def anomaly_calculation(x, y, w, t, rho):
    from fatiando import utils
    from fatiando.mesher import Prism
    from fatiando.gravmag import prism

    ny, nx = np.shape(x)
    prismas = []
    for i in xrange(ny-1):
        for j in xrange(nx-1):
            prisma = Prism(y[i][j], y[i+1][j+1], x[i][j], x[i+1][j+1], t,
                        t + (w[i][j] + w[i+1][j] + w[i][j+1] + w[i+1][j+1])/4.)
            prismas.append(prisma)
    
    gz = -prism.gz(y.ravel(),x.ravel(),np.zeros(nx*ny),prismas,dens=rho)
    gz = utils.mgal2si(np.reshape(gz,(ny,nx)))
    return gz

