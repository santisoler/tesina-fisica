#  Copyright 2015 Santiago R Soler
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  

import numpy as np
from scipy import fftpack
from matplotlib.patches import Rectangle
from lib.basics import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class WindowConstruction(object):
    ## Shape = (ny,nx)
    def __init__(self,x,y,shape,anomaly,ax,ax2,points2pad):
        self.x = x
        self.y = y
        self.shape = shape
        self.ax = ax
        self.ax2 = ax2
        self.points2pad = points2pad
        self.anomaly = anomaly
        self.dx = (max(x)-min(x))/(self.shape[1]-1)
        self.half_width = (min(self.shape)/16)*self.dx
        self.x_center = None
        self.y_center = None
        self.rect = Rectangle((0,0), 1, 1,fc='None')
        self.x1 = None
        self.y1 = None
        self.x2 = None
        self.y2 = None
        self.shape_window = None
        self.l, = self.ax.plot([self.x_center],[self.y_center],'o')
        self.ax.add_patch(self.rect)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.ax.figure.canvas.mpl_connect('key_press_event', self.on_key)
        
        print "\nINSTRUCTIONS:"
        print "Click to select the window center"
        print "Move the center with arrows or click again"
        print "Resize the window with the mouse scroll or with '+' and '-'"
        print "Press 'i' to show information about the window"
        print "Press Enter or Right Click to plot the spectrum of the current window\n"
        
        
    def on_press(self,event):
        if event.inaxes!=self.ax:
            return
        if event.button == 1:
            self.x_center = event.xdata
            self.y_center = event.ydata
            self.x_center, self.y_center = nearest_point(self.x_center,
                                            self.y_center, self.x,self.y)
        elif event.button == 3:
            self.plot_spectrum()
        self.rectangle_construction()
        
        
    def on_scroll(self,event):
        self.half_width += event.step*self.dx
        self.rectangle_construction()
    
    
    def on_key(self,event):
        event_list = ["right","left","up","down"]
        if event.key in event_list:
            if event.key == "right":
                self.x_center += self.dx
            elif event.key == "left":
                self.x_center -= self.dx
            elif event.key == "up":
                self.y_center += self.dx
            elif event.key == "down":
                self.y_center -= self.dx
            self.rectangle_construction()
        if event.key == "enter":
            self.plot_spectrum()
        if event.key == "i":
            print "(x,y)=",(self.x_center,self.y_center),"Width:",self.half_width*2
        if event.key == "+" or event.key == "-":
            if event.key == "+":
                self.half_width += self.dx
            elif event.key == "-":
                self.half_width -= self.dx
            self.rectangle_construction()
            
            
    def rectangle_construction(self):
        self.x1 = self.x_center - self.half_width
        self.x2 = self.x_center + self.half_width
        self.y1 = self.y_center - self.half_width
        self.y2 = self.y_center + self.half_width
        self.shape_window = (2*self.half_width/self.dx+1,
                            2*self.half_width/self.dx+1)
        self.rect.set_width(self.x2 - self.x1)
        self.rect.set_height(self.y2 - self.y1)
        self.rect.set_xy((self.x1, self.y1))
        self.l.set_xdata([self.x_center])
        self.l.set_ydata([self.y_center])
        self.ax.figure.canvas.draw()
        

    def plot_radial_spectrum(self,f,log_spectrum,log_errors):
        self.ax2.cla()
        self.ax2.errorbar(f,log_spectrum,yerr=log_errors,fmt='o')
        self.ax2.set_ylabel(r"$\ln(\Phi_{\Delta T})$")
        self.ax2.set_ylim(-max(abs(log_spectrum)),0.2)
        self.ax2.set_xlim(0,max(f))
        self.ax2.figure.canvas.draw()
    
    
    def plot_spectrum(self):
        area = (self.x1, self.x2, self.y1, self.y2)
        x_cut, y_cut, [anomaly_cut], shape_cut = cut_regular(self.x,self.y,[self.anomaly],self.shape,area)
        n_padding = int((self.points2pad - self.shape_window[1])/2.0)
        if n_padding < 0: n_padding = 0
        x_pad, y_pad, anomaly_pad, shape_pad = padding(x_cut,y_cut,anomaly_cut,shape_cut,n_padding)
        f, spectrum, log_spectrum, errors, log_errors = radial_spectrum(x_pad,y_pad,anomaly_pad,shape_pad)
        self.plot_radial_spectrum(f,log_spectrum,log_errors)

        
        
def cut_regular(x, y, scalars, shape, area):
    """
    Return a subsection of a regular grid.

    The returned subsection is not a copy! In technical terms, returns a slice
    of the numpy arrays. So changes made to the subsection reflect on the
    original grid. Use numpy.copy to make copies of the subsections and avoid
    this.

    Parameters:

    * x, y
        Arrays with the x and y coordinates of the data points.
    * scalars
        List of arrays with the scalar values assigned to the grid points.
    * shape: tuple(ny,nx)
        Shape of the original data points.
    * area
        ``(x1, x2, y1, y2)``: Borders of the subsection

    Returns:

    * ``[subx, suby, subscalars,shape_window]``
        Arrays with x and y coordinates and scalar values of the subsection.
    shape_window: tuple(ny,nx)
        Shape of the cutted grid.

    """
    xmin, xmax, ymin, ymax = area
    if xmin > xmax:
        xmin, xmax = xmax, xmin
    if ymin > ymax:
        ymin, ymax = ymax, ymin

    # Calculates the distance between neighbour points
    if shape[0]*shape[1] != len(x):
        raise ValueError("The shape isn't the shape of the x,y grid.")

    distance_x = abs(max(x)-min(x)) / (shape[1]-1)
    distance_y = abs(max(y)-min(y)) / (shape[0]-1)
    shape_cutted = (int(abs(ymax-ymin)/distance_y + 1),
                    int(abs(xmax-xmin)/distance_x + 1))

    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    inside = [i for i in xrange(len(x))
              if x[i] >= xmin and x[i] <= xmax
              and y[i] >= ymin and y[i] <= ymax]

    return [x[inside], y[inside], [s[inside] for s in scalars], shape_cutted]


def padding(x,y,z,shape,n):
    # padding del z
    n = int(n)
    z = np.reshape(z, shape)
    z_pad = np.pad(z,n,'constant',constant_values=0)
    shape_pad = np.shape(z_pad)
    z_pad = np.reshape(z_pad,shape_pad)
    
    # padding del x,y
    dx = abs(max(x)-min(x))/(shape[1]-1)
    xx = np.linspace(min(x)-n*dx,max(x)+n*dx,shape[1]+2*n)
    yy = np.linspace(min(y)-n*dx,max(y)+n*dx,shape[1]+2*n)
    x_pad, y_pad = np.meshgrid(xx,yy)
    x_pad = np.reshape(x_pad,shape_pad[0]*shape_pad[1])
    y_pad = np.reshape(y_pad,shape_pad[0]*shape_pad[1])
    return x_pad, y_pad, z_pad, shape_pad


def radial_spectrum(x,y,anomaly,shape):
    nodes_distance = (max(x)-min(x))/(shape[1]-1)
    anomaly_2d = np.reshape(anomaly,shape)
    anomaly_fft = fftpack.fftshift(fftpack.fft2(anomaly_2d))
    power_spectrum = abs(anomaly_fft)**2
    rings = rings_construction(shape[0])
    spectrum, errors = radial_profile_w_errors(power_spectrum,rings)
    spectrum_max = max(spectrum)
    spectrum = spectrum/spectrum_max
    errors = errors/spectrum_max
    log_spectrum = np.log(spectrum)
    log_errors = abs(errors/spectrum)
    f = fftpack.fftfreq(shape[0],nodes_distance)
    f = f[:len(f)/2+1]
    return f, spectrum, log_spectrum, errors, log_errors


class AreaSelection(object):
    """
    This objects allows us to select a rectangular area from a plot.
    
    Arguments:
	verbose_points (bool): shows the points coordinates.
	verbose_size (bool): shows the size of the area.
    Properties:
	x1 (float): x coord of the first point
	x2 (float): y coord of the first point
	y1 (float): x coord of the second point
	y2 (float): y coord of the second point
	
    Usage:
	import matplotlib.pyplot as plt
	
	# Plot (you must draw the rectangle)
	area = AreaSelection(True,True)
	plt.plot()
	
	x1 = area.x1
	y1 = area.y1
	x2 = area.x2
	y2 = area.y2
	
	print x1,y1,x2,y2
	
    """
    
    def __init__(self,verbose_points=False,verbose_size=True):
        self.verbose_points = verbose_points
        self.verbose_size = verbose_size
        self.ax = plt.gca()
        self.rect = Rectangle((0,0), 1, 1,fc='r',alpha=0.5)
        self.x1 = None
        self.y1 = None
        self.x2 = None
        self.y2 = None
        self.onpress = False
        self.ax.add_patch(self.rect)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)

    def on_press(self, event):
        if event.inaxes!=self.ax:
            return
        self.onpress = True
        self.x1 = event.xdata
        self.y1 = event.ydata
        if self.verbose_points==True:
            print '(x1,y1):',(self.x1,self.y1)
	    
    def on_motion(self,event):
        if event.inaxes!=self.ax:
            return
        if self.onpress==True:
            self.x2 = event.xdata
            self.y2 = event.ydata
            self.rect.set_width(self.x2 - self.x1)
            self.rect.set_height(self.y2 - self.y1)
            self.rect.set_xy((self.x1, self.y1))
            self.ax.figure.canvas.draw()

    def on_release(self, event):
        if event.inaxes!=self.ax:
            return
        self.onpress = False
        self.x2 = event.xdata
        self.y2 = event.ydata
        self.rect.set_width(self.x2 - self.x1)
        self.rect.set_height(self.y2 - self.y1)
        self.rect.set_xy((self.x1, self.y1))
        self.ax.figure.canvas.draw()
        if self.verbose_points==True:
            print '(x2,y2):',(self.x2,self.y2)
        if self.verbose_size==True:
            print "Area Size: ",abs(self.x2-self.x1),"x", abs(self.y2-self.y1)


def curve_fitting(f_fit,log_spectrum_fit,log_errors_fit,f,log_spectrum,log_errors):
    
    def log_spectrum_function(k,zt,zb,a):
        return a - 2*k*zt + 2*np.log(np.ones(len(k))- np.exp(-k*(zb-zt)))    
    
    res = curve_fit(log_spectrum_function,f_fit,log_spectrum_fit,
        p0=[1000,2500,-1],sigma=log_errors_fit,maxfev=10000*(len(f_fit)+1))
    
    zt = res[0][0]
    zb = res[0][1]
    a = res[0][2]
    varianzas = np.sqrt(np.diag(res[1]))
    
    zt_km = zt/1000.0/(2*np.pi)
    zb_km = zb/1000.0/(2*np.pi)
    varianza_zt = varianzas[0]/1000.0/(2*np.pi)
    varianza_zb = varianzas[1]/1000.0/(2*np.pi)
    varianza_a = varianzas[2]/1000.0/(2*np.pi)
    
    print "\nFitting Results [km]:"
    print "Zt=",zt_km,"+/-",varianza_zt
    print "Zb=",zb_km,"+/-",varianza_zb
    print "Relative error of Zb=",varianza_zb/zb_km, "\n"
    
    if varianza_zb/zb_km > 0.4:
        print "VERY HIGH ERROR!\n"
                
    f_curve = np.linspace(0,max(f_fit),100)
    plt.plot(f_curve[1:],log_spectrum_function(f_curve[1:],zt,zb,a))
    plt.plot(f,log_spectrum,'.')
    plt.errorbar(f_fit,log_spectrum_fit,yerr=log_errors_fit,fmt='o',color='r')
    plt.show()            
    
    return zt,zb,a,varianzas
