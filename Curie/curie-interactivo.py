#  Copyright 2015 Santiago R Soler
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
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


import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from fatiando.vis import mpl
from lib.basics import *
from lib.objects_functions import *
from fatiando import gridder


"""
-------------------------
| Apertura del Proyecto |
-------------------------
"""
proyect_name = "payunia"
grid_name = "curiePayunia.dat"
curie_filename = "curie.dat"
default_spectrum_points=70

grid_name = proyect_name + "/" + grid_name
curie_filename = proyect_name + "/" + curie_filename

print "Project Name:", proyect_name
print "Magnetic Anomaly file:", grid_name
print "Curie Depth out file:", curie_filename


# Pruebo si existe la carpeta del proyecto
try:
    os.stat(proyect_name)
except:
    os.mkdir(proyect_name)
    print "Put the grid file in '" + proyect_name + "' directory."
    sys.exit()
    
# Leo la grilla de anomalias ( shape = (ny,nx) )
#~ x,y,anomaly,shape = gridder.load_surfer(grid_name,fmt='ascii')
x,y,anomaly = np.loadtxt(grid_name,unpack=True)
nx = 254
ny = 300
shape = (ny,nx)
    
# Pruebo si existe el archivo de salida

try:
    os.stat(curie_filename)
    outfile_exists = True
    x_done_points = []
    y_done_points = []
    curie_file = open(curie_filename,'r')
    i = 0
    for line in curie_file:
        i+=1
        if i>2:
            aa,bb,cc,dd,ee,ff,gg,hh = line.split()
            x_done_points.append(float(aa))
            y_done_points.append(float(bb))
    curie_file.close()
except:
    outfile_exists = False
    

## Leo el archivo temporal (ultimo punto seleccionado)
tmp_filename = proyect_name + "/.tmp"
try:
    os.stat(tmp_filename)
    tmp_exists = True
    tmp_file = open(tmp_filename,'r')
    tmp_reading = tmp_file.read()
    tmp_file.close()
    x_last_point = tmp_reading.split()[0]
    y_last_point = tmp_reading.split()[1]
except:
    tmp_exists = False

"""
------------------------------
| Construccion de la Ventana |
------------------------------
"""
points2pad = raw_input("\nAmmount of spectrum points (default = "+
                    str(default_spectrum_points) + ", 0 means no padding): ")
if points2pad == "":
    points2pad = default_spectrum_points
else:
    points2pad = int(points2pad)
points2pad *= 2


# Seleccion de la ventana
bool_plot = False
while bool_plot == False:
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax2.set_aspect('equal')
    if tmp_exists == True:
        mpl.plot(x_last_point,y_last_point,'o',color='w')
    if outfile_exists == True:
        mpl.plot(x_done_points[:-1],y_done_points[:-1],'o',color='k')
        mpl.plot(x_done_points[-1],y_done_points[-1],'s',color='k')
    mpl.contourf(x,y,anomaly,shape,50,interp=False)
    mpl.colorbar()
    window = WindowConstruction(x,y,shape,anomaly,ax2,ax1,points2pad)
    mpl.show()
    
    if window.x_center != None:
        bool_plot = True
    
# Recupero los datos de la ventana
x_center, y_center = window.x_center, window.y_center
width = 2*window.half_width
print "(x,y)=",(x_center,y_center),"Width:",width

# Construccion de la ventana y el espectro de potencias
area = (window.x1, window.x2, window.y1, window.y2)
window_shape = window.shape_window
x_cut, y_cut, [anomaly_cut], shape_cut = cut_regular(x,y,[anomaly],shape,area)
n_padding = int((points2pad - window_shape[1])/2.0)
if n_padding < 0: n_padding = 0
x_pad, y_pad, anomaly_pad, shape_pad = padding(x_cut,y_cut,anomaly_cut,shape_cut,n_padding)
f, spectrum, log_spectrum, errors, log_errors = radial_spectrum(x_pad,y_pad,anomaly_pad,shape_pad)

## Actualizo el archivo temporal
tmp_file = open(tmp_filename,'w')
tmp_file.write(str(int(x_center)) + "\t" + str(int(y_center)))
tmp_file.close()

"""
---------------------
| Fiteo de la Curva |
---------------------
"""
# Seleccion de puntos para fitear
bool_fit = False
while bool_fit == False:
    bool_selection = False
    while bool_selection == False:
        fitting_area = AreaSelection(False,False)
        plt.errorbar(f,log_spectrum,yerr=log_errors,fmt='o')
        plt.ylabel(r"$\ln(\Phi_{\Delta T})$")
        plt.ylim(-max(abs(log_spectrum)),0.2)
        plt.xlim(min(f),max(f))
        plt.show()
        if fitting_area.x1 != None:
            bool_selection = True
    
    area = (fitting_area.x1, fitting_area.x2, fitting_area.y1, fitting_area.y2)
    f_fit, log_spectrum_fit, [log_errors_fit] =  inside_points(f, log_spectrum,
                                                           [log_errors], area)
    zt, zb, a, varianzas = curve_fitting(f_fit,log_spectrum_fit,log_errors_fit,
                                         f,log_spectrum,log_errors)
    if varianzas[1]/zb > 0.4:
        default = True
    else:
        default = False
        
    bool_refit = question("Would you like to redo the fitting?",default=default)
    if bool_refit == False:
        bool_fit = True

"""
-------------------------
| Guardado de los datos |
-------------------------
"""

bool_save = question("Would you like to save this messurement?",default=True)
if bool_save == True:
    if outfile_exists == False:
        curie_file = open(curie_filename,'w')
        curie_file.write("All quantities in meters\n")
        curie_file.write("x_center \t y_center \t Width \t Zt \t Zb \t Zt_err \t Zb_err \t Zb_rel_err \n")
        curie_file.close()

    curie_file = open(curie_filename,'a')
    
    zt_m, zb_m, varianzas_m = zt/(2*np.pi), zb/(2*np.pi), varianzas/(2*np.pi) 
    
    line = '{0:7d} \t {1:7d} \t {2:7d} \t {3:6d} \t {4:6d} \t {5:6d} \t {6:6d}\
            \t {7:.2f} \n'.format(int(x_center), int(y_center), int(width), 
            int(zt_m),int(zb_m),int(varianzas_m[0]), int(varianzas_m[1]), 
            varianzas_m[1]/zb_m) 
    
    curie_file.write(line)


