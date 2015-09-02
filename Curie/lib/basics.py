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

def inside_points(x,y,scalars,area):
    x1, x2, y1, y2 = area
    if x1 > x2:
        x1, x2 = x2, x1
    if y1 > y2:
        y1, y2 = y2, y1
    inside = [i for i in range(len(x))
              if x[i]>x1 and x[i]<x2 and y[i]>y1 and y[i]<y2]    
    
    #~ x_inside, y_inside = [], []
    #~ scalars_inside = [[] for j in range(len(scalars))]
    #~ for i in range(len(x)):
        #~ if x[i]>x1 and x[i]<x2 and y[i]>y1 and y[i]<y2:
            #~ x_inside.append(x[i])
            #~ y_inside.append(y[i])
            #~ for j in range(len(scalars)):
                #~ scalars_inside[j].append(scalars[j][i])
    return x[inside], y[inside], [s[inside] for s in scalars]


def nearest_point(xp,yp,x,y):
    """
    Gives the nearest point of the arrays x,y to the
    point (xp,yp).
    If we have a grid and we want to select a certain point,
    we can find it using this function by giving an near
    point (xp,yp) and calculating the nearest point of
    the grid.
    x must vary first and then y.
    Arguments:
    xp (float): x coordinate of the point
    yp (float): y coordinate of the point
    x (list): 1D array of x coordinates of the grid points.
    y (list): 1D array of y coordinates of the grid points.
    Returns:
    (x_nearest,y_nearest): tuple of the nearest grid points to (xp,yp).
    """
    x_nearest = x[0]
    y_nearest = y[0]
    distance_x = abs(xp-x_nearest)
    distance_y = abs(yp-y_nearest)
    for xi in x:
        if abs(xi-xp)<distance_x:
            x_nearest = xi
            distance_x = abs(xp-x_nearest)
    for yi in y:
        if abs(yi-yp)<distance_y:
            y_nearest = yi
            distance_y = abs(yp-y_nearest)
    return (x_nearest,y_nearest)
    
    
def rings_construction(window_length):
    """
    Calculate points of concentric rings.
    window_length (int) is the number of points per axe.
    It must be an odd number.
    
    """
    assert window_length%2 != 0, "window_length must be odd."
    window_length = int(window_length)
    center = [window_length/2, window_length/2]
    number_of_rings = window_length/2 + 1 
    rings = [[] for i in range(number_of_rings)]
    rings[0].append(center)

    for i in xrange(window_length):
        for j in xrange(window_length):
            distance = (i-center[0])**2 + (j-center[1])**2
            for n in xrange(1,number_of_rings):
                if distance >= (n-0.5)**2 and distance < (n+0.5)**2:
                    rings[n].append([i,j]) ## antes agregabamos el image[j][i], ahora solo los indices
                    break
    return rings


def radial_profile_w_errors(power_spectrum,rings):
    """
    Calculate the average of the spectrum for each ring
    Return a 1-d array containing the radialProfile
    """
    radialProfile = []
    errors = []
    
    for ring in rings:
        radial_average = sum([power_spectrum[x[0]][x[1]] for x in ring])/float(len(ring))
        error = np.sqrt(sum([(radial_average - power_spectrum[x[0]][x[1]])**2 for x in ring])/float(len(ring)))
        radialProfile.append( radial_average ) ## hago un promedio de las image[j][i]
        errors.append(error)
    
    return np.array(radialProfile),np.array(errors)

def question(sentence,default=True):
    answer = False
    yes = ['yes','y']
    no = ['no','n']
    if default == True:
        yes.append('')
        key_indicator = " [Y/n] "
    else:
        no.append('')
        key_indicator = " [y/N] "
    while answer == False:
        choice = raw_input(str(sentence) + key_indicator).lower()
        if choice in yes:
            answer == True
            return True
        elif choice in no:
            answer == True
            return False
