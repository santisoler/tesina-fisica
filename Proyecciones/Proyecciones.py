import pyproj
import numpy as np

def gauss_kruger_argentina(lon, lat, faja, mts2grad=False):
    """
    proj=tmerc  +lat_0=Latitude of natural origin 
                +lon_0=Longitude of natural origin
                +k=Scale factor at natural origin 
                +x_0=False Easting
                +y_0=False Northing
                
    Faja:       Validez:
    1           west of 70.5 W
    2           70.5 W to 67.5 W
    3           67.5 W to 64.5 W
    4           64.5 W to 61.5 W
    5           61.5 W to 58.5 W
    6           58.5 W to 55.5 W
    7           east of 55.5 W
    """
    

    argentina = {'faja1': {'lat_0': -90, 'lon_0': -72, 'k_0': 1.0, 'x_0': 1500000, 'y_0': 0},
                 'faja2': {'lat_0': -90, 'lon_0': -69, 'k_0': 1.0, 'x_0': 2500000, 'y_0': 0},
                 'faja3': {'lat_0': -90, 'lon_0': -66, 'k_0': 1.0, 'x_0': 3500000, 'y_0': 0},
                 'faja4': {'lat_0': -90, 'lon_0': -63, 'k_0': 1.0, 'x_0': 4500000, 'y_0': 0},
                 'faja5': {'lat_0': -90, 'lon_0': -60, 'k_0': 1.0, 'x_0': 5500000, 'y_0': 0},
                 'faja6': {'lat_0': -90, 'lon_0': -57, 'k_0': 1.0, 'x_0': 6500000, 'y_0': 0},
                 'faja7': {'lat_0': -90, 'lon_0': -54, 'k_0': 1.0, 'x_0': 7500000, 'y_0': 0}}
            
    assert faja in argentina, "Faja no valida."
            
    
    lat_0 = argentina[faja]["lat_0"]
    lon_0 = argentina[faja]["lon_0"]
    k_0 = argentina[faja]["k_0"]
    x_0 = argentina[faja]["x_0"]
    y_0 = argentina[faja]["y_0"]
    
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)
    proj = pyproj.Proj(proj="tmerc", lat_0= lat_0, lon_0=lon_0, x_0=x_0,
                       y_0=y_0, k_0=k_0)
    x, y = proj(lon, lat, inverse=mts2grad)
    return x, y



def UTM(lon, lat, zone, hemisphere, mts2deg=False):
    """
    Calculates the projection of lon, lat points to UTM x,y coordinates.
    Also can do the inverse mts2deg.    
    """
    zone = int(zone)
    hemisphere = hemisphere.lower()
    h_list = ["north", "norte", "south", "sur"]
    assert hemisphere in h_list, "Hemisferio no valido (norte o sur)."
            
    proj_argument = '+proj=utm +zone=' + str(zone)
    if hemisphere == "north" or hemisphere == "norte":
        proj_argument += '+north'
    else:
        proj_argument += '+south'
    
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)
    proj = pyproj.Proj(proj_argument)
    x, y = proj(lon, lat, inverse=mts2deg)
    return x, y



def UTM_zone_calculator(lon):
    """
    Calculates the UTM zone for a given latitude.
    """
    if abs(lon) > 180:
        lon = lon - 360
    zones = [[-180 + 6*i,-174 + 6*i] for i in range(60)]
    for i in range(len(zones)):
        if zones[i][0] <= lon < zones[i][1]:
            zone = i+1
            return zone

