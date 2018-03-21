# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap





lon1, lon2 = -69.1, -67.3
lat1, lat2 = -30.9, -27.7


## ------------------------------------------------------------------------
## ARGENTINA GRANDE
## BASEMAP
lon_arg1, lon_arg2 = -72, -55
lat_arg1, lat_arg2 = -40, -21


#~ rect = Rectangle((lon1, lat1), abs(lon2-lon1), abs(lat2-lat1))

#~ bm = Basemap(projection='tmerc', resolution='l', ellps='WGS84',
             #~ lon_0=-69, lat_0=-90, k_0=1,
             #~ llcrnrlon=lon_arg1, llcrnrlat=lat_arg1,
             #~ urcrnrlon=lon_arg2, urcrnrlat=lat_arg2)
bm = Basemap(projection='tmerc', resolution='i', ellps='WGS84',
             lon_0=-69, lat_0=-90, k_0=1,
             llcrnrlon=lon_arg1, llcrnrlat=lat_arg1,
             urcrnrlon=lon_arg2, urcrnrlat=lat_arg2)


bm.drawcoastlines(linewidth=1)
bm.drawcountries(linewidth=2)
bm.drawstates(linewidth=1, color='#3b3b3b')
#~ plt.gca().add_patch(rect)

bm.plot([lon1, lon1, lon2, lon2, lon1], [lat1, lat2, lat2, lat1, lat1], latlon=True, linestyle='--', color='r', linewidth=2)
plt.savefig("argentina-grande-deg.svg",  bbox_inches="tight")
plt.show()
