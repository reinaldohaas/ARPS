#!/usr/bin/env python2.7

#Code: probability_overlay_neighborhood_prcrate.py (Python Script)
#Author: Nate Snook (CASA/CAPS/SoM)
#Written: Mar. 2008
#Last modified:  Jun. 17, 2008
#
#Purpose:
#   Plots a difference for a given ARPS variable field (specified in plot_config.py)
#   between two ARPS history dump files.

from PyNIO import Nio
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap 
from plot_config import *

file1 = Nio.open_file(hdfhistory, mode = 'r', options = None, history='', format='hdf')
file2 = Nio.open_file(hdfhistory2, mode = 'r', options = None, history='', format='hdf')

#varnames = dumpfile.variables.keys()
#print varnames

#read in pt from the two files.
pt1 = file1.variables['pt'][:,:,:]
print 'File #1:'
print 'Min: ', pt1[:,:,:].min(), '      Max: ',pt1[:,:,:].max()

pt2 = file2.variables['pt'][:,:,:]
print 'File #2:'
print 'Min: ', pt2[:,:,:].min(), '      Max: ',pt2[:,:,:].max()

#Calculate difference between pt1 and pt2

ptdiff = pt2 - pt1
print '-----------------------------------------------'
print 'Difference field:'
print 'Minimum difference: ', ptdiff[:,:,:].min(), '     Maximum difference: ', ptdiff[:,:,:].max()
print '-----------------------------------------------'

sfc_pt_diff = ptdiff[0,:,:]
#print pt_sfc
#print pt_sfc.max()

fig=plt.figure()

ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

x = np.arange(0,518000,2000)
y = np.arange(0,518000,2000)

x, y = np.meshgrid(x,y)

#Lambert conformal map projection containing states, counties, rivers, and urban boundaries
map = Basemap(projection='lcc', width=516000, height=516000, lat_1=30.0, lat_2=38.0, lat_0=34.25, lon_0=-99.0, resolution='h', area_thresh=10.) #lambert conformal

thresholds = np.arange(-3.0, 3.0, 0.25)

#create a contour plot
CS = map.contour(x, y, sfc_pt_diff, thresholds, linewidths = 1.5)
#plt.clabel(CS, fontsize=9, inline=1)  #labels the contours

#add a color bar
CB = plt.colorbar(CS, shrink=0.8, extend='both')

#Specifics about the map
map.drawcoastlines()  #Draws coastlines... of course.
map.drawcountries()   #Draws political boundaries of countries.
map.drawrivers(color='blue')  #Draws major rivers
map.drawstates(linewidth='1.5')  #Draws political boundaries of states and provinces.

if draw_counties == 1:
   print 'reading counties from shapefile'
   map.readshapefile('/home/nsnook/shapefiles/county/countyp020','counties',drawbounds=True, linewidth=0.5, color='gray')  #Draws US county boundaries.

if draw_urban == 1:
   print 'reading urban boundaries from shapefile'
   map.readshapefile('/home/nsnook/shapefiles/urban/geo_2000_d_dec_00000_400_py_z1', 'urban areas', drawbounds=True, linewidth=0.5, color='purple')  #Draws urban areas in purple.

if draw_radar == 1:
   print 'reading radar range rings from shapefile'
   map.readshapefile('/home/nsnook/shapefiles/CASA/AllRadarRings', 'radars', drawbounds=True, linewidth=1.5, color='black')  #Draws radar range rings in black, bold lines

map.drawmapboundary()

meridians = np.arange(-110., 90., 2)     #Sets up meridians from 110W to 90W at 2 degree increments.
map.drawmeridians(meridians, labels=[0,0,0,1])  #Draws meridians (lines of constant longitude)

parallels = np.arange(30., 40., 2)       #Sets up parallels from 30N to 40N at 2 degree increments.
map.drawparallels(parallels, labels=[1,0,0,0])    #Draws parallels (lines of constant latitude)

#set plot title
plt.title('Difference in sfc_pt between members 001 and 003 at 0105Z, 9 May 2007')

plt.show()
