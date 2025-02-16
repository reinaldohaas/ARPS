#!/usr/bin/env python2.7

#Code: ensemble_spread.py (Python Script)
#Author: Nate Snook (CASA/CAPS/SoM)
#Written: Apr. 2008
#Last modified:  Jun. 17, 2008
#
#Purpose:
#   Calculates the RMS spread of a given ARPS variable field (defined in 
#   plot_config.py) and plots it on the x-y plane.

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import numpy as np
from numpy.numarray import *
from mpl_toolkits.basemap import Basemap #, shiftgrid
from PyNIO import Nio
from plot_config import *

#---------------------------------------------------------------------------------#
#                         BEGINNING OF EXECUTABLE CODE                            #
#---------------------------------------------------------------------------------#

#---------------------------------------------------#
#  Step 0.5)  Read in geographic data from grdbas   #
#---------------------------------------------------#

#Open the history file containing grid base data
grdbasfile = Nio.open_file(grdbas, mode = 'r', options = None, history='', format='hdf')

#Read in geographic data (center lat, lon; true lat, lon)
ctrlat = float(grdbasfile.attributes['ctrlat'])
ctrlon = float(grdbasfile.attributes['ctrlon'])
trulat1 = float(grdbasfile.attributes['trulat1'])
trulat2 = float(grdbasfile.attributes['trulat2'])
trulon = float(grdbasfile.attributes['trulon'])

#Read in grid data (nx, ny, nz, dx, dy)
nx = int(grdbasfile.attributes['nx'])
ny = int(grdbasfile.attributes['ny'])
nz = int(grdbasfile.attributes['nz'])
dx = float(grdbasfile.attributes['dx'])
dy = float(grdbasfile.attributes['dy'])

#Calculate derived grid data (width_x, width_y)
width_x = (nx - 1) * dx
width_y = (ny - 1) * dy

#Print out grid and geographic data for verification
print 'Finished reading geographic and grid data from grdbas file.'
print 'nx =', nx, ',', 'ny =', ny, ',', 'nz =', nz
print 'dx =', dx, ',', 'dy =', dy
print 'Domain is centered at ( ', ctrlon, ',', ctrlat, ')'
print 'True latitudes are', trulat1, 'and', trulat2
print 'True longitute is', trulon

#---------------------------------------------------#
#  Step 1)  Calculate the ensemble mean.            #
#---------------------------------------------------#

#Create a counter to keep track of how many ensemble members have been used.
membercount = 0

#Create arrays to store ensemble spread, deviation, and ensemble mean
var_spread = zeros((nx,ny))
var_mean = zeros((nx,ny))
dev = zeros((nx,ny))

#Loop over all ensemble members
for historydump in filelist:
   
   #Open the history dump file for reading
   dumpfile = Nio.open_file(historydump, mode = 'r', options = None, history='', format='hdf')
   
   #List the variables that are available 
   #varnames = dumpfile.variables.keys()
   #print varnames
   
   #Read in the data for qr(k,j,i)
   var3d = dumpfile.variables[var][:,:,:]
   print 'Calculating RMS Spread for ' + str(var) + '...'
   print 'global max = ', var3d[:,:,:].max()
   #Get the surface values only
   var_sfc = var3d[lev,:,:]
   #print qr_sfc
   print 'sfc max = ',var_sfc.max()
   print 'sfc min = ',var_sfc.min()
   
   #close the current file
   #Nio.close(history='')
   for ix in range(1,(nx - 1)):   #loop over points in the x-direction
      for jy in range(1,(ny - 1)):   #loop over points in the y-direction
         #sum up the individual members for calculation of an ensemble mean.
	 var_mean[ix,jy] = var_mean[ix,jy] + var_sfc[ix,jy]   
      #end for
   #end for
   
   membercount = membercount + 1    #increment our ensemble member counter.
   print 'Data for ensemble member ', membercount, ' read successfully.'
#end for

for ix in range(1,(nx - 1)):   #loop over points in the x-direction
   for jy in range(1,(ny - 1)):   #loop over points in the y-direction
      var_mean[ix,jy] = (var_mean[ix,jy] / membercount)   #divide by number of members to obtain ensemble mean.
   #end for
#end for

#---------------------------------------------------#
#  Step 2)  Calculate RMS ensemble spread.          #
#---------------------------------------------------#

# Remember -- RMS error is the square Root of the Mean of the Square of the deviations.

print 'calculating deviations...'

#Loop over all ensemble members
for historydump in filelist:
   
   #Open the history dump file for reading
   dumpfile = Nio.open_file(historydump, mode = 'r', options = None, history='', format='hdf')
   
   #Read in the data for qr(k,j,i)
   var3d = dumpfile.variables[var][:,:,:]
   var_sfc = var3d[lev,:,:]
   for ix in range(1,(nx - 1)):   #loop over points in the x-direction
      for jy in range(1,(ny - 1)):   #loop over points in the y-direction
         dev[ix,jy] = abs(var_sfc[ix,jy] - var_mean[ix,jy])  #calculate deviation for this member
	 var_spread[ix,jy] = var_spread[ix,jy] + dev[ix,jy] * dev[ix,jy]  #sum the squares of the deviations
      #end for
   #end for
#end for

print 'calculating RMS spread...'
for ix in range(1,(nx - 1)):   #loop over points in the x-direction
   for jy in range(1,(ny - 1)):   #loop over points in the y-direction
      var_spread[ix,jy] = var_spread[ix,jy] / membercount  #divide by number of members
      var_spread[ix,jy] = sqrt(var_spread[ix,jy])   #...and get the final value by taking the square root.
   #end for
#end for

print '     Maximum RMS spread', var_spread.max()
print '     Minimum RMS spread', var_spread.min()

#---------------------------------------------------#
#  Step 3)  Plot a figure!                          #
#---------------------------------------------------#

print 'creating the basemap figure...'

#Create a figure
fig=plt.figure()

ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

x = np.arange(0,width_x+1,dx)  #the '+1' is so that arange() will include the final gridpoint...
y = np.arange(0,width_y+1,dy)  #...at width_x and width_y.  Trust me, it's necessary.
                               #If you're using grid spacing less than 1 meter, change it to '+dx'.

x, y = np.meshgrid(x,y)

#Lambert conformal map projection containing states, counties, rivers, and urban boundaries
map = Basemap(projection='lcc', width=width_x, height=width_y, lat_1=trulat1, lat_2=trulat2, lat_0=ctrlat, lon_0=ctrlon, resolution='h', area_thresh=10.) #lambert conformal -- should match parameters in arps.input!

thresholds = np.arange(0.0, var_spread.max()*1.15 , round(var_spread.max() * 0.1, 5))

#create a contour plot
CS = map.contourf(x, y, var_spread, thresholds, linewidths = 1.5)
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
   map.readshapefile(shapefile_dir + '/county/countyp020','counties',drawbounds=True, linewidth=0.5, color='#333333')  #Draws US county boundaries.

if draw_urban == 1:
   print 'reading urban boundaries from shapefile'
   if color_plots == 1:
      map.readshapefile(shapefile_dir + '/urban/geo_2000_d_dec_00000_400_py_z1', 'urban areas', drawbounds=True, linewidth=0.5, color='purple')  #Draws urban areas in purple.
   else:
      map.readshapefile(shapefile_dir + '/urban/geo_2000_d_dec_00000_400_py_z1', 'urban areas', drawbounds=True, linewidth=0.5, color='#333333')  #Draws urban areas in gray.
   #end if
   
if draw_radar == 1:
   print 'reading radar range rings from shapefile'
   map.readshapefile(shapefile_dir + '/CASA/CASA_radars_40km_range', 'radars', drawbounds=True, linewidth=1.0, color='black')  #Draws radar range rings in black, bold lines

map.drawmapboundary()

meridians = np.arange(-120., 80., 2)     #Sets up meridians from 110W to 90W at 2 degree increments.
map.drawmeridians(meridians, labels=[0,0,0,1])  #Draws meridians (lines of constant longitude)

parallels = np.arange(20., 50., 2)       #Sets up parallels from 30N to 40N at 2 degree increments.
map.drawparallels(parallels, labels=[1,0,0,0])    #Draws parallels (lines of constant latitude)

#set plot title
plt.title('40 Member Ensemble -- RMS Ensemble Spread of ' + str(var) + ' -- ' + str(arpstime))

#Overlay the arbitrary contour field on the map.
#map.contour(self, X, Y, Z)

print 'rendering plot...'

plt.show()

