#!/usr/bin/env python2.7

#Code: plot_var.py (Python Script)
#Author: Nate Snook (CASA/CAPS/SoM)
#Written: Mar. 2008
#Last modified:  Jun. 17, 2008
#
#Purpose:
#   Plots a single ARPS variable field (specified in plot_config.py) 
#   in the x-y plane.

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import sys
from numpy import *
from subprocess import call
from mpl_toolkits.basemap import Basemap, shiftgrid
from PyNIO import Nio
from plot_config import *

#---------------------------------------------------------------------------------#
#       BEGINNING OF EXECUTABLE CODE                                              #
#---------------------------------------------------------------------------------#

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

if len(sys.argv) == 2:  #if the user entered a variable on the command line...
   var = sys.argv[1] #...then use it.

if len(sys.argv) == 3:  #if the user entered 2 arguments on the command line...
   var = sys.argv[1]       #...the first should be the variable...
   hdfhistory = sys.argv[2]  #...and the second should be the argument.

#Initialize reflectivity arrays
var_sfc_filtered = zeros((nx,ny))

#Open the history dump file for reading
dumpfile = Nio.open_file(hdfhistory, mode = 'r', options = None, history='', format='hdf')

#List the variables that are available
#varnames = dumpfile.variables.keys()
#print varnames

if dimensions == 3:
   try:
      #Read in the data for qr(k,j,i)
      var_3D = dumpfile.variables[var][:,:,:]
   except:
      if (var == 'ptprt') or (var == 'ptpert'): # If the user wants perturbation potential temperature (ptprt)...
         pt = dumpfile.variables['pt'][:,:,:]   #...then calculate it from potential temperature (pt) and ptbar.
         ptbar = dumpfile.variables['ptbar'][:,:,:]
         ptprt = zeros((len(pt[:]),len(pt[1][:]), len(pt[1][1][:])))
         ptprt[:,:,:] = pt[:,:,:] - ptbar[:,:,:]    #perturbation = actual value - average value
         var_3D = ptprt
      if (var == 'AGL'):
         sfc_height = dumpfile.variables['zpsoil'][0,:,:]  #...calculate it from zpsoil[sfc] and zp.
         zp = dumpfile.variables['zp'][:,:,:]
         agl = zp - sfc_height
         var_3D = agl
      else:
         print ' '
         print '*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
         print '*   Variable "' + var + '" is not supported.  Stop goofing around and try again!   *'
         print '*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
         print ' '
         print 'Valid options are:'
         call('./varlist.py')         
      #end if/else
   #end try/except
      
   #Get the surface values only
   var_sfc = var_3D[lev,:,:]

   #----------------------------------------
   # Output information on height of surface
   #----------------------------------------
   sfc_height = grdbasfile.variables['zpsoil'][0,:,:]  #...calculate it from zpsoil[sfc] and zp.
   zp = dumpfile.variables['zp'][:,:,:]
   agl = zp - sfc_height   
   print '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   print 'Plotting at grid level ' + str(lev + 1) + '.'
   print 'Level ' + str(lev + 1) + ' is between ' + str(agl[lev,:,:].min()) + ' and ' + str(agl[lev,:,:].max()) + 'm AGL.'
   print '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   
   #Read in wind for plotting, if necessary
   if windbarbs == 1:
      u = dumpfile.variables['u'][:,:,:]
      v = dumpfile.variables['v'][:,:,:]
      u_sfc = u[lev,:,:]
      v_sfc = v[lev,:,:]
      u_sfc_kts = u_sfc * 1.94384449
      v_sfc_kts = v_sfc * 1.94384449

elif dimensions == 2:
   var_sfc = dumpfile.variables[var][:,:]
elif dimensions == 1: 
   print '1-dimensional variables are not supported in this script'
else:
   print "I don't know what crazy world you're living in, but", dimensions, "dimensions is NOT a viable option here."
   exit

#Filter bad data...
if filt == 1:
   for ix in arange(0,nx):
      for jy in arange(0,ny):
         if (var_sfc[ix,jy] > filt_max or var_sfc[ix,jy] < filt_min):
            var_sfc_filtered[ix,jy] = -1.0
	 else:
            var_sfc_filtered[ix,jy] = var_sfc[ix,jy]
         #end if
      #end for
   #end for
#end if

print 'Plotting', var
print 'Maximum value = ' + str(var_sfc.max()) + '  ||  Minimum value = ' + str(var_sfc.min())
if filt == 1:
   print 'Filtered Max =', var_sfc_filtered.max(), ' ||  Filtered Min =', var_sfc_filtered.min()

#Create a figure
fig=plt.figure()

ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

x = arange(0,width_x+dx,dx)  #the '+dx' or '+dy' is so that arange() will include the final...
y = arange(0,width_y+dy,dy)  #...gridpoint at width_x and width_y.  Trust me, it's necessary.

x, y = meshgrid(x,y)

#Cyllindrical map projection containing states, counties, rivers, and urban boundaries
map = Basemap(projection='lcc', width=width_x, height=width_y, lat_1=trulat1, lat_2=trulat2, lat_0=ctrlat, lon_0=ctrlon, resolution='h', area_thresh=10.) #lambert conformal -- uses parameters from grdbas file!

#Define default thresholds for some common variables.  Feel free to change these as you see fit.
thresholds = arange(var_sfc.min(), var_sfc.max() + (var_sfc.max() - var_sfc.min())/6, (var_sfc.max() - var_sfc.min())/6)   #by default, six equal partitions.
if var == 'pt':    #Potential temperature contour thresholds
   print 'Variable is potential temperature (pt) -- contouring from 290 to 310K'
   thresholds = [290, 292, 294, 296, 298, 300, 302, 304, 306, 308, 310]
if color_plots == 0:
   colormap = gray_darktolight[::2]
else:
   colormap = None
#end if/else (color_plots)
   
#create a contour plot
if filt == 1:
   CF = map.contourf(x, y, var_sfc_filtered, thresholds, colors = colormap)
   CS = map.contour(x, y, var_sfc_filtered, thresholds, colors = 'black', linewidths = 1.5, linestyles = 'solid')
else:
#   CF = map.contourf(x, y, var_sfc, thresholds, colors = colormap)
#   CS = map.contour(x, y, var_sfc, thresholds, colors = 'black', linewidths = 1.5, linestyles = 'solid')

   CF = map.contourf(x, y, var_sfc, thresholds, colors = colormap)
#   CS = map.contour(x, y, var_sfc, colors = 'black', linewidths = 1.5, linestyles = 'solid')
#end if/else (filt)

if windbarbs == 1:
   #plot every {skip}th wind barb
   WND = map.barbs(x[::skip,::skip], y[::skip,::skip], u_sfc_kts[::skip,::skip], v_sfc_kts[::skip,::skip])
   print 'Windbarbs will be plotted at level ' + str(lev + 1) + '...'
#end if (windbarbs)

#add a color bar
CB = plt.colorbar(CF, shrink=0.8, extend='both')

#Gridpoint labels... it's sad that it's come to this.
gridpoint_labels = False  #Do we want gridpoint labels? True or False
if gridpoint_labels == True:
   for ix in arange(100, 200, 5):
      for jy in arange(100, 200, 5):
         plt.text((ix * dx), (jy * dy), '(' + str(ix) + ',' + str(jy) + ')', color = 'red', fontsize = 9)
      #end for
   #end for
#end if

#Specifics about the map
map.drawcoastlines()  #Draws coastlines... of course.
map.drawcountries()   #Draws political boundaries of countries.
# map.drawrivers(color='blue')  #Draws major rivers
map.drawstates(linewidth='1.5', color = "#333333")  #Draws political boundaries of states and provinces.

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
   
#Gridpoint labels... it's sad that it's come to this.
gridpoint_labels = 0  #Do we want gridpoint labels?  1 = yes, 0 = no.
if gridpoint_labels == 1:
   for ix in arange(147, nx-55, 1):
      for jy in arange(134, 188, 1):
         plt.text((ix * dx), (jy * dy), '(' + str(ix) + ',' + str(jy) + ')', color = 'gray', fontsize = 9)
      #end for
   #end for
#end if

map.drawmapboundary()
meridians = arange(-120., 80., 2)     #Sets up meridians from 120W to 80W at 2 degree increments.
map.drawmeridians(meridians, labels=[0,0,0,1])  #Draws meridians (lines of constant longitude)

parallels = arange(20., 48., 2)       #Sets up parallels from 20N to 46N at 2 degree increments.
map.drawparallels(parallels, labels=[1,0,0,0])    #Draws parallels (lines of constant latitude)

#set plot title
plt.title('Plot of ' + str(var) + ' in ' + str(hdfhistory) + ' -- Grid level : ' + str(lev))

#Overlay the arbitrary contour field on the map.
#map.contour(self, X, Y, Z)

plt.show()

