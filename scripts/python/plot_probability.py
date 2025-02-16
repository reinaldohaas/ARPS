#!/usr/bin/env python2.7

#Code: plot_probability.py (Python Script)
#Author: Nate Snook (CASA/CAPS/SoM)
#Written: Jun. 2008
#Last modified:  Jun. 17, 2008
#
#Purpose:
#   Creates a plot of probability of a given variable (defined in plot_config.py) 
#   greater than a given threshold (also defined in plot_config.py) using an ensemble
#   neighborhood method (neighborhood radius is also defined in plot_config.py).

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from numpy import *
from math import sqrt
from mpl_toolkits.basemap import Basemap #, shiftgrid
from PyNIO import Nio
from plot_config import *

#---------------------------------------------------------------------------------#
#                         BEGINNING OF EXECUTABLE CODE                            #
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
deltax = float(grdbasfile.attributes['dx'])
deltay = float(grdbasfile.attributes['dy'])

#Calculate derived grid data (width_x, width_y)
width_x = (nx - 1) * deltax
width_y = (ny - 1) * deltay

#Create a counter to keep track of how many ensemble members have been used...
membercount = 0
#...and another one to keep track of how many gridpoints are in our neighborhood.
neighborcount = 0

#We want to speed up our code, so we'll limit the area over which we search to
#determine if a point is in the neighborhood.  We'll do this by specifying a
#search radius.
searchradius_x = int(round((neighborhood_radius/deltax), 0))
searchradius_y = int(round((neighborhood_radius/deltay), 0))

#Set up a 'map' of the neighborhood to speed up execution
neighborhood_map = zeros((searchradius_x + searchradius_x + 1, searchradius_y + searchradius_y + 1))

for i in xrange(0, searchradius_x + searchradius_x + 1):
   xpts = searchradius_x - i    #determine how many points from the center we are in the x-dir.
   for j in xrange(0, searchradius_y + searchradius_y + 1):
      ypts = searchradius_y - j #determine how many points from the center we are in the y-dir.
      distance = sqrt( (((abs(ypts))*deltay) * (abs(ypts))*deltay) + 
                       (((abs(xpts))*deltax) * (abs(xpts))*deltax) )
      neighborhood_map[i, j] = distance
   #end for
#end for

#Check to make sure that the variable chosen is one that the script is setup to handle.
#as of the lastest update, supported variables include:
#       w  (vertical velocity, m/s)
#       qh (hail mixing ratio, unitless (g/g))
#       qr (rain mixing ratio, unitless (g/g))
#       prcrate4 (total precipitation rate (kg/(m^2)s))

#  Vertical vorticity and radar reflectivity are currently handled as special cases
#  by their own scripts as they are not directly stored in ARPS history dump files.

if var == 'w' or var == 'qh' or var == 'qr' or var == 'prcrate4':
   pass
else:
   print 'The variable', var, 'is not among the valid choices currently supported.'
   print 'You will need to edit plot_config.py and plot_probability.py to include it.'
   exit

#---------------------------------------------------------------------------------#
#                      PARAMETER VERIFICATION PRINT BLOCK                         #
#---------------------------------------------------------------------------------#
#Print out grid and geographic data for verification
print '-----------------------------------------------------------'
print 'Finished reading geographic and grid data from grdbas file.'
print 'nx =', nx, ',', 'ny =', ny, ',', 'nz =', nz
print 'dx =', deltax, ',', 'dy =', deltay
print 'Domain is centered at ( ', ctrlon, ',', ctrlat, ')'
print 'True latitudes are', trulat1, 'and', trulat2
print 'True longitute is', trulon

print '-----------------------------------------------------------'
print 'Running neighborhood verification for', var, '.'
print ' '
print 'Using a neighborhood radius of', neighborhood_radius, 'meters for fuzzy verification.'
print 'This requires a search radius of', searchradius_x, 'points in the x-direction...'
print '...and', searchradius_y, 'points in the y-direction.'
print 'Now beginning calculations...'
print '-----------------------------------------------------------'
#---------------------------------------------------------------------------------#

#Create an array to store probability
prob_field = zeros((nx,ny))
field_sfc = zeros((nx,ny))
field_sf = zeros((nx,ny))

#A temporary variable is needed to store probability at each point as it is calculated
current_prob = 0.

#Loop over all ensemble members
for historydump in filelist:
   
   #Open the history dump file for reading
   dumpfile = Nio.open_file(historydump, mode = 'r', options = None, history='', format='hdf')
   
   #List the variables that are available 
   #varnames = dumpfile.variables.keys()
   #print varnames
   
   print 'Now reading data from', historydump

# The following 'if' block handles reading in data for various variable choices.
   if var == 'prcrate4':
      #Read in the data for prcrate4(k,j,i) -- it is *already* a surface variable!
      field_sf = dumpfile.variables[var][:,:]
      print 'sfc max =', field_sfc.max(), '      sfc min = ', field_sfc.min()
      
      var_threshold = prcrate_threshold  #set the appropriate thresholds for plotting
      var_units = 'kg/(m^2)s'   #and note what the units are.
   
   elif var == 'qh':
      #Read in the data for qh(k,j,i)
      if (membercount + 1 <= num_lin or membercount + 1 > num_lin + num_wsm):
         field = dumpfile.variables['qh'][:,:,:]
      else:
         field = dumpfile.variables['qg'][:,:,:]
      print 'global max = ', qh[:,:,:].max()
      
      #Get the surface values only
      field_sf = field[0,:,:]
      #print qh_sfc
      print 'sfc max = ',field_sf.max(), '       sfc min = ',field_sf.min()
      
      var_threshold = qh_threshold  #set the appropriate thresholds for plotting
      var_units = 'g/g'
   
   else: # (for w and qr)
      #Read in the data for var(k,j,i)
      field = dumpfile.variables[var][:,:,:]
      print 'global max = ', field[:,:,:].max()
      
      #Get the surface values only
      field_sf = field[0,:,:]
      #print var_sfc
      print 'sfc max = ',field_sf.max(), '       sfc min = ',field_sf.min()
      
      if var == 'qr':
         var_threshold = qr_threshold
         var_units = 'g/g'
      #end if
      if var == 'w': 
         var_threshold = w_threshold
	 var_units = 'm/s'
	 print 'threshold and units set!'
      #end if
   
# If filtering is turned on, filter the data.
   if filt == 1:
      print 'Filtering is ACTIVE... now filtering data...'
      for ix in arange(0,nx):
         for jy in arange(0,ny):
            if field_sf[ix,jy] > filt_max or field_sf[ix,jy] < filt_min:
               field_sfc[ix,jy] = filtered
	    else:
               field_sfc[ix,jy] = field_sf[ix,jy]
            #end if
         #end for
      #end for
      print 'Filtering complete.'
   else:
      field_sfc = field_sf
   #end if/else
   
   #close the current file
   #Nio.close(history='')
   
   #Now, loop over the values in the variable array, getting a probability for each point
   #based on the surrounding points in its neighborhood.  Neighborhoods are circular with 
   #a radius defined by 'neighborhood_radius'.
   
   #Remember:  In Python, range(a, b, c) gives all numbers from a to b-1 with a step of c.
   #           For example, range(4, 8, 1) gives [4, 5, 6, 7].  Keep this in mind!
   
   for ix in xrange(0,nx):   #loop over points in the x-direction
      for jy in xrange(0,ny):   #loop over points in the y-direction
	 #Now, begin neighborhood calculations.  Iterate over each point within the search 
	 #radius to determine if it is in the neighborhood.  If it is, then verify a '1' if
	 #threshold is met, and '0' otherwise.
	 
         #reset the current probability and neighborhood counter to 0 for use with this point.
         neighborcount = 0
	 current_prob = 0
	 
	 for i in xrange((ix - searchradius_x),(ix + searchradius_x + 1),1):
	    for j in xrange((jy - searchradius_y),(jy + searchradius_y + 1),1):
	        if ((i < 0) or (i > nx - 1)):
		   continue  #We don't want to include points outside the domain.
		#end if
	        if ((j < 0) or (j > ny - 1)):
		   continue  #No points outside the domain in the y-direction either.
		#end if
	        
		#Go into the neighborhood map and see if our gridpoint is within the radius.
		coordx = searchradius_x - (ix - i)
		coordy = searchradius_y - (j - jy)
		if(neighborhood_map[coordx, coordy] <= neighborhood_radius):
		   neighborcount = neighborcount + 1      #This point is in the neighborhood.
	 	   if(field_sfc[i, j] > var_threshold):    #If field exceeds our threshold...
	    	      current_prob = (current_prob + 1)   #...then increment probability by 1.
		   #end if
		#end if
            #end for loop (j)
         #end for loop (i)
	 current_prob = (current_prob / neighborcount)   #Now we have a probability for this member.
	 prob_field[ix,jy] = prob_field[ix,jy] + current_prob   #Add this member's probability value to the accumulator.
	 
	 #The following is for debugging and can be omitted if desired
         if((i == 10 + searchradius_x) and (j == 10 + searchradius_y)):
            print 'The average neighborhood for this member contains', neighborcount, 'gridpoints.'
	 #end if
      #end for loop (jy)
   #end for loop (ix)
   
   membercount = membercount + 1    #increment our ensemble member counter.
   print 'Calculations for ensemble member ', membercount, ' completed.'
   print '=--=--=--=--=--=--=--=--=--=--='
#end for

prob_field = (prob_field / membercount)   #divide by number of members to obtain ensemble probability.

#Create a figure
fig=plt.figure()

ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

x = arange(0,width_x+deltax,deltax)  #the '+dx' or '+dy' is so that arange() will include the final...
y = arange(0,width_y+deltay,deltay)  #...gridpoint at width_x and width_y.  Trust me, it's necessary.

x, y = meshgrid(x,y)

#Lambert conformal map projection containing states, counties, rivers, and urban boundaries
map = Basemap(projection='lcc', width=width_x, height=width_y, lat_1=trulat1, lat_2=trulat2, lat_0=ctrlat, lon_0=ctrlon, resolution='h', area_thresh=10.) #lambert conformal -- should match parameters in arps.input!

#create a contour plot
CS = map.contour(x, y, prob_field, probability_contours, linewidths = 1.5, colors = prob_colors)
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
   map.readshapefile(shapefile_dir + '/county/countyp020','counties',drawbounds=True, linewidth=0.5, color='gray')  #Draws US county boundaries.

if draw_urban == 1:
   print 'reading urban boundaries from shapefile'
   if color_plots == 1:
      map.readshapefile(shapefile_dir + '/urban/geo_2000_d_dec_00000_400_py_z1', 'urban areas', drawbounds=True, linewidth=0.5, color='purple')  #Draws urban areas in purple.
   else:
      map.readshapefile(shapefile_dir + '/urban/geo_2000_d_dec_00000_400_py_z1', 'urban areas', drawbounds=True, linewidth=0.5, color='gray')  #Draws urban areas in gray.
   #end if
   
if draw_radar == 1:
   print 'reading radar range rings from shapefile'
   map.readshapefile(shapefile_dir + '/CASA/AllRadarRings', 'radars', drawbounds=True, linewidth=1.0, color='black')  #Draws radar range rings in black, bold lines

map.drawmapboundary()

meridians = arange(-110., 90., 2)     #Sets up meridians from 110W to 90W at 2 degree increments.
map.drawmeridians(meridians, labels=[0,0,0,1])  #Draws meridians (lines of constant longitude)

parallels = arange(30., 40., 2)       #Sets up parallels from 30N to 40N at 2 degree increments.
map.drawparallels(parallels, labels=[1,0,0,0])    #Draws parallels (lines of constant latitude)

#set plot title
plt.title('Ensemble-based probability of {%s > %s %s} -- 0330Z, 9 May 2007' %(var, var_threshold, var_units))

#Overlay the arbitrary contour field on the map.
#map.contour(self, X, Y, Z)

plt.show()

