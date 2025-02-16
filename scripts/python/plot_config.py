
#------------------------------------------------------------------------#
#                                                                        #
#           ARPS PYTHON PLOTTING SCRIPTS CONFIGURATION BLOCK             #
#                                                                        #
#------------------------------------------------------------------------#
#    Author:  Nate Snook (Apr. 2009)                                     #
#             nsnook@ou.edu                                              #
#------------------------------------------------------------------------#
#    Modification History:                                               #
#       --> Record all modifications you make to the code HERE.          #
#           It will make your life easier -- trust me.                   #
#------------------------------------------------------------------------#

# This is the CONFIGURATION FILE for the ARPS Python Plotting Scripts.

# Please edit this file BEFORE you attempt to run any of the plotting
# scripts.  Do not add new variables or delete current ones without also 
# making the corresponding changes to the individual plotting utilities; 
# failing to do so will likely break them.

# You can get updated versions of these utilities from Nate Snook (nsnook@ou.edu)

#------------------------------------------------------------------------#
#                       Standard Python Libraries                        #
#------------------------------------------------------------------------#

#DO NOT change the standard library declarations unless you know what you're doing!
from numpy import arange   #needed to define the contour intervals

#------------------------------------------------------------------------#
#                           Plotting Resources                           #
#------------------------------------------------------------------------#

#Locations of needed shapefiles and other resources are declared here.
#DO NOT change these declarations unless you know what you're doing!

#Where are the shapefiles for use with basemap located?  Please keep them all in one place.
shapefile_dir = '/data5/nsnook/python_scripts/shapefiles/'


#------------------------------------------------------------------------#
#                          General Input Files                           #
#------------------------------------------------------------------------#

# NOTE:  When setting file names, use the FULL PATH for best results.  
#        If you use a relative path, it may cause problems.

# Set your grid base filename.  Lat, lon, and domain data will be read from here.
grdbas = '/data4/nsnook/9may2007_2km/CNTL/enf001.hdfgrdbas'
	    
# For plots of a single forecast or ensemble member, list the history dump file to use data from:
hdfhistory = '/data4/nsnook/9may2007_2km/CNTL/enf007.hdf025200'

# If you're plotting a difference between two forecasts or ensemble members (field_difference)
# then list the member to compare to here -- when you run field_difference
# it will show values of (hdfhistory - hdfhistory2)
hdfhistory2 = '/data4/nsnook/9may2007_2km/CNTL/enf003.hdf025200'

#------------------------------------------------------------------------#
#                         Microphysics Options                           #
#------------------------------------------------------------------------#
	    
# If you're working with an ensemble, then the number 
# of members from various microphysics types are...

#For 'CNTL' and '88D'/'NoCASA' experiments
#num_lin = 16     #Number of 88D (in ARPS: mphyopt = 2) members
#num_wsm = 16     #Number of WSM6 (in ARPS: mphyopt = ?) members
#num_schultz = 8  #Number of Schultz (in ARPS: mphyopt = ?) members

#For 'Lin'/'NoMMP' experiment
num_lin = 40
num_wsm = 0
num_schultz = 0
# NOTE:  If not using an ensemble, these values will be ignored.


#------------------------------------------------------------------------#
#               Probabilistic/Ensemble Input and Options                 #
#------------------------------------------------------------------------#

#How many total ensemble members are you using?
n_ens = num_lin + num_wsm + num_schultz

#What forecast time are you using?  Feel free to add a local time to the string, if desired
utc_time = '0100 UTC' #Set this string to the UTC time associated with your file.
hdftime = '025200'  #This is the corresponding 'tttttt' time used in the hdf file suffix.
arpstime =  str(hdftime) + 's (' + str(utc_time) + ')'  #This string will appear in your plot titles.

#What is your experiment name?
exp_name = '9may2007_2km_CNTL'

# For probabilistic forecasts and ensemble mean plots, you will need to
# specify which files to use for EACH member of the ensemble 
# You should EXPLICITLY LIST ALL history files you want to use HERE.
# If including an ensemble mean, list that file LAST.

# For all variables directlystored in history dumps, give those files (e.g. 'arpsrun.hdf007200')
# For radar reflectivity, give truref files generated by ossedata (e.g. 'arpsrun.hdftruref007200')

# In all these lists, you can give the complete filename (e.g. '<path_to_data>/arpsrun.hdf003600' )
# or make use of variables from this file (e.g. '<path_to_data>/arpsrun.hdf' + str(hdftime) )

filelist = ['/data4/nsnook/9may2007_6km_relax/enf001.hdf' + str(hdftime),
            '/data4/nsnook/9may2007_6km_relax/enf002.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf003.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf004.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf005.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf006.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf007.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf008.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf009.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf010.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf011.hdf' + str(hdftime),
            '/data4/nsnook/9may2007_6km_relax/enf012.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf013.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf014.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf015.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf016.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf017.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf018.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf019.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf020.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf021.hdf' + str(hdftime),
            '/data4/nsnook/9may2007_6km_relax/enf022.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf023.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf024.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf025.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf026.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf027.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf028.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf029.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf030.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf031.hdf' + str(hdftime),
            '/data4/nsnook/9may2007_6km_relax/enf032.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf033.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf034.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf035.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf036.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf037.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf038.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf039.hdf' + str(hdftime),
	    '/data4/nsnook/9may2007_6km_relax/enf040.hdf' + str(hdftime),
            ]
#If you have an ensemble mean, put it here:
#            '/data4/nsnook/9may2007_6km_relax/enmean.hdf' + str(hdftime)]

#If you have an observation file(s), put it here
obsfiles = ['/data5/nsnook/radar/KTLX/enmean.hdfrefl3d' + str(hdftime),
            '/data5/nsnook/radar/KDYX/enmean.hdfrefl3d' + str(hdftime),
            '/data5/nsnook/radar/KFWS/enmean.hdfrefl3d' + str(hdftime),
            '/data5/nsnook/radar/KAMA/enmean.hdfrefl3d' + str(hdftime),
            '/data5/nsnook/radar/KLBB/enmean.hdfrefl3d' + str(hdftime)]
           
#If you have model gridded radar data generated by osse_data, list the file(s) here.
#osse_files = ['/data5/nsnook/kraken_data/may20.hdftruref018600']
osse_files = ['/data5/nsnook/2km_CNTL025/CNTL_NEM_deterministic.hdftruref' + str(hdftime)]

# Set probabilistic parameters
neighborhood_radius = 5000.  #radius of neighborhood to use for 'fuzzy' verification (in METERS)
verif_radius = 25000.  #Verification radius for 'probability of <feature> within <distance> of a point', (in METERS)
prob_plot_obs = True   #Do you want to plot a verification (observation) on the plot?  'True' or 'False'.

# Set thresholds for 'probability of <feature> within <distance> of a point'.
vort_verif_thres = 0.000

# Set the threshold for probabilstic forecasting for point variables
# (i.e. the value at a point must exceed this value to verify as 'present')
#                          #VARIABLE             -- UNITS
#---------------------------------------------------------------------
qh_threshold = 0.000025    #hail mixing ratio    -- unitless (multiply by 1000 to get g/kg)
qr_threshold = 0.0030      #rain mixing ratio    -- unitless (multiply by 1000 to get g/kg)
dBZ_threshold = 40.0       #radar reflectivity   -- dBZ
vort_threshold = 0.00025   #vertical vorticity   -- s^-1
prcrate_threshold = 0.0017639 #precip rate       -- kg/(m^2*s)  (default of 0.0017639 = 0.25 in/hr (approx.))
w_threshold = 10.0         #vertical velocity    -- m/s

#------------------------------------------------------------------------#
#                           Plotting Options                             #
#------------------------------------------------------------------------#

# -----Basic options----- #

casaon = 0       # Do you want to include the CASA radars in this figure?  (1 = yes, 0 = no)
colorcode = 0    # Do you want to color-code the figure background for CASA vs. WSR-88D radars? ( 1 = yes, 0 = no)
color_plots = 1       # Do you want color or grayscale plots?  (1 = color, anything else = grayscale)
add_title = 1         # Do you want to include a title on your plot?  (1 = yes, 0 = no)
obs_on = True        # For a radar observation plot, do you have observed radar data on the model grid?
osse_on = False        # For a radar observation plot, do you have model radar data generated using 'ossedata'?
dBZ_composite = False # For radar observations, do you want composite reflectivity?  'True' or 'False'
transparency = 0.8    # Between 0 and 1, what level of transparency do you want?  (1 = opaque, 0 = fully transparent)
y_min = 0.00     # For line, scatter, or bar graphs, what should be the y-axis minimum?
y_max = 16.00     # For line, scatter, or bar graphs, what should be the y-axis maximum?
y_step = 2.00    # For line, scatter, or bar graphs, what should be the y-axis step?

# ------------------------ #

# -----Domain options----- #
use_subdomain = False   #Would you like to use a subdomain for skill score calculations?

# WARNING -- THE CUSTOM DOMAIN OPTIONS HAVE NOT YET BEEN IMPLEMENTED!  DO NOT TURN THEM ON!
custom_domain = False  #Do you want to plot something other than the entire domain?  "True" or "False".
x_begin = -1  #x-coordinate (i) of the left (western) edge of the custom domain area to be plotted
y_begin = -1  #y-coordinate (j) of the bottom (southern) edge of the custom domain area to be plotted
x_end = -1    #x-coordinate (i) of the right (eastern) edge of the custom domain area to be plotted
y_end = -1    #y-coordinate (j) of the top (northern) edge of the custom domain area to be plotted


# ------------------------ #

# For a plot of a single field, enter the variable you want to plot here.
# It also sets which variable to use for 'vardump'
# Options include 'vort', 'qr', 'qh', 'ptprt', etc... (for a full list of choices, run 'varlist')
# 'Vr' and 'Z' are also options for radar-related plots!
var = 'v'

# If you want a single x-y slice, what vertical level (i.e. k = ?) do you want to plot at?
# Keep in mind that Python counts starting from ZERO, so '0' is the lowest model level.
lev = 1   #9

#DO NOT turn on both windbarbs and windvects -- you cannot have both!
# Do you want to plot wind barbs as well as your variable field?  (1 = yes, anything else = no)
windbarbs = 0

# Do you want to plot wind vectors as well as your variable field?  (1 = yes, anything else = no)
windvects = 0

# If windbarbs = 1 or windvects = 1, how often should we plot wind barbs or vectors?  
# (i.e. '1' = every gridpoint, '5' = every 5th gridpoint, etc...)
skip = 8

#What number of spatial dimensions does the variable listed above have? (should be 1, 2, or 3)
dimensions = 3

timeheight_top = 12000  #maximum height to show in time-height figures (in meters)

#The following options determine what maps are overlayed on your plot.
draw_radar = 1  #Do you want radar range rings for CASA plotted? (1 = yes, anything else = no)
draw_counties = 1  #Do you want county borders plotted?  (1 = yes, anything else = no)
draw_urban = 1  #Do you want urban boundaries plotted?  (1 = yes, anything else = no)

# Set the contour increments for various plotting options
# The convention is:  arange(<minimum contour value>, <maximum contour value>, <contour increment>)
probability_contours = [0.0, 0.05, 0.20, 0.40, 0.60, 0.80, 0.95, 1.00] #contours for probability (default from 0.0 to 1.0 by 0.1)
refl_contours = arange(10., 75., 5.)  #contours for reflectivity (from 10 dBZ to 70 dBZ by 5 dBZ)
refl_contours_tbl = arange(10.,90.,10.)       #Contours for reflectivity to mimic "refl_grayscale.tbl" from arpsplt
vort_contours = arange(0.000, 0.013, 0.002)   #contours for vertical vorticity (from 0.001 to 0.012 by 0.001)
pt_contours = arange(290.0, 301.0, 1.0)       #Contours for potential temperature (from 290 K to 300 K by 1)
w_contours = arange(0.0, 36.0, 4.0) #Contours for vertical velocity (from 0 m/s to 32 m/s by 4)
ptprt_contours = arange(-7.0, 4.0, 1.0)       #Contours for potential temp. perturbation (From -7 K to 4 K by 1)

#For spaghetti plots, we need to define which contour (or contours) to use.
spag_contour_refl = arange(40., 41., 1.)         #contours *only* the 40 dBZ contour line
spag_contour_vort = arange(0.004, 0.005, 0.001)  #contours *only* the 0.004 s^-1 contour line
spag_contour_qr = arange(0.003, 0.004, 0.001) #contours *only* the 3.0 g/kg contour line
spag_contour_qh = arange(0.000025, 0.000030, 0.000005)  #contours *only* the 0.03 g/kg contour line
spag_contour_prcrate = arange(0.0017639, 1.0017639, 1.)  #contours *only* the 0.25 in/hr contour line
spag_en_mean = True  #Do we want to show an ensemble mean line on the spaghetti plot? "True" or "False".
spag_obs = True      #Do we want to show an observed data contour on the spaghetti plot?  "True" or "False".

#If plotting a single field, should we filter bad data?  (1 = yes, anything else = no)
filt = 0

#What minimum and maximum values should be used to filter? 
#   NOTE: (values outside this range will be thrown out if filt = 1)
filt_min = -998.   #units are the same as those of 'var'
filt_max = 0.099   #units are the same as those of 'var'
filtered = -1.     #value to flag bad data with (data that fall outside the filter threshold)

#---------------------------------------------------------------------------------#
#                               Colormap options                                  #
#---------------------------------------------------------------------------------#

#The colormap is defined using official HTML colornames, which matplotlib supports.
#You can change the colormap to use any set of colors you want, but use official 
#accepted HTML color names or you may find that your code breaks.

#Approximates the usual <<blue --> green --> yellow --> red --> pink --> white>> radar colormap
#Designed for use from 5dBZ to 70dBZ
refl_colors = ['#00FFFF', '#6495ED', '#000090', '#00FF00', '#00BB00', '#008800', '#FFFF00',
               '#FFDD00', '#DAA520', '#FF0000', '#B21111', '#990000', '#FF00FF', '#BB55DD']

#A colormap of 16 shades of gray, with the first 4 shades extremely light
grays = ['#ffffff', '#eeeeee', '#dddddd', '#cccccc', '#bbbbbb', '#aaaaaa', '#999999', '#888888',
         '#777777', '#666666', '#555555', '#444444', '#181818', '#111111', '#0a0a0a', '#050505']

#A colormap of 16 shades of gray, from dark to light
gray_darktolight = ['#000000', '#222222', '#444444', '#666666', '#777777', '#888888',
         '#999999', '#aaaaaa', '#bbbbbb', '#cccccc', '#dddddd', '#eaeaea', '#ffffff']
         
#A colormap designed for probability plots:  green denotes low probabilities, yellow denotes
#moderate probabilities, and red denotes high probabilities.
prob_colors = ['#00FF88', '#00FF00', '#33FF00', '#CCFF00', '#FFFF00', '#FFCC00', '#FF3300',
               '#FF0000', '#CC0000']

#A colormap for probability plots in grayscale, for use in AMS Journal publications.
#prob_grays =  ['#CCCCCC', '#BBBBBB', '#999999', '#888888', '#666666', '#555555', '#333333',
#               '#222222', '#111111']                                  #For intervals of 0.1
prob_grays = ['#FFFFFF', '#D2D2D2', '#BDBDBD', '#AAAAAA', '#949494', '#737373', '#525252']   #For intervals of 0.2

#A colormap for probability plots in grayscale with lighter colors, for use in AMS Journal publications.
prob_light_grays =  ['#D9D9D9', '#CCCCCC', '#C0C0C0', '#B2B2B2', '#A5A5A5', '#929292', '#777777',
               '#5C5C5C', '#444444']

#A colormap to mimic "refl_grayscale.tbl" in arpsplt
refl_grayscale_tbl = [(.9,.9,.9), (.8,.8,.8), (.7,.7,.7), (.6,.6,.6), (.5,.5,.5),
                      (.4,.4,.4), (.25,.25,.25), (.1,.1,.1), (0,0,0), (0,0,0), (0,0,0)]


