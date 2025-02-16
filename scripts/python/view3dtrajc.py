#################################################################################
# FILE:		view3dtrajc.py
# CREATED:      2010-07-30 by Brett Roberts/Nick Engerer
# MODIFIED:	2011-07-01 by Brett Roberts
# REQUIRES:	Mayavi
#
# Renders a 3D interactive plot of trajectories from an arpstrajc output file.
#
# Usage:
#	python /abspath/view3dtrajc.py runname.trajc_xxx output_dir
#
#	NOTES:
#	- You must explicitly call python at the command line, rather than
#	   simply executing the script (e.g., ./view3dtrajc.py will not work)
#	- You must specify the absolute path of this script at the command
#	   line, even in the event that you are already in the directory
#	   which contains it
#################################################################################


import sys
import os
from numpy import *
from enthought.mayavi import mlab

filepath = os.path.dirname(sys.argv[0])
sys.path.append(filepath + '/../modules/')

import arpstrajc_lib as atl
import mayavi_plot_lib as mpl

n_args = len(sys.argv)

if(n_args < 1):
    print "ERROR: You must specify an arpstrajc output file as input."
    exit
else:
    infile = sys.argv[1]

# Image size
imgsize = (1000,800)
# Set of colors for painting trajectories (iterates through and repeats if necessary)
color_table = [(0,0,0),(0.5,0,0),(0,0.5,0),(0,0,0.5),(0.5,0.5,0),(0.5,0,0.5),(0,0.5,0.5)]

# Open trajectory file
try:
        # Trajectory data
        trajc_data = atl.get_trajc_array(infile)
except IOError:
        print "ERROR: The input file you specified either is not a valid arpstrajc file or does not exist."
        exit

# Trajectory domain
trajc_domain = atl.get_domain_range(trajc_data)

# Create Mayavi plot object
plot = mlab.figure(figure='plot3dtracj', size=imgsize)

# Render horizontal grid on xy-plane
mpl.plot_horiz_grid(plot, [trajc_domain[0],trajc_domain[1]])

#Render each trajectory
i = 0
while i < len(trajc_data):
    this_trajc = trajc_data[i]
    mpl.plot_3d_trajc(plot, this_trajc, str(i), trajc_color=color_table[i%len(color_table)])
    i += 1

mlab.orientation_axes(figure = plot)

mlab.show()
