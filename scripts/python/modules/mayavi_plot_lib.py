#################################################################################
# FILE:         mayavi_plot_lib.py
# CREATED:      2010-07-30 by Brett Roberts
# MODIFIED:     2011-06-30 by Brett Roberts
# REQUIRES:	Mayavi
#
# Contains functions to aid in producing three-dimensional Mayavi plots,
# particularly related to trajectories.
#
#
# LIST OF FUNCTIONS:
#
# plot_horiz_grid( plot_obj, horiz_domain[, grid_color[, tubesize]])
# plot_3d_trajc( plot_obj, trajc_timeseries[, label[, trajc_color[, label_color[, tubesize]]]])
#
#################################################################################

from numpy import *
from enthought.mayavi import mlab

#############################################################
# DEFAULT VALUES                                            #
#############################################################

#Grid color
default_grid_color = (1, 1, 1)
#Radius of tubes comprising grid
default_grid_tubesize = 15

#Trajectory color
default_trajc_color = (1, 0, 0)
#Trajectory initial-point label color
default_trajc_label_color = (1, 1, 1)
#Radius of tubes comprising trajectories
default_trajc_tubesize = 50


#############################################################
# FUNCTION:                                                 #
#   plot_horiz_grid( plot_obj, horiz_domain[, grid_color[,  #
#       tubesize]])                                         #
#                                                           #
# Purpose:                                                  #
#   Plots a horizontal grid on the xy-plane                 #
#                                                           #
# Arguments:                                                #
#   plot_obj : mlab object to which grid is drawn           #
#   horiz_domain : 2D array giving horizontal domain, i.e.  #
#       [[xmin,xmax],[ymin,ymax]]                           #
#                                                           #
# Output:                                                   #
#   none                                                    #
#############################################################
def plot_horiz_grid( plot_obj, horiz_domain, z = 0, grid_color = default_grid_color, tubesize = default_grid_tubesize ):
	"""Plots a horizontal grid on the xy-plane"""

	domain_xrange = horiz_domain[0][1] - horiz_domain[0][0]
	domain_yrange = horiz_domain[1][1] - horiz_domain[1][0]
	domain_hrange = max([domain_xrange,domain_yrange])

	domain_scale = math.pow(10, floor(math.log10(domain_hrange)))
	grid_inc = ceil(domain_hrange / domain_scale) * (domain_scale / 10)

	grid_xmin = horiz_domain[0][0] - (horiz_domain[0][0] % grid_inc)
	grid_xmax = horiz_domain[0][1] + (grid_inc - (horiz_domain[0][1] % grid_inc))
	grid_ymin = horiz_domain[1][0] - (horiz_domain[1][0] % grid_inc)
	grid_ymax = horiz_domain[1][1] + (grid_inc - (horiz_domain[1][1] % grid_inc))

	i = grid_xmin
	while i <= grid_xmax:
		mlab.plot3d([i,i], [grid_ymin,grid_ymax], [z,z], figure=plot_obj, color=grid_color, tube_radius=tubesize)
		mlab.text3d(i, grid_ymin, 0, str(int(i)) ,figure=plot_obj, color=grid_color, scale=grid_inc/10)
		i += grid_inc

	i = grid_ymin
	while i <= grid_ymax:
		mlab.plot3d([grid_xmin,grid_xmax], [i,i], [z,z], figure=plot_obj, color=grid_color, tube_radius=tubesize)
		mlab.text3d(grid_xmin, i, 0, str(int(i)), figure=plot_obj, color=grid_color, scale=grid_inc/10)
		i += grid_inc

	return


#############################################################
# FUNCTION:                                                 #
#   plot_3d_trajc( plot_obj, trajc_timeseries[, label[,     #
#       trajc_color[, label_color[, tubesize]]]])           #
#                                                           #
# Purpose:                                                  #
#   Plots a 3D trajectory as a continuous line.             #
#                                                           #
# Arguments:                                                #
#   plot_obj : mlab object to which grid is drawn           #
#   trajc_timeseries : dictionary of arrays of the format   #
#       {'x':[x0,x1,...],'y':[y0,y1,...],'z':[z0,z1,...]}   #
#                                                           #
# Output:                                                   #
#   none                                                    #
#############################################################
def plot_3d_trajc( plot_obj, trajc_timeseries, label = '', trajc_color = default_trajc_color, label_color = default_trajc_label_color, tubesize = default_trajc_tubesize):
	"""Plots a 3D trajectory as a continuous line."""
    
	mlab.plot3d(trajc_timeseries["x"], trajc_timeseries["y"], trajc_timeseries["z"], figure=plot_obj, tube_radius=tubesize, color=trajc_color)
	mlab.text3d(trajc_timeseries["x"][0], trajc_timeseries["y"][0], trajc_timeseries["z"][0], label, figure=plot_obj, color=label_color, scale=tubesize*2)
