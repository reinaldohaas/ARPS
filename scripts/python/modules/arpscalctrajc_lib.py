#################################################################################
# FILE:         arpscalctrajc_lib.py
# CREATED:      2011-06-29 by Brett Roberts
# MODIFIED:     2011-08-23 by Brett Roberts
#
# Contains functions to read and manipulate data from
# arpscalctrajc output files.
#
#
# LIST OF CLASSES:
#
# calctrajcdata
#
# LIST OF FUNCTIONS:
#
# get_calctrajc_array( traj_file )
# get_header( traj_file )
# get_run_name( traj_file )
# get_ntrajcs( traj_file )
# get_npnts( traj_file )
# get_domain_range( traj_file ) 
#

import sys
import math


#############################################################
# FUNCTION:                                                 #
#   get_calctrajc_array( traj_file )                        #
#                                                           #
# Purpose:                                                  #
#   Reads in an arpscalctrajc output file and creates a     #
#   system of Python lists containing the data.             #
#                                                           #
# Arguments:                                                #
#   traj_file   : filename (absolute path) of the           #
#                 .trajcxxx.data file                       #
#                                                           #
# Output:                                                   #
#   Trajectory array containing an array of data arrays for #
#   each timestep. The structure is as follows:             #
#      [traj1[[t1,x1,y1,...],[t2,x2,y2,...],...],...]       #
#############################################################
def get_trajc_array( traj_file ):
    """Reads in an arpscalctrajc output file and creates a system of Python lists containing the data."""
    
    f = open(traj_file)

    # First line contains run name; skip
    run_name = f.readline()
    run_name = run_name.strip()

    # Skip over lines containing initial specifications
    tmp = f.readline()
    while(len(tmp.split()) > 1):
          tmp = f.readline()

    # Set ntimes to number of time steps
    ntimes = int(tmp)

    # Set ntrajcs to number of trajectories
    tmp = f.readline()
    ntrajcs = int(tmp)

    # Retrieve header array
    header = get_header(traj_file)
    f.readline()

    # Array of trajectories (to be returned)
    all_trajcs = []

    # Populate array with blank sub-arrays for each trajectory
    i = 1
    while i <= ntrajcs:
        all_trajcs.append([])
        i += 1

    # Loop through all time steps
    i = 1
    while i <= ntimes:
        # First line contains time
        time = f.readline()
        time = time.rstrip('\r\n')
        time = time.lstrip('     ')

        # Loop through each trajectory at the current time step
        j = 1
        
        while j <= ntrajcs:
            tmp = f.readline()
            tlist = tmp.split()

            all_trajcs[j-1].append(tlist)

            j += 1

        i += 1

    return all_trajcs


#############################################################
# FUNCTION:                                                 #
#   get_run_name( traj_file )                               #
#                                                           #
# Purpose:                                                  #
#   Retrieves the run name from a trajc file.               #
#                                                           #
# Arguments:                                                #
#   traj_file   : filename (absolute path) of .trajc_xxx    #
#                   file                                    #
#                                                           #
# Output:                                                   #
#   Name of ARPS run specified in calctrajc file.           #
#############################################################
def get_run_name( traj_file ):
    """Retrieves the run name from a calctrajc file."""
    
    f = open(traj_file)

    #First line contains run name
    run_name = f.readline()
    run_name = run_name.strip()

    return run_name


#############################################################
# FUNCTION:                                                 #
#   get_npnts( traj_file )                                  #
#                                                           #
# Purpose:                                                  #
#   Retrieves the number of time steps from a trajc file.   #
#                                                           #
# Arguments:                                                #
#   traj_file : filename (absolute path) of .trajc file     #
#                                                           #
# Output:                                                   #
#   Number of time steps in the trajc file.                 #
#############################################################
def get_npnts( traj_file ):
    """Retrieves the number of time steps from a calctrajc file."""
    
    f = open(traj_file)

    #First line contains run name; skip
    run_name = f.readline()

    #Skip over lines containing initial specifications
    tmp = f.readline()
    while(len(tmp.split()) > 1):
          tmp = f.readline()

    #This line contains time step count
    npnts = int(tmp)

    return npnts


#############################################################
# FUNCTION:                                                 #
#   get_ntrajcs( traj_file )                                #
#                                                           #
# Purpose:                                                  #
#   Retrieves the number of trajectories from a trajc file. #
#                                                           #
# Arguments:                                                #
#   traj_file : filename (absolute path) of .trajc file     #
#                                                           #
# Output:                                                   #
#   Number of trajectories in the trajc file.               #
#############################################################
def get_ntrajcs( traj_file ):
    """Retrieves the number of trajectories from a calctrajc file."""
    
    f = open(traj_file)

    #First line contains run name; skip
    run_name = f.readline()

    #Skip over lines containing initial specifications
    tmp = f.readline()
    while(len(tmp.split()) > 1):
          tmp = f.readline()

    #This line contains number of trajectories
    tmp = f.readline()
    ntrajcs = int(tmp)
    
    return ntrajcs


#############################################################
# FUNCTION:                                                 #
#   get_header( traj_file )                                 #
#                                                           #
# Purpose:                                                  #
#   From a calctrajc output file, retrieves the header,     #
#   which specifies the variables written out and their     #
#   order.                                                  #
#                                                           #
# Arguments:                                                #
#   traj_file : filename (absolute path) of .trajc file     #
#                                                           #
# Output:                                                   #
#   Header from the calctrajc file.                         #
#############################################################
def get_header( traj_file ):
    """Retrieves the header from a calctrajc file."""
    
    f = open(traj_file)

    # First line contains run name; skip
    run_name = f.readline()

    # Skip over lines containing initial specifications
    tmp = f.readline()
    while(len(tmp.split()) > 1):
          tmp = f.readline()

    # Skip over line containing number of trajectories
    tmp = f.readline()

    # This line contains the header string
    header_txt = f.readline()

    # Remove all whitespace from the header string
    header_txt = header_txt.replace(' ','')
    header_txt = header_txt.replace('\t','')
    header_txt = header_txt.strip(' \r\n')

    # Split the header string at commas
    header = header_txt.split(',')
    
    return header


#############################################################
# FUNCTION:                                                 #
#   get_domain_range( traj_array )                          #
#                                                           #
# Purpose:                                                  #
#   Finds the x, y, and z coordinates which encompass the   #
#   entire domain covered by a set of trajectories.         #
#                                                           #
# Arguments:                                                #
#   traj_array   : array returned from get_trajc_array()    #
#                                                           #
# Output:                                                   #
#   A two-dimensional array of the following format:        #
#       [[xmin,xmax],[ymin,ymax],[zmin,zmax]]               #
#############################################################
def get_domain_range( traj_array ):
    """Finds the x, y, and z coordinates which encompass the entire domain covered by a set of trajectories."""
    
    xmin = min(traj_array[0]['x'])
    xmax = max(traj_array[0]['x'])
    ymin = min(traj_array[0]['y'])
    ymax = max(traj_array[0]['y'])
    zmin = min(traj_array[0]['z'])
    zmax = max(traj_array[0]['z'])

    for traj,arr in enumerate(traj_array):
        this_xmin = min(arr['x'])
        this_xmax = max(arr['x'])
        this_ymin = min(arr['y'])
        this_ymax = max(arr['y'])
        this_zmin = min(arr['z'])
        this_zmax = max(arr['z'])

        if(this_xmin < xmin):
            xmin = this_xmin
        if(this_xmax > xmax):
            xmax = this_xmax
        if(this_ymin < ymin):
            ymin = this_ymin
        if(this_ymax > ymax):
            ymax = this_ymax
        if(this_zmin < zmin):
            zmin = this_zmin
        if(this_zmax > zmax):
            zmax = this_zmax
        
    return [[xmin,xmax],[ymin,ymax],[zmin,zmax]]


#################################################################
# CLASS calctrajcdata						#
# 								#
# An interface for the data contained within an arpscalctrajc	#
# file. In addition to providing accessor methods to basic	#
# information such as the run name, number of trajectories, and	#
# number of points, it also provides methods to access data in	#
# more useful formats than are possible by reading sequentially	#
# from the .data file.						#
#								#
#								#
# ACCESSOR METHODS						#
#								#
# get_ntrajcs() : returns the number of trajectories		#
#								#
# get_npnts() : returns the number of points per trajectory	#
#								#
# get_dt() : returns the timestep size				#
#								#
# get_header() : returns an array containing the calculated	#
#	fields, in order					#
#								#
# get_run_name() : returns the run name				#
#								#
# get_timeseries(trajc_id, field) : returns an array containing	#
#	the timeseries of the specified field for the		#
#	specified trajectory					#
#								#
# get_averaged_timeseries(field) : returns an array		#
#	containing the timeseries of the specified field for	#
#	the mean of all trajectories in the file		#
#								#
# get_summed_timeseries(trajc_id, fields) : returns an array	#
#	containing a timeseries for the specified trajectory	#
#	constructed by taking the sum of all specified fields	#
#	at each timestep. For example, if each RHS term of a	#
#	time-tendency equation is diagnosed as a field in the	#
#	.data file, you might want to sum all of them to	#
#	determine the time tendency at each timestep.		#
#								#
# get_difference_timeseries(trajc_id, field1, field2) : returns	#
#	an array containing a timeseries for the specified	#
#	trajectory constructed by taking the difference between	#
#       the two fields at each timestep. May be useful for	#
#	plotting errors, for example.				#
#								#
# get_derivative_timeseries( timeseries ) : returns an array	#
#       containing a timeseries that represents the time	#
#       derivative of the input array.				#
#								#
# get_integrated_timeseries( timeseries, init_value ) : returns	#
#	an array containing a timeseries that represents an	#
#	integration of the input array starting from the	#
#	specified initial value.				#
#################################################################
class calctrajcdata:
	"""An interface for the data contained within an arpscalctrajc file."""

	def __init__( self, _filename ):
		self.filename = _filename

		self.ntrajcs = get_ntrajcs(self.filename)
		self.npnts = get_npnts(self.filename)
		self.header = get_header(self.filename)
		self.run_name = get_run_name(self.filename)
		self.trajc_array = get_trajc_array(self.filename)

	def get_ntrajcs( self ):
	        """Returns the number of trajectories in the dataset."""
		return self.ntrajcs

	def get_npnts( self ):
		"""Returns the number of points (timesteps) per trajectory in the dataset."""
		return self.npnts

	def get_header( self ):
		"""Returns a list of the variables in the dataset."""
		return self.header

	def get_run_name( self ):
		"""Return the name of this run."""
		return self.run_name

	def get_timeseries( self, trajc_id, field ):
	       	"""Returns a timeseries of the specified field for the specified trajectory."""
		if( (trajc_id < 0) or (trajc_id >= self.ntrajcs) ):
			return "ERROR"

		field_index = self.get_field_index(field)

		if(field_index >= 0):
			ts = []
			i = 0
			while(i < self.npnts):
				ts.append(float(self.trajc_array[trajc_id][i][field_index]))
				i+=1

			return ts
		else:
			return "ERROR"

        def get_averaged_timeseries( self, field ):
                """Returns a timeseries of the mean of all trajectories for specified field."""
                field_index = self.get_field_index(field)

                if(field_index >= 0):
                        ts = []
                        i = 0
                        while(i < self.npnts):
				j = 0
				isum = 0.0
				while(j < self.ntrajcs):
					isum = isum + float(self.trajc_array[j][i][field_index])

					j+=1

				imean = isum / self.ntrajcs
				ts.append(imean)

                                i+=1

                        return ts
                else:
                        return "ERROR"

	def get_summed_timeseries( self, timeseries_set ):
		"""Returns a timeseries representing the sum of multiple constituent series."""

		ts = []

		i = 0
		while(i < self.npnts):
			value = 0.0
			
			for cts in timeseries_set:
				value = value + float(cts[i])
                
			ts.append(value)
			i+=1

		return ts

	def get_difference_timeseries( self, field1, field2 ):
		"""Returns a timeseries of (field1 - field2)."""

		ts = []

		i = 0
		while(i < self.npnts):
			value = float(field1[i]) - float(field2[i])
			ts.append(value)
			i+=1

		return ts

	def get_derivative_timeseries( self, timeseries ):
		"""Returns a timeseries of the time derivative of the provided series."""

		dt = self.get_dt()
		dts = []

		# Note: use forward difference for first point, backward difference for last,
		#  and centered for every point in between.
		dt0 = (timeseries[1] - timeseries[0]) / dt
		dtf = (timeseries[-1] - timeseries[-2]) / dt

		dts.append(dt0)

		i = 1
		while(i < (len(timeseries) - 1)):
			dtcur = (timeseries[i+1] - timeseries[i-1]) / (2*dt)
			dts.append(dtcur)
			i+=1

		dts.append(dtf)

		return dts

	def get_integrated_timeseries( self, timeseries, init_value ):
		"""Returns the time integration of the provided timeseries, starting from a user-defined initial value."""

		curval = float(init_value)
		its = [curval]

		i = 1
		while(i < (len(timeseries))):
			avgdif = (timeseries[i-1] + timeseries[i]) / 2.0
			curval = curval + (self.get_dt() * avgdif)
			its.append(curval)
			i+=1

		return its

	def get_vectorsum_timeseries( self, xfield, yfield ):
                """Returns a timeseries of the vector magnitude of a set of x,y components."""

                ts = []

                i = 0
                while(i < self.npnts):
                        value = sqrt((float(xfield[i]))**2 + (float(yfield[i]))**2)
                        ts.append(value)
                        i+=1

                return ts

	def get_relerror_timeseries( self, basefield, compfield ):
		"""Returns a timeseries of relative error, defined as abs((base - comp)/base)."""

		ts = []

		i = 0
                while(i < self.npnts):
                        value = abs((float(basefield[i]) - float(compfield[i]))/float(basefield[i]))
                        ts.append(value)
                        i+=1

		return ts

	def average_value( self, timeseries ):
		"""Returns the average value over a single timeseries."""

		sum = 0.0

		for tsval in timeseries:
			sum+=float(tsval)

		return float(sum/len(timeseries))

	def average_abs_value( self, timeseries ):
                """Returns the average absolute value over a single timeseries."""

                sum = 0.0

                for tsval in timeseries:
                        sum+=float(abs(tsval))

                return float(sum/len(timeseries))

	def get_field_index( self, field ):
		"""Returns the array index for the specified field."""
		field_index = -1

		i = 0
		while(i < len(self.header)):
			if(self.header[i] == field):
				field_index = i
				return i
			i = i + 1

		return -999

	def get_dt( self ):
		"""Returns the timestep size for the dataset."""

		ti = self.get_field_index('t')
		return float(self.trajc_array[0][1][ti]) - float(self.trajc_array[0][0][ti])
