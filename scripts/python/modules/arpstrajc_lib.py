# arpstrajc_lib.py
#
# Contains functions to read and manipulate data from
# arpstrajc output files.
#
#
# LIST OF FUNCTIONS:
#
# get_trajc_array( traj_file )
# get_runname( traj_file )
# get_timestep_count( traj_file )
# get_trajc_count( traj_file )
# get_domain_range( traj_array )
#

import sys


#############################################################
# FUNCTION:                                                 #
#   get_trajc_array( traj_file )                            #
#                                                           #
# Purpose:                                                  #
#   Reads in an arpstrajc output file and creates a system  #
#   of Python lists containing the trajectory coordinates.  #
#                                                           #
# Arguments:                                                #
#   traj_file   : filename (absolute path) of .trajc_xxx    #
#                   file                                    #
#                                                           #
# Output:                                                   #
#   Array containing a dictionary containing arrays for the #
#   x, y, and z coordinates of each trajectory. The         #
#   structure is as follows:                                #
#       traj_array["x":[coords],"y":[coords],"z":[coords]]  #
#############################################################
def get_trajc_array( traj_file ):
    """Reads in an arpstrajc output file and creates a system
    of Python lists containing the trajectory coordinates."""
    
    f = open(traj_file)

    #First line contains run name; skip
    run_name = f.readline()
    run_name = run_name.strip()

    #Skip over lines containing initial specifications
    tmp = f.readline()
    while(len(tmp.split()) > 1):
          tmp = f.readline()

    #Set ntimes to number of time steps
    ntimes = int(tmp)

    #Set ntrajcs to number of trajectories
    tmp = f.readline()
    ntrajcs = int(tmp)

    #Array of trajectories (to be returned)
    all_trajcs = []

    #Populate array with blank sub-dictionaries for each trajectory
    i = 1
    while i <= ntrajcs:
        all_trajcs.append({"x":[],"y":[],"z":[],"t":[]})
        i += 1

    #Loop through all time steps
    i = 1
    while i <= ntimes:
        #First line contains time
        time = f.readline()
        time = time.rstrip('\r\n')
        time = time.lstrip('     ')

        #Skip second line (number of trajectories)
        tmp = f.readline()

        lines_per_timestep = (ntrajcs / 2)
        if((ntrajcs % 2) != 0):
            lines_per_timestep += 1.

        #Loop through each trajectory at the current time step
        j = 1
        cnt = 1
        
        while j <= lines_per_timestep:
            tmp = f.readline()
            tlist = tmp.split()

            #Retrieve data for trajectory(ies) on current line
            #   note: there can be either 1 or 2 trajectories per line
            if len(tlist) == 3:
                all_trajcs[cnt-1]["x"].append(float(tlist[0]))
                all_trajcs[cnt-1]["y"].append(float(tlist[1]))
                all_trajcs[cnt-1]["z"].append(float(tlist[2]))
                all_trajcs[cnt-1]["t"].append(float(time))
                cnt += 1
            elif len(tlist) == 6:
                all_trajcs[cnt-1]["x"].append(float(tlist[0]))
                all_trajcs[cnt-1]["y"].append(float(tlist[1]))
                all_trajcs[cnt-1]["z"].append(float(tlist[2]))
                all_trajcs[cnt-1]["t"].append(float(time))
                cnt += 1
                all_trajcs[cnt-1]["x"].append(float(tlist[3]))
                all_trajcs[cnt-1]["y"].append(float(tlist[4]))
                all_trajcs[cnt-1]["z"].append(float(tlist[5]))
                all_trajcs[cnt-1]["t"].append(float(time))
                cnt += 1
            else:
                sys.exit()

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
#   Name of ARPS run specified in trajc file.               #
#############################################################
def get_run_name( traj_file ):
    """Retrieves the run name from a trajc file."""
    
    f = open(traj_file)

    #First line contains run name
    run_name = f.readline()
    run_name = run_name.strip()

    return run_name


#############################################################
# FUNCTION:                                                 #
#   get_timestep_count( traj_file )                         #
#                                                           #
# Purpose:                                                  #
#   Retrieves the number of time steps from a trajc file.   #
#                                                           #
# Arguments:                                                #
#   traj_file : filename (absolute path) of .trajc file     #                             #
#                                                           #
# Output:                                                   #
#   Number of time steps in the trajc file.                 #
#############################################################
def get_timestep_count( traj_file ):
    """Retrieves the number of time steps from a trajc file."""
    
    f = open(traj_file)

    #First line contains run name; skip
    run_name = f.readline()

    #Skip over lines containing initial specifications
    tmp = f.readline()
    while(len(tmp.split()) > 1):
          tmp = f.readline()

    #This line contains time step count
    ntimes = int(tmp)

    return ntimes


#############################################################
# FUNCTION:                                                 #
#   get_trajectory_count( traj_file )                       #
#                                                           #
# Purpose:                                                  #
#   Retrieves the number of trajectories from a trajc file. #                                                   #
#                                                           #
# Arguments:                                                #
#   traj_file : filename (absolute path) of .trajc file     #
#                                                           #
# Output:                                                   #
#   Number of trajectories in the trajc file.               #
#############################################################
def get_trajectory_count( traj_file ):
    """Retrieves the number of trajectories from a trajc file."""
    
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
    """Finds the x, y, and z coordinates which encompass the
    entire domain covered by a set of trajectories."""
    
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
