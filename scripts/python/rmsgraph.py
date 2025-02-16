#!/usr/bin/env python2.7

#Code: rmsgraph.py (Python Script)
#Author: Nate Snook (CASA/CAPS/SoM)
#Written: Jun. 2008
#Last modified:  Jun. 30, 2008
#
#Purpose:
#   'rmsgraph.py' creates a multi-panel plot of the RMS statistics for a each radar in
#   an ARPS EnKF data assimilation experiment.  Plotted are the RMS innovation of Vr and
#   Z, as well as the RMS spread of those same two variables.

import os
import subprocess as sp
from numpy import *
import matplotlib.pyplot as plt
#from plot_config import *

#Initialize rms arrays
anavr = []  #RMS innovation (or error) of analysis for radial velocity (Vr)
anaz = []   #RMS innovation (or error) of analysis for reflectivity (Z)
bgvr = []   #RMS innovation (or error) of background (forecast) for radial velocity (Vr)
bgz =[]     #RMS innovation (or error) of background (forecast) for reflectivity (Z)
sprdanavr = [] #RMS spread of Vr
sprdanaz = []  #RMS spread of Z
sprdbgvr = []  #RMS spread of Vr in background
sprdbgz = []   #RMS spread of Z in background
allvr = []  #Combined anavr and bgvr for plotting purposes
allz = []   #Combined anaz and bgz for plotting purposes
allsprdvr = [] #Combined sprdanavr and sprdbgvr for plotting purposes
allsprdz = []  #Combined sprdanaz and sprdbgz for plotting purposes
times = []  #A listing of times for use in plotting

#
mult_inflat = 30  #what percentage was used for covariance inflation?

#List the home directory
homedir = '/home/nsnook/python_scripts/rms/CNTL0' + str(mult_inflat) + '/'

#List all the radars to analyze RMS stats for
radarlist = ['KAMA', 'KDYX', 'KLBB', 'KTLX', 'KVNX', 'KCYR', 'KLWE', 'KRSP', 'KSAO']

#List all the times that you wish to analyze stats for
timelist = ['003900', '004200', '004500', '004800', '005100', '005400', '005700', 
            '006000', '006300', '006600', '006900', '007200']

#initialize radar counter variable
nradar = 0

for radar in radarlist:
   curanavr = []     #values of the RMS variables (Vr and Z) for the...
   curanaz = []      #...analysis and background for the current radar.
   
   curbgvr = []      #Same thing for the background RMS innovation of Vr
   curbgz = []       #Background RMS innovation of Z
   
   curallvr = []     #Background and analysis innovation of Vr
   curallz = []      #Background and analysis innovation of Z
   
   curallsprdvr = [] #Background and analysis spread of Vr
   curallsprdz = []  #Background and analysis spread of Z
   
   cursprdbgvr = []  #RMS spread of Vr for the background
   cursprdbgz = []   #RMS spread of Z for the background
   
   cursprdanavr = [] #RMS spread of Vr for the analysis
   cursprdanaz = []  #RMS spread of Z for the analysis
   
   for time in timelist:
      anafile = homedir + radar + 'rmsdan' + time  #define analysis and...
      bgfile = homedir + radar + 'rmsdbg' + time   #...background filenames for RMS innovation...
      sprdanafile = homedir + radar + 'spreadana' + time  #and also define analysis and...
      sprdbgfile = homedir + radar + 'innov' + time  #...background filenames for RMS spread.
      
      thefile = open(bgfile, 'r')
      vr_in = float(thefile.read(15))   #Read in RMS innovation of Vr...
      z_in = float(thefile.read(15))    #...and Z for this radar for the background.
      curbgvr.append(vr_in)
      curbgz.append(z_in)
      curallvr.append(vr_in)
      curallz.append(z_in)
      if nradar == 0:    #We only need to collect the list of times once, thanks.
         times.append(float(time))
      #end if
      
      thefile = open(anafile, 'r')
      vr_in = float(thefile.read(15))   #Read in RMS innovation of Vr...
      z_in = float(thefile.read(15))    #...and Z for this radar for the analysis.
      curanavr.append(vr_in)
      curanaz.append(z_in)
      curallvr.append(vr_in)
      curallz.append(z_in)
      if nradar == 0:    #We only need to collect the list of times once, thanks.
         times.append(float(time))
      # end if
      
      thefile = open(sprdbgfile, 'r')
      vr_in = float(thefile.read(15))   #Read in RMS spread of Vr...
      z_in = float(thefile.read(15))    #...and Z for this radar for the analysis.
      cursprdbgvr.append(vr_in)
      cursprdbgz.append(z_in)
      curallsprdvr.append(vr_in)
      curallsprdz.append(z_in)
   
      thefile = open(sprdanafile, 'r')
      vr_in = float(thefile.read(15))   #Read in RMS spread of Vr...
      z_in = float(thefile.read(15))    #...and Z for this radar for the analysis.
      cursprdanavr.append(vr_in)
      cursprdanaz.append(z_in)
      curallsprdvr.append(vr_in)
      curallsprdz.append(z_in)
   #end for (time)
   
   anavr.append([])  #append a new array element for the current radar for all 12 variables
   anaz.append([])
   bgvr.append([])
   bgz.append([])
   sprdanavr.append([])
   sprdanaz.append([])
   sprdbgvr.append([])
   sprdbgz.append([])
   allvr.append([])
   allz.append([])
   allsprdvr.append([])
   allsprdz.append([])
   
   for item in curanavr:
      anavr[nradar].append(item)   #Store Vr data into the anavr array
   #end for
   for item in curanaz:
      anaz[nradar].append(item)    #...same for Z data
   #end for
   for item in curbgvr:
      bgvr[nradar].append(item)    #...and so forth...
   #end for
   for item in curbgz:
      bgz[nradar].append(item)
   #end for
   for item in cursprdanavr:
      sprdanavr[nradar].append(item)
   #end for
   for item in cursprdanaz:
      sprdanaz[nradar].append(item)
   #end for
   for item in cursprdbgvr:
      sprdanavr[nradar].append(item)
   #end for
   for item in cursprdbgz:
      sprdanaz[nradar].append(item)
   #end for
   for item in curallvr:
      allvr[nradar].append(item)
   #end for
   for item in curallz:
      allz[nradar].append(item)
   #end for
   for item in curallsprdvr:
      allsprdvr[nradar].append(item)
   #end for
   for item in curallsprdz:
      allsprdz[nradar].append(item)
   #end for
   nradar = nradar + 1
#end for (radar)

#-------------------------------------------#
#Print output block -- use for debugging!
#for counter in arange(0, nradar, 1):
#   print 'RMS Analysis Vr innovation for', radarlist[counter]
#   print anavr[counter]
#   print 'RMS Background Vr innovation for', radarlist[counter]
#   print bgvr[counter]
#   print 'RMS Analysis Z innovation for', radarlist[counter]
#   print anaz[counter]
#   print 'RMS Background Z innovation for', radarlist[counter]
#   print bgz[counter] 
#   print '-----------------------------------------'
#end for (counter)
#print 'The times are:'
#print times
print 'Plotting is underway'
#--------------------------------------------#

pt = plt.figure()
#Finally, create a multi-panel plot of the data
pt.text(.5, .95, 'RMS Innovation for CNTL (using '+ str(mult_inflat) + '% covariance inflation)', horizontalalignment = 'center', size = 24)
for counter in arange(0, 9, 1):
   #each iteration of the 'for' loop creates a single subplot
   plt.subplot(3,3,counter + 1)
   plt.xticks(arange(3600,8400,1200))
   plt.yticks(arange(0,21,3))
   plt.plot(times, allvr[counter], 'go-', label = 'Vr', markersize = 4)
   plt.plot(times, allz[counter], 'cv-', label = 'Z', markersize = 4)
   plt.plot(times, allsprdvr[counter], 'r--', label = 'spread (Vr)')
   plt.plot(times, allsprdz[counter], 'y--', label = 'spread (Z)')
   plt.axis([3600,7200,0,18])
   if counter == 2:
      lgnd = plt.legend(loc=(0.85, 0.75))  #We only need one legend, thanks.
      for t in lgnd.get_texts():
         t.set_fontsize('small') #make legend text small
      #end for
   #end if
   plt.title(radarlist[counter])
#end for

plt.show()
