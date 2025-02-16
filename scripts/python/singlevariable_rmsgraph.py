#!/usr/bin/env python2.7

#Code: singlevariable_rmsgraph.py (Python Script)
#Author: Nate Snook (CASA/CAPS/SoM)
#Written: Jun. 2008
#Last modified:  Aug. 6, 2008
#
#Purpose:
#   This script creates a multi-panel plot of the RMS statistics for a each radar in
#   an ARPS EnKF data assimilation experiment.  Plotted are the RMS innovation of Vr or Z,
#   for multiple experiments (for example, using different covariance inflation factors).

import os
import subprocess as sp
from numpy import *
import matplotlib.pyplot as plt
from plot_config import *

#Initialize rms arrays
ana = []  #RMS innovation (or error) of analysis for the variable of choice
bg = []   #RMS innovation (or error) of background (forecast) for the variable of choice
sprdana = [] #RMS spread of the variable in the analysis
sprdbg = []  #RMS spread of the variable in the background
allinnov = []  #Combined ana and bg for plotting purposes
allsprd = [] #Combined sprdana and sprdbg for plotting purposes
times = []  #A listing of times for use in plotting

#NOTE: Be sure to enclose your numbers in quotes to make them strings for 'dir_names'!
#In what directories are your experiment files stored?  One directory per experiment!

#dir_names = ['2km_CNTL', '2km_88d']
#experiment_names = ['CNTL', '88D']
dir_names = ['lin_members/rms', 'wsm_members/rms', 'nem_members/rms']  #What experiments do you want to include in the plot?
experiment_names = ['Lin Members', 'WSM Members', 'NEM Members'] #And what are their names?
casadomain = False  #We are using the special CASA domain only data!
show88donly = True #We *only* want to show KTLX and KVNX for the CASA domain data!


#List the directory where RMS innovation data is stored
homedir = '/data5/nsnook/2km_CNTL025/'

#Switches for Vr and Z plotting:  1 = on, 0 = off, anything else = stop goofing around!
z_on = 0     #Shall we plot a line for radar reflectivity (Z)?
vr_on = 1    #Shall we plot a line for radial velocity (Vr)?

#List all the radars to analyze RMS stats for
if (casadomain == True and show88donly != True):
   radarlist = ['KTLX', 'KVNX', 'KCYR', 'KLWE', 'KRSP', 'KSAO']
   if colorcode == 1:
      plotbgs = ['papayawhip', 'papayawhip', 'lightcyan', 'lightcyan', 'lightcyan', 'lightcyan']
   else:
      plotbgs = ['white', 'white', 'white', 'white', 'white', 'white']
   numradars = 6
elif (casadomain == True and show88donly == True):
   numradars = 2
   radarlist = ['KTLX', 'KVNX']
   plotbgs = ['white', 'white']
elif casaon == 1:
   radarlist = ['KAMA', 'KDYX', 'KLBB', 'KTLX', 'KVNX' , 'KCYR', 'KLWE', 'KRSP', 'KSAO']
   if colorcode == 1:
      plotbgs = ['papayawhip', 'papayawhip', 'papayawhip', 'papayawhip', 'papayawhip', 'lightcyan', 'lightcyan', 'lightcyan', 'lightcyan']
   else:
      plotbgs = ['white','white','white','white','white','white','white','white','white']
   numradars = 9
else:
   radarlist = ['KAMA', 'KDYX', 'KLBB', 'KTLX', 'KVNX']
   if colorcode == 1:
      plotbgs = ['papayawhip','papayawhip','papayawhip','papayawhip','papayawhip']
   else:
      plotbgs = ['white','white','white','white','white']
   numradars = 5
#end if (casaon)

#------------------------------------------------#
# Special code section for response to reviewers #
#------------------------------------------------#

#Use KTLX and KVNX only!
#radarlist = ['KTLX', 'KVNX']
#if colorcode == 1:
#   plotbgs = ['papayawhip','papayawhip']
#else:
#   plotbgs = ['white','white']
#numradars = 2

#----------END SPECIAL CODE SECTION--------------#

#List all the times that you wish to analyze stats for
timelist = ['003900', '004200', '004500', '004800', '005100', '005400', '005700', 
            '006000', '006300', '006600', '006900', '007200']

#What value should we use for missing data?
missing = NaN

#Are we using data with differing assimilation cycle windows?  1 = yes, anything else = no
differing_cycles = 0

#initialize experiment counter variable
npct = 0

for dirname in dir_names:

   #(re-)initialize radar counter variable
   nradar = 0

   #Create a new dimension in the matricies to store data for this experiment.
   #Each array will be of the form: variable[experiment,radar,times]
   ana.append([])
   bg.append([])
   sprdana.append([])
   sprdbg.append([])
   allinnov.append([])
   allsprd.append([])
   
   for radar in radarlist:
      curana = []         #values of the RMS innovation of the analysis
      curbg = []          #values of the RMS innovation of the background
      curallinnov = []    #Background and analysis innovation of variable
      cursprdbg = []      #RMS spread of the variable for the background
      cursprdana = []     #RMS spread of the variable for the analysis
      curallsprd = []     #Background and analysis spread of variable
      
      for time in timelist:
         #First we must define the filenames for the analysis and background for innovation and spread.
         anafile = homedir + dirname + '/' + radar + 'rmsdan' + time
         bgfile = homedir + dirname + '/'+ radar + 'rmsdbg' + time 
         sprdanafile = homedir + dirname + '/'+ radar + 'spreadana' + time
         sprdbgfile = homedir + dirname + '/'+ radar + 'spreadfcs' + time
         
         try:
            thefile = open(bgfile, 'r')
            print 'OK for ' + radar + ' at ' + time + ' for ' + dirname
            vr_in = float(thefile.read(15))  #Read in RMS innovation of Vr...
            z_in = float(thefile.read(15))   #...and Z for this radar for the background.
         except IOError:
            print 'File not found: ' + bgfile + '!'
            vr_in = missing                  #If this file does not exist, set the value...
            z_in = missing                   #...to 'missing' for both Vr and Z.
         #end try
         if var == 'Vr':
            curbg.append(vr_in)
            curallinnov.append(vr_in)
         elif var == 'Z':
            curbg.append(z_in)
            curallinnov.append(z_in)
         else:
            print 'Variable type ' + var + ' is not valid.  Please choose either "Vr" or "Z".'
            exit
         #end if
         if nradar == 0 and npct == 0:    #We only need to collect the list of times once, thanks.
            times.append(float(time))
         #end if
      
         try:
            thefile = open(anafile, 'r')
            vr_in = float(thefile.read(15))   #Read in RMS innovation of Vr...
            z_in = float(thefile.read(15))    #...and Z for this radar for the analysis.
         except IOError:
            vr_in = missing                   #If this file does not exist, set the value...
            z_in = missing                    #...to 'missing' for both Vr and Z.
         #end try
         if var == 'Vr':
            curana.append(vr_in)
	    curallinnov.append(vr_in)
	 elif var == 'Z':
            curana.append(z_in)
	    curallinnov.append(z_in)
	 else:
            print 'Variable type ' + var + ' is not valid.  Please choose either "Vr" or "Z".'
	    exit
	 # end if
         if nradar == 0 and npct == 0:    #We only need to collect the list of times once, thanks.
            times.append(float(time))
         # end if

         try:
            thefile = open(sprdbgfile, 'r')

            vr_in = float(thefile.read(15))   #Read in RMS innovation of Vr...
            z_in = float(thefile.read(15))    #...and Z for this radar for the analysis.

            
            print 'BG Spread:  vr_in = ' + str(vr_in) + '|| z_in = ' + str(z_in)
         except IOError:
            try:  #If using 'postinnov' rather than 'arpsenkf', we need to use a different file...
               print 'No "spreadfcs" files found, trying "innov"...'
               sprdfcsfile = homedir + dirname + '/'+ radar + 'innov' + time
               thefile = open(sprdfcsfile, 'r')
               
               print 'For file: ' + sprdbgfile + ':'
               vr_in = float(thefile.read(15))   #Read in RMS spread of Vr...
               z_in = float(thefile.read(15))    #...and Z for this radar for the analysis.
               print 'BG Spread:  vr_in = ' + str(vr_in) + '|| z_in = ' + str(z_in)
            except IOError:
               print 'ERROR -- values are MISSING for fcst spread!'
               print 'File ' + sprdfcsfile + ' not found!'
               vr_in = missing                   #If this file does not exist, set the value...
               z_in = missing                    #...to 'missing' for both Vr and Z.
            # end nested try
         #end outer try
         if var == 'Vr':
            cursprdbg.append(vr_in)
	    curallsprd.append(vr_in)
	 elif var == 'Z':
            cursprdbg.append(z_in)
	    curallsprd.append(z_in)
	 else:
            print 'Variable type ' + var + ' is not valid.  Please choose either "Vr" or "Z".'
	    exit
	 # end if

         try:
            thefile = open(sprdanafile, 'r')
            vr_in = float(thefile.read(15))   #Read in RMS spread of Vr...
            z_in = float(thefile.read(15))    #...and Z for this radar for the analysis.
         except IOError:
            vr_in = missing                   #If this file does not exist, set the value...
            z_in = missing                    #...to 'missing' for both Vr and Z.
         #end try
         if var == 'Vr':
            cursprdana.append(vr_in)
	    curallsprd.append(vr_in)
	 elif var == 'Z':
            cursprdana.append(z_in)
	    curallsprd.append(z_in)
	 else:
            print 'Variable type ' + var + ' is not valid.  Please choose either "Vr" or "Z".'
	    exit
	 # end if
      # end for (time)
      print '-  -  -  -  -  -  -  -  -  -  -  -  -  -  -'
       
      ana[npct].append([])
      bg[npct].append([])
      sprdana[npct].append([])
      sprdbg[npct].append([])
      allinnov[npct].append([])
      allsprd[npct].append([])
   
      for item in curana:
         ana[npct][nradar].append(item)   #Store Vr data into the anavr array
      #end for
      for item in curbg:
         bg[npct][nradar].append(item)    #...and so forth...
      #end for
      for item in cursprdana:
         sprdana[npct][nradar].append(item)
      #end for
      for item in cursprdbg:
         sprdbg[npct][nradar].append(item)
      #end for
      for item in curallinnov:
         allinnov[npct][nradar].append(item)
      #end for
      for item in curallsprd:
         allsprd[npct][nradar].append(item)
      #end for
      nradar = nradar + 1   #iterate radar counter
   #end for (radar)
   npct = npct + 1   #iterate experiment counter
#end for (dirname)
#print 'len(timelist) = ' + str(len(timelist))
print 'allinnov[1][0][:] = ' + str(allsprd[1][0][:])
#print 'allsprd[:][0][0] = ' + str(allsprd[:][0][0])
#print 'allsprd[1][:][0] = ' + str(allsprd[1][:][0])

#Now, deal with 'missing' data for the case of uneven assimilation periods.
#We'll use a simple averaging algorithm to do this.

if differing_cycles == 1:

   nexp = 0  #reset experiment iterator
   for dirname in dir_names:
      nradar = 0   #reset radar iterator
      for radar in radarlist:
         ntime = 0  #reset time iterator
         for itime in times:
            if ntime < 2:
               ntime = ntime + 1
               continue
            #end if
            if ntime > len(times) - 2:
               ntime = ntime + 1
               continue
            #end if
            if not (allinnov[nexp][nradar][ntime] > -99999.0 and allinnov[nexp][nradar][ntime] < 99999.0):
               print 'allinnov[' + str(nexp) + '][' + str(nradar) + '][' + str(ntime) + '] = ' + str(allinnov[nexp][nradar][ntime])
               allinnov[nexp][nradar][ntime] = (allinnov[nexp][nradar][ntime - 1] + allinnov[nexp][nradar][ntime + 2]) / 2
            
               allinnov[nexp][nradar][ntime + 1] = allinnov[nexp][nradar][ntime]
            #end if
            if not (allsprd[nexp][nradar][ntime] > -99999.0 and allinnov[nexp][nradar][ntime] < 99999.0):
               print 'allsprd[' + str(nexp) + '][' + str(nradar) + '][' + str(ntime) + '] = ' + str(allsprd[nexp][nradar][ntime])
               allsprd[nexp][nradar][ntime] = (allsprd[nexp][nradar][ntime - 1] + allsprd[nexp][nradar][ntime + 2]) / 2
               allsprd[nexp][nradar][ntime + 1] = allsprd[nexp][nradar][ntime]
            #end if
            ntime = ntime + 1
         #end for (time)
         nradar = nradar + 1
      #end for (radar)
      nexp = nexp + 1
   #end for (dirname)
#end if (differing_cycles)

#Mask the data for the missing values:
allsprd_masked = ma.masked_where(isnan(allsprd),allsprd)
allinnov_masked = ma.masked_where(isnan(allinnov),allinnov)

print 'allinnov[1][0][:] = ' + str(allsprd[1][0][:])

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
print 'The times are:'
print times
print 'Plotting is underway'
#--------------------------------------------#

#Define colors and styles of lines used in the plots (see pyplot documentation for more details)
if color_plots == 1:
   plt_colors = ['blue', '#4444FF', 'red', '#FF4444', 'green', '#44FF44', 'black', '#777777']
   plt_linestyles = ['-', '--', '-', '--', '-', '--', '-', '--']
#   lineformat = ['b-', 'c--', 'r-', 'm--', 'g-', 'y--', 'k-', 'b--']
else:
   plt_colors = ['#000000', '#000000', '#888888', '#888888', '#bbbbbb', '#bbbbbb', '#dddddd', '#dddddd']
   plt_linestyles = ['-', '--', '-', '--', '-', '--', '-', '--']
#end if/else

#Adds support for transparency effects -- transparency is the transparency level (0 = fully transparent, 1 = opaque)
plt_transparency = [1.0, 1.0, transparency, transparency, transparency, transparency, transparency, transparency]

pt = plt.figure()
#Finally, create a multi-panel plot of the data

#Add a title, if desired.
if add_title == 1:
   if vr_on == 1 and z_on == 1:
      pt.text(.5, .95, 'RMS Innovation and Spread of ' + var, horizontalalignment = 'center', size = 24)
   elif vr_on != 1 and z_on == 1:
      pt.text(.5, .95, 'RMS Innovation and Spread of ' + var, horizontalalignment = 'center', size = 24)
   elif vr_on == 1 and z_on != 1:
      pt.text(.5, .95, 'RMS Innovation and Spread of ' + var, horizontalalignment = 'center', size = 24)
   else:
      pt.text(.5, .95, 'RMS Innovation and Spread of absolutely nothing -- try plotting a variable!', horizontalalignment = 'center', size = 24)
#end if

for counter in range(0, numradars, 1):
   #each iteration of the 'for' loop creates a single subplot
   if casadomain == True and show88donly == False:
      plt.subplot(3,2,counter + 1, axisbg = plotbgs[counter])
   elif casadomain == True and show88donly == True:
      plt.subplot(2,1,counter + 1, axisbg = plotbgs[counter])
   elif casaon == 0:
      plt.subplot(3,2,counter + 1, axisbg = plotbgs[counter])
   else:
      plt.subplot(3,3,counter + 1, axisbg = plotbgs[counter])
   #end if
   plt.xticks(arange(3600,8400,1200))
   plt.yticks(arange(0,y_max + y_step,y_step))
   innercount = 0
   linecount = 0
   for experiment in experiment_names:
      plt.plot(times, allinnov_masked[innercount][counter], color = plt_colors[linecount], linewidth = 1.75,
               linestyle = plt_linestyles[linecount], alpha = plt_transparency[linecount], label = 'Innovation' + ' (' + experiment + ')' )
      plt.plot(times, allsprd_masked[innercount][counter], color = plt_colors[linecount + 1], linewidth = 1.75,
               linestyle = plt_linestyles[linecount + 1], alpha = plt_transparency[linecount + 1], label = 'Spread' + ' (' + experiment + ')' )
#      plt.plot(times, allinnov_masked[innercount][counter], lineformat[linecount], 
#               label = 'Innovation' + ' (' + experiment + ')' )
#      plt.plot(times, allsprd_masked[innercount][counter], lineformat[linecount + 1], 
#               label = 'Spread ' + ' (' + experiment + ')')
      innercount = innercount + 1
      linecount = linecount + 2
   #end for
   plt.axis([3600,7200,y_min,y_max])

   if casadomain == True and show88donly == False:
      #For 6 radars covering the CASA domain
      if counter == 1:
         lgnd = plt.legend(loc=(0.64,0.58))  #Useful for plotting just inside upper-right plot
#         lgnd = plt.legend(loc=(0.65,0.57))   #Useful if we need to lower the top of the legend a bit.
#         lgnd = plt.legend(loc=(0.65, 0.70))  #We only need one legend, thanks.  (DEFAULT LEGEND POSITION)
         for t in lgnd.get_texts():
            t.set_fontsize('small') #make legend text small
         #end for
      #end if (counter)
      
   if casadomain == True and show88donly == True:
      #For 2 KTLX and KVNX covering the CASA domain
      if counter == 0:
         lgnd = plt.legend(loc=(0.74,0.68))  #Useful for plotting just inside upper-right plot
#         lgnd = plt.legend(loc=(0.65,0.57))   #Useful if we need to lower the top of the legend a bit.
#         lgnd = plt.legend(loc=(0.65, 0.70))  #We only need one legend, thanks.  (DEFAULT LEGEND POSITION)
         for t in lgnd.get_texts():
            t.set_fontsize('small') #make legend text small
         #end for
      #end if (counter)
#-----------------------------
   elif casaon == 1:
      #For 9 radars
      if counter == 2:
         lgnd = plt.legend(loc=(0.64,0.58))  #Useful for plotting just inside upper-right plot
#         lgnd = plt.legend(loc=(0.65,0.57))   #Useful if we need to lower the top of the legend a bit.
#         lgnd = plt.legend(loc=(0.65, 0.70))  #We only need one legend, thanks.  (DEFAULT LEGEND POSITION)
         for t in lgnd.get_texts():
            t.set_fontsize('small') #make legend text small
         #end for
      #end if (counter)
#-----------------------------
   else:
      #For 5 radars
      if counter == 4:
         lgnd = plt.legend(loc=(1.12,0.35))  #We only need one legend, thanks.
#         for t in lgnd.get_texts():
#            t.set_fontsize('small') #make legend text small
         #end for
      #end if (counter)
   #end if (casaon)
   plt.title(radarlist[counter])
#end for

#For 5 radars
plt.subplots_adjust(left=0.04, right=0.97, bottom=0.04, top=0.91, wspace=0.12, hspace=0.20)

plt.show()
