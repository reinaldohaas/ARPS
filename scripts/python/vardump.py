#!/usr/local/bin/python2.4

#Code: vardump.py (Python Script)
#Author: Nate Snook (CASA/CAPS/SoM)
#Written: May 2009
#Last modified:  Dec. 18, 2009
#
#Purpose:
#   This simple code dumps the contents of a given variable from an ARPS history file.

from PyNIO import Nio
from plot_config import *
from numpy import zeros
from subprocess import call
import sys

if len(sys.argv) == 2:  #if the user entered a variable on the command line...
   varname = sys.argv[1]
elif len(sys.argv) == 3:
   hdfhistory = sys.argv[1]
   varname = sys.argv[2]
elif len(sys.argv) > 3:
   print '***FATAL ERROR: Too many arguments given!'
else:  #otherwise...
   varname = var
#end if/elif/else

#x and y coordinates for 1D output (used with zp, zpsoil)
xcoord = 88
ycoord = 55

dumpfile = Nio.open_file(hdfhistory, mode = 'r', options = None, history='', format='hdf')

try:  #If the variable is one of the exact ones stored in the history dump file...
   try:
      thevariable = dumpfile.variables[varname][:] #...simply gather data from that variable.
   except:
      thevariable = dumpfile.attributes[varname] #It might also be an attribute we can read.
except:  #If not, we may still be able to do something...
   if (varname == 'AGL') or (varname == 'agl'):  #If the user wants height above ground level (AGL)...
      sfc_height = dumpfile.variables['zpsoil'][0,:,:]  #...calculate it from zpsoil[sfc] and zp.
      zp = dumpfile.variables['zp'][:,:,:]
      agl = zp - sfc_height   
      print '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      for level in arange(0,len(zp[:]),1):
         print 'Level ' + str(level + 1) + ' is between ' + str(agl[level,:,:].min()) + ' and ' + str(agl[level,:,:].max()) + 'm AGL.'
      #end for
      print '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      thevariable = 'Um... look up.'
   #end if
   elif (varname == 'ptprt') or (varname == 'ptpert'):  #If the user wants potential temperature perturbation...
      pt = dumpfile.variables['pt'][:,:,:]  #...then calculate it from potential temperature (pt) and ptbar.
      ptbar = dumpfile.variables['ptbar'][:,:,:]
      ptprt = zeros((len(pt[:]),len(pt[1][:]), len(pt[1][1][:])))
      ptprt[:,:,:] = pt[:,:,:] - ptbar[:,:,:]    #perturbation = actual value - average value
      thevariable = ptprt
   else:
      print ' '
      print '*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
      print '*   Variable "' + varname + '" is not supported.  Stop goofing around and try again!   *'
      print '*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
      print ' '
      print 'Valid options are:'
      call('./varlist.py')
   #end if/else
#end try/except
      
if varname == 'zp':  #special case for zp to give 1D vertical output
   print "dimension 1 -- " + str(len(thevariable[:]))
   print "dimension 2 -- " + str(len(thevariable[1][:]))
   print "dimension 3 -- " + str(len(thevariable[1][1][:]))
   print "Providing 1D VERTICAL output at x,y = (" + str(xcoord) + ',' + str(ycoord) + ')...'
   tempvar = thevariable
   thevariable = zeros((len(thevariable[:])))
   thevariable = tempvar[:,xcoord,ycoord]   
#end if

if varname == 'zpsoil':  #special case for zpsoil to give 1D vertical output
   print "dimension 1 -- " + str(len(thevariable[:]))
   print "dimension 2 -- " + str(len(thevariable[1][:]))
   print "dimension 3 -- " + str(len(thevariable[1][1][:]))
   print "Providing 1D VERTICAL output at x,y = (" + str(len(thevariable[1][:])) + ',' + str(len(thevariable[1][1][:])) + ')...'
   tempvar = thevariable
   thevariable = zeros((len(thevariable[:])))
   thevariable = tempvar[:,xcoord,ycoord]   
#end if
    
print 'Contents of the variable "' + str(varname) + '" in "' + str(hdfhistory) + '":' 
print '-----------------------------'
print thevariable[:]
print '-----------------------------'
