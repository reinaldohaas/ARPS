#!/usr/local/bin/python2.4

#Code: varlist.py (Python Script)
#Author: Nate Snook (CASA/CAPS/SoM)
#Written: Mar. 2008
#Last modified:  Jun. 17, 2008
#
#Purpose:
#   This simple script lists the attributes and variables stored in an ARPS history dump.

from PyNIO import Nio
from plot_config import *
import sys

if len(sys.argv) > 1:  #if the user entered a variable on the command line...
   hdfhistory = sys.argv[1]

testfile = Nio.open_file(hdfhistory, mode = 'r', options = None, history='', format='hdf')

varnames = testfile.variables.keys()
print '--------------------------------------'
print 'Listing of variables stored in', hdfhistory
print varnames
print '--------------------------------------'

attrnames = testfile.attributes.keys()
print 'Listing of attributes stored in', hdfhistory
print attrnames
print '--------------------------------------'
