#!/bin/env python
## ---------------------------------------------------------------------
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## ---------------------------------------------------------------------
##
## This is a python program that converts ARPS history files to VDF (VAPOR
## data format).
##
## It will generate a .vdf metafile and a VAPOR data collection
## (VDC) from a given set ARPS history files.
##
## ---------------------------------------------------------------------
##
## HISTORY:
##   Yunheng Wang (09/15/2011)
##   Upgraded to support VAPOR v2.0.2.
##
##   Yunheng Wang (10/24/2011)
##   Fixed a bug with dimension order for terrain height. Detected by Alex Schenkman.
##
##############################################################################
##
## Requirements:
##   numpy            1.2.0 or later
##   ScientificPython 2.9.0 or later
##
##############################################################################

import sys, os, re, getopt, subprocess, time
import math

## Import the necessary numpy and netCDF modules
from numpy import *
from Scientific.IO.NetCDF import *

##======================================================================
## USAGE
##======================================================================
def usage( istatus = 0 ) :
  '''-------------------------------------------------------------------
  Print Usage and then exit
  -------------------------------------------------------------------'''

  version  = '1.1'
  lastdate = '2011.10.24'

  print >> sys.stderr, """\n  Usage: %s [options] ARPS_runname Start_time End_time Time_interval\n
     \tARPS_runname\t ARPS runname for the file name base
     \tStart_time\t ARPS history file start time in seconds
     \tEnd_time\t ARPS history file end time in seconds
     \tTime_interval\t Time interval between ARPS files
     \n  OPTIONS:\n
     \tOption \t\tDefault \tInstruction
     \t-------\t\t--------\t-----------
     \t-v, --verbose      \t \tVerbose
     \t-g, --geo          \t \tPrint georeference information for inserting images
     \t-h, --help         \t \tPrint this help
     \t-l, --level\t0     \t \tMaximum refinement level
     \t-m, --map_name\tBMNG \t \tPredefined map from a well-known server

     \t                BMNG       \tNASA BlueMarble, the default
     \t                landsat    \tLandsat imagery
     \t                USstates   \tUS state boundaries
     \t                UScounties \tUS state and county boundaries
     \t                world      \tworld political boundaries
     \t                rivers     \tmajor rivers
     \t                none       \tDo not fetch terrain image

     \t-t, --transparent     \t \tRequest a transparent background
     """ % os.path.basename(cmd)
  print >> sys.stderr, '  For questions or problems, please contact yunheng@ou.edu (v%s, %s).\n' % (version,lastdate)

  sys.exit(istatus)
#enddef

##======================================================================
## Parse command line arguments
##======================================================================
def parseArgv(argv) :
  '''-------------------------------------------------------------------
  Parse command line arguments
  -------------------------------------------------------------------'''
  global _debug

  argsdict = dict()
  fetchdict = dict()

  try:
    opts, args = getopt.getopt(argv,'hvgl:tm:',           \
           ['help','verbose','geo','level','transparent','map_name'])

    for opt, arg in opts :
      if opt in ('-h','--help'):
        usage(0)
      elif opt in ( '-v', '--verbose'):
        _debug = True
      elif opt in ( '-g', '--geo'):
        argsdict['geofile']=True
      elif opt in ( '-l', '--level'):
        if int(arg) in range(0,100) :
          argsdict['reflevel'] = int(arg)
        else :
          print >> sys.stderr, 'ERROR: refinement level (%s) out of range.' % arg
          usage(3)

      elif opt in ( '-t'):
        fetchdict['-t'] = ''
      elif opt in ( '-m', '--map_name'):
        if (arg.lower() == 'none') :
          argsdict['fetchterrain'] = False
        elif arg in ('BMNG','landsat','USstates','UScounties','world','rivers') :
          fetchdict['-m']=arg
        else :
          print >> sys.stderr, 'ERROR: Unknown argument (%s) for option "-m"' % arg
          usage(3)

  except getopt.GetoptError:
    print >> sys.stderr, 'ERROR: Unknown option.'
    usage(2)

  if (len(args) != 4 ) :
    print >> sys.stderr, 'ERROR: wrong number of command line arguments.'
    usage(1)
  else :
    argsdict['arpsrunname']=args[0]

    if re.match('^\d+$',args[1]) :
      argsdict['istart'] = int(args[1])
    else :
      print >> sys.stderr, 'ERROR: wrong command line argument %s.' % args[1]
      usage(2)

    if re.match('^\d+$',args[2]) and args[2] >= args[1] :
      argsdict['iend'] = int(args[2])
    else :
      print >> sys.stderr, 'ERROR: wrong command line argument %s.' % args[2]
      usage(2)

    if re.match('^\d+$',args[3]) :
      argsdict['istep'] = int(args[3])
    else :
      print >> sys.stderr, 'ERROR: wrong command line argument %s.' % args[3]
      usage(2)

  argsdict['fetchOpts']=fetchdict
  return argsdict
#enddef

##======================================================================
## Main entry point
##======================================================================
def main(arpsrunname='arps25may1998',istart=0,iend=0,istep=0,
         reflevel=0,geofile=False,
         fetchterrain=True,fetchOpts=dict()):
  '''-------------------------------------------------------------------
    Main entry point
  ------------------------------------------------------------------'''
  numts = (iend-istart)/istep+1
  ts = 0
  for itime in range(istart,iend+istep,istep) :
    arpsfile = '%s.net%06d' %(arpsrunname,itime)

    if (ts < 1) :
      print '\n=== Preparing the VDF metafile ...'

      filein = NetCDFFile(arpsfile,'r')

      ##----------------------------------------------------------------
      ## Get file dimensions
      dimstr = getDimStr(filein)

      ##----------------------------------------------------------------
      ## Get map projection
      arpsmap = getARPSMap(filein)
      projstr = arpsmap.getPROJstr()
      ##print 'projstr = %s' % projstr

      ##----------------------------------------------------------------
      ## Get coordinates
      arpsrange = getARPSRange(filein,arpsmap)
      extstr = '%.2f:%.2f:%.2f:%.2f:%.2f:%.2f' % (arpsrange.xyzmass) ## mass domain

      ##----------------------------------------------------------------
      ## Get times
      arpstime = getARPSTime(filein)
      ##print arpstime.getstr()

      ##----------------------------------------------------------------
      ## Get variable names and dimension string for each variable
      var2dList = []
      var3dList = []
      vardimdict = dict()
      varstr = getVarLists(filein,var2dList,var3dList,vardimdict)

      ##----------------------------------------------------------------
      ## fetch terrain image
      fetchTerrainMap(fetchterrain,arpsmap,arpsrange.llvector,arpsrunname,
                      fetchOpts)

      ##----------------------------------------------------------------
      ## Generate steady terrain field
      hgtfname = arpsrunname+ '_HGT.nc'
      print '\n--- Preparing terrain height in %s ...' %hgtfname
      derive_hgt(filein,hgtfname)
      var2dList.append('HGT')
      vardimdict['HGT'] = 'x:y:dummy'

      filein.close()

      ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ## run vdfcreate
      ##
      vdfcreate(arpstime,numts,istart,istep,dimstr,projstr,extstr,
                var2dList,var3dList,arpsrunname,reflevel)

      if geofile :
        geocreate(arpstime,istart,iend,istep,projstr,arpsrunname,arpsrange)

    print '\n--- Processing file "%s" ... ' % arpsfile

    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ## run ncdf2vdf
    ##
    ncdf2vdf(varstr,vardimdict,ts,arpsrunname,arpsfile,reflevel)

    ncdf2vdf('HGT',vardimdict,ts,arpsrunname,hgtfname,reflevel,
             contain_time = False,delete_file = False)

    ts += 1

  ##--------------------------------------------------------------------
  ## Clear temporary files
  os.unlink(hgtfname)

  return
#enddef main

########################################################################

def vdfcreate (arpstime,numts,istart,istep,dimstr,projstr,extstr,
               var2dList,var3dList,runname,reflevel) :
    '''----------------------------------------------------------------
     run vdfcreate

    vdfcreate -dimension 66x66x42 -mapprojection "+proj=lcc +lat_1=30 +lat_2=60 lon_0=-100 +lat_0=38 +ellps=sphere"
    -gridtype layered -extents -16:-16:-250:2074:2074:20250 -vars3d U:V:W:PT arps25may1998.vdf
    ----------------------------------------------------------------'''

    tmpvdf = '%s_%d.vdf' % (runname,os.getpid())

    vdfcreatecmd = 'vdfcreate -numts %d -startt %d -deltat %d ' %       \
                   (numts,istart,istep)
    vdfcreatecmd += '-dimension %s -mapprojection %s -gridtype layered -extents %s -vars3d %s -vars2dxy %s' % \
                   (dimstr,projstr,extstr,':'.join(var3dList),':'.join(var2dList))
    if reflevel > 0 : vdfcreatecmd += ' -level %d' % reflevel
    vdfcreatecmd += ' %s' % tmpvdf

    print vdfcreatecmd
    vdfP = subprocess.Popen(vdfcreatecmd,shell=True)
    vdfP.communicate()
    ##os.system(vdfcreatecmd)

    ##
    ## Now, modify the UserTimeStampString in the .vdf file
    ##
    if vdfP.returncode == 0 :
      tmpf = open(tmpvdf,'r')
      vdfF = open('%s.vdf' % runname,'w')
      ropen = False
      iopen = False
      itime = istart
      for line in tmpf :
        if line.find('<UserTimeStampString Type="String">') >= 0:
          ropen = True
        if line.find('</UserTimeStampString>') >=0 :
          ropen = False

        if line.find('<UserTime Type="Double">') >= 0:
          iopen = True
        if line.find('</UserTime>') >=0 :
          iopen = False

        if ropen and line.strip().isdigit() :
          line = '  %s\n' % arpstime.getstr(itime)
          itime += istep
        elif iopen and line.strip().isdigit() :
          line = '  %d\n' % arpstime.getsec(itime)
        vdfF.write(line)


      vdfF.close()
      tmpf.close()

    os.unlink(tmpvdf)

#enddef vdfcreate

def ncdf2vdf(varstr,dimstrs,ts,arpsrunname,arpsfile,reflevel,
             contain_time = True,delete_file = False) :
    '''----------------------------------------------------------------
     run ncdf2vdf

     ncdf2vdf -ts 0 -varname ELEVATION -ncdfvar ZP -dimnames "x:y:z_stag"
     -cnstnames Time -cnstvals 0 arps25may1998.vdf arps25may1998.net003600
    ----------------------------------------------------------------'''

    global _debug

    for varName in varstr.split(':') :
      if varName == 'ZP' : vapName = 'ELEVATION'
      else : vapName = varName
      ##var = file.variables[varName]
      ##dimstr = ':'.join(reversed(var.dimensions[1:]))

      ncdf2vdfcmd = 'ncdf2vdf -ts %d -varname %s -ncdfvar %s -dimnames "%s" ' % \
                    (ts,vapName, varName, dimstrs[varName])
      if contain_time : ncdf2vdfcmd += '-cnstnames Time -cnstvals 0 '
      if reflevel > 0 : ncdf2vdfcmd += "-level %d " % reflevel
      ncdf2vdfcmd += '%s.vdf %s' % ( arpsrunname, arpsfile)

      print ncdf2vdfcmd

      ##os.system(ncdf2vdfcmd)
      ncdf2vdfP = subprocess.Popen(ncdf2vdfcmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      (n2vout,n2verr) = ncdf2vdfP.communicate()

      if ncdf2vdfP.returncode != 0:
        sys.stdout.flush()
        print >> sys.stderr, '%s : %s %s' % (arpsfile, varName, n2verr)
      else :
        for outline in n2verr.splitlines() :
          if outline.find('staggered') > 0 or outline.find('not match') > 0 :
            print '%s : %s %s ' % (arpsfile, varName, outline)
            sys.stdout.flush()

          if _debug and outline.find('min and max values') > 0 :
            print '%s : %s %s ' % (arpsfile, varName, outline)
            sys.stdout.flush()

    if delete_file : os.unlink(arpsfile)

    return
#enddef ncdf2vdf

def geocreate(arpstime,istart,iend,istep,projstr,arpsrunname,arpsrange) :

  geocmd = 'tiff2geotiff -4 %s -m %s_geo.txt %s.tiff %s_geo.tiff' %(
           projstr,arpsrunname,arpsrunname,arpsrunname)

  gfile = open('%s_geo.txt' % arpsrunname,'w')

  for itime in range(istart,iend+istep,istep) :
    line = arpstime.getstr(itime)
    line += ' %.2f %.2f %.2f %.2f'% arpsrange.llmass
    line += ' 0.2 0.315 0.8 0.73\n'
    gfile.write(line)
  gfile.close()

  print '\n*** 2DXY images should be converted with the following command:'
  print '%s' % geocmd

#enddef geocreate

def derive_hgt(arpsfile,filename) :
  '''
  derive terrain heights in an intermediate file for applying image to
  in VAPOR terrain
  '''
  hgtfile = NetCDFFile(filename,'w')

  hgtfile.createDimension('x',arpsfile.dimensions['x'])
  hgtfile.createDimension('y',arpsfile.dimensions['y'])

  hgtvar = hgtfile.createVariable('HGT','f',('y','x'))

  zpvar = arpsfile.variables['ZP']

  hgtvar.assignValue(zpvar[2,:,:])

  hgtfile.sync()

  hgtfile.close()

#enddef derive_hgt


########################################################################

##def reversed(aList) :
##    '''-----------------------------------------------------------------
##      To replace the reversed builtin that is only available after
##      Python 2.4
##    -----------------------------------------------------------------'''
##    size = len(aList)
##    bList=[]
##    for i in range(size) :
##      #print i,size,size-i-1
##      bList.append(aList[size-i-1])
##
##    return bList

def getDimStr(file) :
    '''-----------------------------------------------------------------
     Get file dimensions
    -----------------------------------------------------------------'''
    fileDims = file.dimensions
    dimstr = '%dx%dx%d' % (fileDims['x'],fileDims['y'],fileDims['z'])
    ##print dimstr
    return dimstr
#enddef getDimStr

def getARPSMap(filein) :
    '''-----------------------------------------------------------------
     Get map projection

     Note that the origin is set in the central lat/lon of the model domain.
    ----------------------------------------------------------------'''
    mapproj = getattr(filein,'MAPPROJ')
    lat_1   = getattr(filein,'TRUELAT1')
    lat_2   = getattr(filein,'TRUELAT2')
    trulon  = getattr(filein,'TRUELON')
    scale   = getattr(filein,'SCLFCT')

    arpsmap = mapprojection(mapproj[0],scale[0],lat_1[0],lat_2[0],trulon[0])

    ctrlat  = getattr(filein,'CTRLAT')
    ctrlon  = getattr(filein,'CTRLON')
    ##dx = getattr(filein,'DX')
    ##dy = getattr(filein,'DY')
    ##nx = filein.dimensions['x_stag']
    ##ny = filein.dimensions['y_stag']

    ##(ctrx,ctry) = arpsmap.lltoxy(ctrlat[0],ctrlon[0])
    ##swx = ctrx - (nx-3)/2.0 * dx[0]
    ##swy = ctry - (ny-3)/2.0 * dy[0]
    ##
    ##arpsmap.setorig(1,swx,swy)
    arpsmap.setorig(2,ctrlat[0],ctrlon[0])

    return arpsmap

#enddef getARPSMap

def getVarLists(filein,var2dList,var3dList,varDim_dict) :
    '''----------------------------------------------------------------
     Get variable names and dimension string
    ---------------------------------------------------------------'''

    for varName in filein.variables.keys():
      var = filein.variables[varName]
      if (len(var.shape) == 3) :
        ##print "%s => %s" % (varName, var.shape)
        if (var.dimensions[0] == 'Time') :        ## 2D variables
          var2dList.append(varName)
          varDim_dict[varName] = ':'.join(reversed(var.dimensions[1:]))
          varDim_dict[varName] += ':dummy'
        elif (var.dimensions[0] == 'z_stag') :    ## ZP
          var3dList.append(varName)
          varDim_dict[varName] = ':'.join(reversed(var.dimensions))

      if (len(var.shape) == 4) :
        ##print "%s => %s" % (varName, var.shape)
        if (var.dimensions[0] == 'Time') and (var.dimensions[1] != 'nstyp_total') :
          var3dList.append(varName)
          varDim_dict[varName] = ':'.join(reversed(var.dimensions[1:]))

    varStr = ':'.join(var3dList+var2dList)
    return varStr

#enddef getVarLists

def getARPSRange(filein,arpsmap) :

    ##varx = file.variables['x_stag']
    ##vary = file.variables['y_stag']
    ##varz = file.variables['z_stag']
    ##x = varx.getValue()
    ##y = vary.getValue()
    ##z = varz.getValue()
    ##xs = (x[0]+x[1])/2.0
    ##ys = (y[0]+y[1])/2.0
    ##zs = (z[0]+z[1])/2.0
    ##xe = (x[-1]+x[-2])/2.0
    ##ye = (y[-1]+y[-2])/2.0
    ##ze = (z[-1]+z[-2])/2.0

    dx = getattr(filein,'DX')
    dy = getattr(filein,'DY')
    nx = filein.dimensions['x_stag']
    ny = filein.dimensions['y_stag']

    xe = (nx-3)/2.0 * dx[0] - dx[0]/2.0
    ye = (ny-3)/2.0 * dy[0] - dy[0]/2.0
    xs = -1*xe
    ys = -1*xe

    varz = filein.variables['z_stag']
    z = varz.getValue()
    zs = (z[0]+z[1])/2.0
    ze = (z[-1]+z[-2])/2.0

    arpsrange = domainrange(arpsmap,xs,xe,ys,ye,zs,ze,dx[0],dy[0])

    print
    arpsrange.printxyzmass()
    arpsrange.printllmass()
    arpsrange.printllvector()
    print

    return arpsrange

#enddef getARPSRange

def getARPSTime(filein) :

    timestr = getattr(filein,'INITIAL_TIME')
    ##print timestr
    arpstime = modeltime(timestr)

    return arpstime

#enddef getARPSTime

def fetchTerrainMap(action,arpsmap,llrange,arpsname,opts) :
    '''-----------------------------------------------------------------
     Get terrain image from well-known server
    ----------------------------------------------------------------'''
    ##(maxlat,rlon) = arpsmap.xytoll(0.0, ye)
    ##(rlat,minlon) = arpsmap.xytoll(xs,  ye)
    ##(rlat,maxlon) = arpsmap.xytoll(xe,  ye)
    ##
    ##(minlat,rlon) = arpsmap.xytoll(xs,  ys)
    ##(rlat,  rlon) = arpsmap.xytoll(xe,  ys)
    ##minlat = min(minlat,rlat)
    ##
    ##print '  ARPS domain latitude range is %.2f to %.2f \n  and longitude range is %.2f to %.2f.' % \
    ##      (minlat,maxlat,minlon,maxlon)

    if action :
      fetchcmd = 'getWMSImage.sh -o %s_terrain.tiff ' % arpsname
      for opt in opts.keys():
        fetchcmd += ' %s %s' % (opt,opts[opt])
      fetchcmd += ' %.2f %.2f %.2f %.2f' % llrange
      print '%s\n' % fetchcmd
      subprocess.Popen(fetchcmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()

#enddef fetchTerrainMap

#################### class mapprojection ###############################

class mapprojection (object) :
  '''
  General Information

  This is written from the corresponding Fortran program in the ARPS package.
  it requires to import math module.

  '''

  ##from math import *

  def __init__(self,iproj,scale,lat1,lat2,orient):
    '''
    Set constants for map projections.

    iproj        Map projection number
                 1=North Polar Stereographic   (-1 South Pole)
                 2=Northern Lambert Conformal  (-2 Southern)
                 3=Mercator
                 4=Lat,Lon

    scale        Map scale factor,  at latitude=latnot
                 Distance on map = (Distance on earth) * scale
                 For ARPS model runs, generally this is 1.0
                 For ARPS plotting this will depend on window
                 size and the area to be plotted.

    lat1/lat2    Real "True" latitude(s) of map projection
                 (degrees, positive north)
                 Except for iproj=1, only latnot(1) is used

    orient       Longitude line that runs vertically on the map.
                 (degrees, negative west, positive east)
    '''

    ## Constants

    self.ERADIUS = 6371000.             ## mean earth radius in m

    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    self.xorig=0.
    self.yorig=0.
    self.jproj=abs(iproj)
    self.jpole= 1
    if iproj < 0 : self.jpole = -1

    if  (self.jproj == 0 ) :

      print 'No map projection will be used.'

    elif ( self.jproj == 1 ) :
      ##
      ##----------------------------------------------------------------
      ##
      ##  Polar Stereographic projection
      ##  For this projection:
      ##      projc1 is the scaled earth's radius, scale times eradius
      ##      projc2 is the numerator of emfact, the map image scale factor.
      ##      projc3 is projc2 times the scaled earth's radius.
      ##
      ##----------------------------------------------------------------
      ##
      self.trulat1 = lat1
      self.rota    = orient
      self.scmap   = scale
      self.projc1  = scale*self.ERADIUS
      self.projc2  = (1. + math.sin(self.jpole*math.radians(self.trulat1)) )
      self.projc3  = projc1*projc2
      if (self.jpole > 0) :
        print '  Map projection set to Polar Stereographic\n',
        '  X origin, Y origin set to 0.,0. at the North Pole.'
      else :
        print '  Map projection set to Polar Stereographic\n',
        '  X origin, Y origin set to 0.,0. at the South Pole.'

    elif ( self.jproj == 2 ):
      ##
      ##----------------------------------------------------------------
      ##
      ##  Lambert Conformal Conic Projection.
      ##  For this projection:
      ##      projc1 is the scaled earth's radius, scale times eradius
      ##      projc2 is cos of trulat(1)
      ##      projc3 is tan (45. - trulat/2) a const for local map scale
      ##      projc4 is the cone constant, n
      ##
      ##-----------------------------------------------------------------
      ##
      self.trulat1=lat1
      self.trulat2=lat2
      self.rota=orient
      self.scmap=scale
      self.projc2=math.cos(math.radians(self.trulat1))
      self.projc3=math.tan(math.radians((45.-0.5*self.jpole*self.trulat1)))

      denom1=math.cos(math.radians(self.trulat2))
      denom2=math.tan(math.radians((45.-0.5*self.jpole*self.trulat2)))

      if (abs(self.trulat1-self.trulat2) > 0.01 and denom2 != 0.) :
        denom3=log( self.projc3/denom2 )
      else :
        denom3=0.

      if(denom1 != 0 and denom3 != 0.) :
        self.projc4=log( self.projc2/denom1 ) / denom3

        if ( self.projc4 < 0.) :
          print '  Warning in SETMAPR for Lambert Projection\n',
          '  For the true latitudes provided, %f and %f.\n' %(self.trulat1,self.trulat2),
          '  projection must be from opposite pole ... changing pole.'

          self.jpole=-1*self.jpole
          self.projc3=math.tan(math.radians((45.-0.5*self.jpole*self.trulat1)))
          denom2=math.tan(math.radians((45.-0.5*self.jpole*self.trulat2)))
          if(denom2 != 0.) :
            denom3=log( self.projc3/denom2 )
          else :
            denom3=0.

          if(denom1 != 0 and denom3 != 0.) :
            self.projc4=log( self.projc2/denom1 ) / denom3
          else :
            print '  Error (1) in SETMAPR for Lambert Projection',
            '  Illegal combination of trulats one: %f and two: %f.\n' %(self.trulat1,self.trulat2)
            sys.exit(1)


        self.projc1=scale*self.ERADIUS/self.projc4
      elif(denom3 == 0. and denom2 != 0.) :   ## tangent
        self.projc4=math.sin(math.radians(self.jpole*self.trulat1))
        if( self.projc4 < 0.) :
          print '  Warning in SETMAPR for Lambert Projection\n',
          '  For the true latitudes provided, %f and %f.\n' %(self.trulat1,self.trulat2),
          '  projection must be from opposite pole...changing pole.'
          self.jpole=-1*self.jpole
          self.projc4=math.sin(math.radians(self.jpole*self.trulat1))

        self.projc1=scale*self.ERADIUS/self.projc4
      else :
        print '  Error (1) in SETMAPR for Lambert Projection',
        '  Illegal combination of trulats one: %f and two: %f.\n' %(self.trulat1,self.trulat2)
        sys.exit(1)

      if(self.jpole > 0) :
        print '  Map projection set to Lambert Conformal\n',            \
              '  X origin, Y origin set to 0.,0. at the North Pole.'
      else :
        print '  Map projection set to Lambert Conformal\n',            \
              '  X origin, Y origin set to 0.,0. at the South Pole.'

    elif( self.jproj == 3 ) :
     ##
     ##-----------------------------------------------------------------
     ##
     ##  Mercator Projection.
     ##  For this projection:
     ##      self.projc1 is the scaled earth's radius, scale times self.ERADIUS
     ##      self.projc2 is cos of self.trulat1
     ##      self.projc3 is self.projc1 times self.projc2
     ##
     ##-----------------------------------------------------------------
     ##
      self.trulat1=lat1
      self.rota=orient
      self.scmap=scale
      self.projc1=scale*self.ERADIUS
      self.projc2=math.cos(math.radians(self.trulat1))
      self.projc3=self.projc1*self.projc2
      if(self.projc2 <= 0.) :
        print '  Error (1) in SETMAPR for Mercator Projection',
        '  Illegal true latitude provided: %f. ' % self.trulat1
        sys.exit(1)

      print '  Map projection set to Mercator',
      '  X origin, Y origin set to 0.,0. at the equator, %f.',self.rota,
      '  Y positive toward the North Pole.'

    elif( self.jproj == 4 ) :

      ##
      ##----------------------------------------------------------------
      ##
      ##  Lat, Lon Projection.
      ##  For this projection:
      ##      self.projc1 is the scaled earth's radius, scale times self.ERADIUS
      ##      self.projc2 is cos of self.trulat1
      ##      self.projc3 is self.projc1 times self.projc2 times 180/pi
      ##
      ##----------------------------------------------------------------
      ##
      self.trulat1=latnot(1)
      self.rota=orient
      self.scmap=scale
      self.projc1=scale*self.ERADIUS
      self.projc2=math.cos(math.radians(self.trulat1))
      if(self.projc2 <= 0.) :
        print '  Error (1) in SETMAPR for Lat,Lon Projection',
        '  Illegal true latitude provided: %f.' % self.trulat1
        sys.exit(1)

      self.projc3=self.projc1*self.projc2/d2rad
      print '  Map projection set to Lat, Lon',
      '  X origin, Y origin set to 0.,0. at the equator, 0. long',
      '  Y positive toward the North Pole.'

    elif( self.jproj == 5 ) :
      ##
      ##-----------------------------------------------------------------------
      ##
      ##  WDT mapproj
      ##
      ##  Approximate flat earth projection (using approximate great circle to
      ##  compute distances).
      ##
      ##                For a 512 km box
      ##                at a latitude of 75 it is off by 0.8 km
      ##                                 70              0.4 km
      ##                                 65              0.3 km
      ##                                 55              0.1 km
      ##                                 45              0.05 km
      ##                                 35              0.02 km
      ##                                  0              0.01 km
      ##                at the corners of the box.
      ##
      ##  For this projection:
      ##      self.projc1 = lat0
      ##      self.projc2 = COS(radians(lat0)
      ##      self.projc3 = lon0
      ##      self.projc4 = radians(scale*self.ERADIUS ## deg_to_km
      ##
      ##-----------------------------------------------------------------------
      ##
      self.trulat1=lat1
      self.rota=orient
      self.scmap=scale
      self.projc1 = self.trulat1
      self.projc2 = math.cos(math.radians(self.trulat1))
      self.projc3 = orient
      self.projc4 = math.radians(scale*self.ERADIUS) ## deg_to_km
      if(self.trulat1 <= 0.) :
        print '  Error (1) in SETMAPR for Lat,Lon Projection',
        '  Illegal true latitude provided: %f.' % self.trulat1
        sys.exit(1)
    else :
      print ' projection (%d) is not supported.' %iproj
      sys.exit(2)

  #enddef __init__

  ######################################################################

  def setorig(self,iopt,x0,y0) :
    '''
    Set the origin for the map projection.
    This is call after subroutine mapproj if the origin
    must be moved from the original position, which is the
    pole for the polar stereographic projection and the
    Lambert conformal, and the equator for Mercator.

    iopt        origin setting option
                1: origin given in corrdinate x,y
                2: origin given in lat,lon on earth

    x0          first coordinate of origin
    y0          second coordinate of origin

    '''
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    if( iopt == 1 ) :
      ##-----------------------------------------------------------------------
      ##
      ##  iopt=1 origin is given in x,y in absolute coordinates.
      ##
      ##-----------------------------------------------------------------------
      ##
      self.xorig=x0
      self.yorig=y0
      (rlat,rlon) = self.xytoll(0.,0.)

      print '  Coordinate origin set to absolute x = %f, y = %f,\n  Latitude = %f, longitude= %f.' % \
            (self.xorig,self.yorig,rlat,rlon)

    elif( iopt == 2 ) :

      ##
      ##-----------------------------------------------------------------------
      ##
      ##  iopt=2 origin is given in lat,lon on earth
      ##
      ##-----------------------------------------------------------------------
      ##
      ##
      self.xorig=0.
      self.yorig=0.
      (xnew,ynew) = self.lltoxy(x0,y0)
      self.xorig=xnew
      self.yorig=ynew
      print '  Coordinate origin set to absolute x = %f, y = %f,\n  Latitude = %f, longitude= %f.' % \
            (xnew,ynew,x0,y0)

    else :
      (rlat,rlon) = self.xytoll(0.,0.)
      print ' Setorig option %d not supported.\n' % iopt,     \
            '    Coordinate origin unchanged at x = %f, y = %f,\n   Latitude = %f, longitude= %f.' % \
            (self.xorig, self.yorig,rlat,rlon)

  #enddef setorig

  ####################################################################

  def xytoll(self,x,y) :
    '''
    Determine latitude and longitude given X,Y coordinates on
    map projection.

    INPUT:

      x        x in map coordinates
      y        y in map coordinates
               Units are meters unless the scale parameter is
               not equal to 1.0

    OUTPUT:

      rlat     latitude.
               (degrees, negative south, positive north)
      rlon     longitude.
               (degrees, negative west, positive east)
    '''

    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    if ( self.jproj == 0 ) :
      ratio= math.degrees(1/self.ERADIUS)
      rlat = ratio*(y+self.yorig)
      rlon = ratio*(x+self.xorig)

    elif( self.jproj == 1 ) :
      ##
      ##----------------------------------------------------------------
      ##
      ##  Polar Stereographic projection
      ##  For this projection:
      ##      self.projc1 is the scaled earth's radius, scale times self.ERADIUS
      ##      self.projc2 is the numerator of emfact, the map image scale factor.
      ##      self.projc3 is self.projc2 times the scaled earth's radius.
      ##
      ##----------------------------------------------------------------
      ##
      yabs=y+self.yorig
      xabs=x+self.xorig
      radius=sqrt( xabs*xabs + yabs*yabs )/self.projc3
      rlat = self.jpole*(90. - 2.*math.degrees(math.atan(radius)))
      rlat = min(rlat, 90. )
      rlat = max(rlat,-90. )

      if((self.jpole*yabs) > 0.) :
        dlon=180. + math.degrees(math.atan(-xabs/yabs))
      elif((self.jpole*yabs) < 0.) :
        dlon=math.degrees(math.atan(-xabs/yabs))
      elif (xabs > 0.) :     ## y=0.
        dlon=90.
      else :
        dlon=-90.

      rlon = self.rota + self.jpole*dlon
      if(rlon > 180)  : rlon=rlon-360.
      if(rlon < -180) : rlon=rlon+360.
      rlon=min(rlon, 180.)
      rlon=max(rlon,-180.)

    elif ( self.jproj == 2 ) :

      ##
      ##----------------------------------------------------------------
      ##
      ##  Lambert Conformal Conic Projection.
      ##  For this projection:
      ##      self.projc1 is the scaled earth's radius, scale times self.ERADIUS/n
      ##      self.projc2 is cos of self.trulat1
      ##      self.projc3 is tan (45. - trulat/2) a const for local map scale
      ##      self.projc4 is the cone constant, n
      ##
      ##----------------------------------------------------------------
      ##
      yabs=y+self.yorig
      xabs=x+self.xorig
      radius=sqrt( xabs*xabs+ yabs*yabs )
      ratio=self.projc3*((radius/(self.projc1*self.projc2))**(1./self.projc4))
      rlat=self.jpole*(90. -2.*math.degrees((math.atan(ratio))))
      rlat=min(rlat, 90.)
      rlat=max(rlat,-90.)

      yjp=self.jpole*yabs
      if(yjp > 0.) :
        dlon=180. + math.degrees(math.atan(-xabs/yabs))/self.projc4
      elif(yjp < 0.) :
        dlon=math.degrees(math.atan(-xabs/yabs))/self.projc4
      elif (xabs > 0.) :     ## y=0.
        dlon=90./self.projc4
      else :
        dlon=-90./self.projc4

      rlon= self.rota + self.jpole*dlon
      if(rlon > 180)  : rlon=rlon-360.
      if(rlon < -180) : rlon=rlon+360.
      rlon=min(rlon, 180.)
      rlon=max(rlon,-180.)

    elif( self.jproj == 3 ) :

      ##
      ##----------------------------------------------------------------
      ##
      ##  Mercator Projection.
      ##  For this projection:
      ##      self.projc1 is the scaled earth's radius, scale times self.ERADIUS
      ##      self.projc2 is cos of self.trulat1
      ##      self.projc3 is self.projc1 times self.projc2
      ##
      ##----------------------------------------------------------------
      ##
      yabs=y+self.yorig
      xabs=x+self.xorig
      rlat=(90. - 2.*math.degrees(math.atan(exp(-yabs/self.projc3))))
      rlat=min(rlat, 90.)
      rlat=max(rlat,-90.)
      dlon=math.degrees((xabs/self.projc3))
      rlon=self.rota + dlon
      if(rlon > 180)  : rlon=rlon-360.
      if(rlon < -180) : rlon=rlon+360.

    elif( self.jproj == 4 ) :

      ##
      ##----------------------------------------------------------------
      ##
      ##  Lat, Lon Projection.
      ##  For this projection:
      ##      self.projc1 is the scaled earth's radius, scale times self.ERADIUS
      ##      self.projc2 is cos of self.trulat1
      ##      self.projc3 is self.projc1 times self.projc2 times 180/pi
      ##
      ##----------------------------------------------------------------
      ##
      rlon=x+self.xorig
      rlat=y+self.yorig

    elif( self.jproj == 5 ) :

      ##
      ##----------------------------------------------------------------
      ##
      ##  WDT mapproj
      ##
      ##  Approximate flat earth projection (using approximate great circle to
      ##  compute distances).
      ##
      ##  For this projection:
      ##      self.projc1 = lat0
      ##      self.projc2 = COS(radians(lat0)
      ##      self.projc3 = lon0
      ##      self.projc4 = radians(scale*self.ERADIUS ## deg_to_km
      ##
      ##----------------------------------------------------------------
      ##
      yabs = y + self.yorig
      xabs = x + self.xorig
      rlat = self.projc1 + yabs/self.projc4
      rlat = min(rlat, 90.)
      rlat = max(rlat,-90.)
      rlon = xabs/self.projc4/(0.5*(self.projc2+math.cos(rlat*d2rad))) + self.projc3
      if(rlon > 180)  : rlon=rlon-360.
      if(rlon < -180) : rlon=rlon+360.
    else :
      print ' projection %d is not supported.' % self.jproj
      sys.exit(1)

    return (rlat,rlon)

  #enddef xytoll

  ####################################################################

  def lltoxy(self,rlat,rlon) :

    '''
    Determine x, y coordinates on map projection from the given latitude
    and longitude.

    INPUT:

      rlat     latitude.
               (degrees, negative south, positive north)

      rlon     longitude.
               (degrees, negative west, positive east)

    OUTPUT:

      xloc     x in map coordinates
      yloc     y in map coordinates

    '''

    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    if( self.jproj == 0 ) :
      ratio=math.radians(self.ERADIUS)
      xloc = ratio*rlon - self.xorig
      yloc = ratio*rlat - self.yorig

    elif( self.jproj == 1 ) :
      ##
      ##----------------------------------------------------------------
      ##
      ##  Polar Stereographic projection
      ##  For this projection:
      ##      self.projc1 is the scaled earth's radius, scale times self.ERADIUS
      ##      self.projc2 is the numerator of emfact, the map image scale factor.
      ##      self.projc3 is self.projc2 times the scaled earth's radius.
      ##
      ##----------------------------------------------------------------
      ##
      denom=(1. + math.sin(math.radians(self.jpole*rlat)))
      if(denom == 0.) : denom=1.0E-10
      radius=self.jpole*self.projc3*math.cos(math.radians(rlat))/denom
      dlon=self.jpole*math.radians((rlon-self.rota))
      xloc= radius*math.sin(dlon) - self.xorig
      yloc=-radius*math.cos(dlon) - self.yorig

    elif( self.jproj == 2 ) :

      ##
      ##----------------------------------------------------------------
      ##
      ##  Lambert Conformal Conic Projection.
      ##  For this projection:
      ##      self.projc1 is the scaled earth's radius, scale times self.ERADIUS/n
      ##      self.projc2 is cos of self.trulat1
      ##      self.projc3 is tan (45. - trulat/2) a const for local map scale
      ##      self.projc4 is the cone constant, n
      ##
      ##----------------------------------------------------------------
      ##

      ## Handle opposite pole
      if (self.jpole*rlat < -89.9) :
        lat = -89.9 * self.jpole
      else :
        lat = rlat

      radius=self.projc1*self.projc2*(math.tan(math.radians((45.-0.5*self.jpole*lat)))/self.projc3)**self.projc4
      tem = rlon-self.rota
      if( tem < -180.0) : tem = 360.0+tem
      if( tem > 180.0)  : tem = tem-360.0
      dlon=self.projc4*math.radians(tem)
      xloc=            radius*math.sin(dlon) - self.xorig
      yloc=-self.jpole*radius*math.cos(dlon) - self.yorig

    elif( self.jproj == 3) :
      ##
      ##----------------------------------------------------------------
      ##
      ##  Mercator Projection.
      ##  For this projection:
      ##      self.projc1 is the scaled earth's radius, scale times self.ERADIUS
      ##      self.projc2 is cos of self.trulat1
      ##      self.projc3 is self.projc1 times self.projc2
      ##
      ##----------------------------------------------------------------
      ##
      dlon=rlon-self.rota
      if(dlon < -180.) : dlon=dlon+360.
      if(dlon > 180.)  : dlon=dlon-360.
      xloc=self.projc3*math.radians(dlon) - self.xorig
      denom=math.tan(math.radians((45. - 0.5*rlat)))
      if( denom <= 0. ) : denom=1.0E-10
      yloc=-self.projc3*log(denom) - self.yorig

    elif(self.jproj == 4) :

      ##
      ##----------------------------------------------------------------
      ##
      ##  Lat, Lon Projection.
      ##  For this projection:
      ##      self.projc1 is the scaled earth's radius, scale times self.ERADIUS
      ##      self.projc2 is cos of self.trulat1
      ##      self.projc3 is self.projc1 times self.projc2 times 180/pi
      ##
      ##----------------------------------------------------------------
      ##
      xloc=rlon-self.xorig
      if (xloc < -180.) : xloc=xloc+360.
      if (xloc >  180.) : xloc=xloc-360.
      yloc=rlat-self.yorig

    elif( self.jproj == 5 ) :
      ##
      ##----------------------------------------------------------------
      ##
      ##  WDT mapproj
      ##
      ##  Approximate flat earth projection (using approximate great circle to
      ##  compute distances).
      ##
      ##  For this projection:
      ##      self.projc1 = lat0
      ##      self.projc2 = cos(radians(lat0)
      ##      self.projc3 = lon0
      ##      self.projc4 = radians(scale*self.ERADIUS ## deg_to_km
      ##
      ##----------------------------------------------------------------
      ##
      xloc = self.projc4*(rlon - self.projc3)* 0.5*(self.projc2 + math.cos(rlat*d2rad)) - self.xorig
      yloc = self.projc4*(rlat - self.projc1) - self.yorig
    else :
      print ' projection %d is not supported.' % (self.jproj)
      sys.exit(1)

    return (xloc, yloc)

  #enddef lltoxy

  def getPROJstr(self) :
    ''' Get ARPS map project in PROJ4 projection string.'''

    if (self.jproj == 0) :
      projstr = '"+porj=latlong +ellps=sphere"'
    elif (self.jproj == 1):    ## Polar stereographic (pure north or south)
      projstr = '"+proj=stere +lat_0=%.2f +lat_ts=%.2f +lon0=%.2f +ellps=sphere"' % \
                (self.jpole*90.0,self.trulat1,self.rota)
    elif (self.jproj == 2):    ## Lambert
      (ctrlat,ctrlon) = self.xytoll(0.,0.)
      projstr='"+proj=lcc +lat_1=%.2f +lat_2=%.2f +lat_0=%.2f +lon_0=%.2f +ellps=sphere"' % \
              (self.trulat1,self.trulat2,ctrlat,ctrlon)
    elif (self.jproj == 3):    ## Mercator
      projstr='"+proj=ob_tran +lat_ts=%.2f +lon_0=%.2f +ellps=sphere"' % \
              (self.trulat1,self.rota)
    ##elif (self.jproj == 6):  ## cassini or rotated lat/lon
    ##  projstr='"+proj=merc +o_proj=latlong +o_lat_p=%.2f +o_lon_p=%.2f +lon_0=%.2f +ellps=sphere"' % \
    ##          (self.ctrlat,self.ctrlon,self.rota)
    else:
      print "Unsupported map projection %d." % self.jproj
      sys.exit(-1)

    return projstr

  #enddef

#endclass mapprojection

#################### class domainrange ###############################

class domainrange(object) :

  def __init__(self,arpsmap,xs,xe,ys,ye,zs,ze,dx,dy) :
    self.xyzmass = (xs,ys,zs,xe,ye,ze)
    (maxlat,rlon) = arpsmap.xytoll(0.0, ye)
    (rlat,minlon) = arpsmap.xytoll(xs,  ye)
    (rlat,maxlon) = arpsmap.xytoll(xe,  ye)

    (minlat,rlon) = arpsmap.xytoll(xs,  ys)
    (rlat,  rlon) = arpsmap.xytoll(xe,  ys)

    minlat = min(minlat,rlat)
    self.llmass = (minlon,minlat,maxlon,maxlat)

    xs -= dx/2.0
    xe += dx/2.0
    ys -= dy/2.0
    ye += dy/2.0

    (maxlat,rlon) = arpsmap.xytoll(0.0, ye)
    (rlat,minlon) = arpsmap.xytoll(xs,  ye)
    (rlat,maxlon) = arpsmap.xytoll(xe,  ye)

    (minlat,rlon) = arpsmap.xytoll(xs,  ys)
    (rlat,  rlon) = arpsmap.xytoll(xe,  ys)

    minlat = min(minlat,rlat)
    self.llvector = (minlon,minlat,maxlon,maxlat)

  #enddef __init__

  def printllmass(self) :
    print '  Box lat/lon range are: %.2f %.2f %.2f %.2f' % self.llmass
  #enddef

  def printllvector(self) :
    print '  Vector Box lat/lon range are: %.2f %.2f %.2f %.2f' % self.llvector
  #enddef

  def printxyzmass(self) :
    print '''  Box extends are: (x0:y0:z0) = %.2f:%.2f:%.2f
                   (x1:y1:z1) =  %.2f: %.2f: %.2f ''' % self.xyzmass
  #enddef

#endclass domainrange

class modeltime(object) :

  def __init__(self,timestr) :
    ##timestr = '%04d-%02d-%02d_%02d:%02d:%02d' % (iyr, imon, iday, ihr, imin, isec)
    mtime = time.strptime(timestr+' UTC',"%Y-%m-%d_%H:%M:%S %Z")    ## struct_time
    esec = time.mktime(mtime)
    self.abstsec = esec - time.timezone                ## converted to UTC time
  #enddef __init__

  def getstr(self,itime=0) :

    newtime = time.gmtime(self.abstsec + itime)
    timestr = time.strftime('%Y-%m-%d_%H:%M:%S',newtime)

    return timestr
  #enddef getstr

  def getsec(self,itime=0) :

    newtime = self.abstsec + itime

    return newtime
  #enddef getsec

#endclass modeltime

############################# testNetcdf ###############################
def testNetcdf() :
  ## Open the file
  file = NetCDFFile('testFile.cdf', 'w')

  ## Create some global_ attribute using a constant
  setattr(file, 'versionNumber', 1)

  ## Create a global_ attribute using a variable
  magicNum = 42
  setattr(file, 'magicNumber', magicNum)

  ## Get the value of these global_ variables
  val1 = getattr(file, 'versionNumber')
  val2 = getattr(file, 'magicNumber')
  print "versionNumber =", val1, "magicNumber =", val2

  ## Make some dimensions
  file.createDimension('smallDim', 4)
  file.createDimension('mediumDim', 25)
  file.createDimension('largeDim', 100)

  ## Make a new variable
  varDims = ('smallDim', 'mediumDim')
  var1 = file.createVariable('varOne', 'i', varDims)

  ## Get the size of a variable
  var1Shape = var1.shape
  print "The size of varOne is:", var1Shape

  ## Put some data in the variable
  for i in arange(var1Shape[0]):
      for j in arange(var1Shape[1]):
          var1[i, j] = i * j

  #Make another variable
  varDims = ('mediumDim', 'largeDim')
  file.createVariable('varTwo', 'f', ('mediumDim', 'largeDim'))

  ## Get the new variable (could have done this when we created var2)
  var2 = file.variables['varTwo']

  ## Initialize the new variable to 0.0
  for i in arange(var2.shape[0]):
      for j in arange(var2.shape[1]):
          var2[i, j] = 0.0

  ## Get the data from first variable
  data = var1.getValue()

  ## Print out the data
  print data

  ## Get the dimension names of var1
  dimNames = var1.dimensions
  print "Dimension names of var1:", dimNames

  ## Create some attributes for var1
  setattr(var1, 'units', 'Degrees C')
  setattr(var1, 'precision', 2)
  setattr(var1, 'maxValue', 19.999)

  ## Read the variable attributes we just created
  att1 = getattr(var1, 'units')
  att2 = getattr(var1, 'precision')
  att3 = getattr(var1, 'maxValue')

  ## Print the value of these variable attributes
  print "units =", att1, 'precision = ', att2, 'maxValue = ', att3

  ## Close the netCDF file
  file.close()
#enddef testNetcdf

#############################  Portral   ###############################
if __name__ == "__main__":

  _debug = False

  cmd  = sys.argv[0]
  argsdict = parseArgv(sys.argv[1:])

  main(**argsdict)
