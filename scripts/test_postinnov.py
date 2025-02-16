#!/usr/bin/python
## This is user script program for postinnov, 
## which computes innovation statistics.
##
## Youngsun Jung (8/13/2009)
##
## Env for sooner
#BSUB -q normal
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -x
#BSUB -W 02:00
#BSUB -J "test"
#BSUB -u youngsun.jung@ou.edu

## Import utilities
import sys, string, time, os, shutil, stat, pty, dircache
from operator import *

## Set parameters (User defined variables)
nen = 60
tsta = int(60)            ## initial time to compute innovation
tinter = int(60)          ## time interval btw successive analaysis
cycle = int(5)            ## number of assimilation cycle
rmsfcst = int(2)          ## See the same option in arpsenkf.input

## set absolute path
workpath = '/scratch/yjung/EnKF/REALCASE/exp4'
binpath = "/home/yjung/arps5.3/bin"   ## binary dir

############################################
## Start excution
############################################
print time.strftime('%X %x %Z')

os.chdir(workpath)

## Purge old links
M = "rm -rf %s/*.hdfenkfi %s/*.hdfenkfo"%(workpath,workpath)
os.system(M)

## repeat computation within DA windows
for k in range (1,cycle+1):

  time1 = tsta + (k - 1) * tinter

  for i in range (1,nen+1):
    numen = '%03d'%i
    if(rmsfcst == 1 or rmsfcst == 2):
      file1 = '%s/enf%s.hdf%06d'%(workpath,numen,time1)
      file2 = '%s/enf%s.hdfenkfi'%(workpath,numen)
      while(not(os.path.exists(file1))):
        sys.exit("No enf file exists to link")
      ln = 'ln -s %s %s'%(file1,file2)
      os.system(ln)

    if(rmsfcst == 0 or rmsfcst == 2):
      file1 = '%s/ena%s.hdf%06d'%(workpath,numen,time1)
      file2 = '%s/enf%s.hdfenkfo'%(workpath,numen)
      while(not(os.path.exists(file1))):
        sys.exit("No enf file exists to link")
      ln = 'ln -s %s %s'%(file1,file2)
      os.system(ln)

  M = '%s/postinnov %s/arps.input < %s/arpsenkf.input'%(binpath,workpath,workpath)
  os.system(M)

  M = "rm -rf %s/*.hdfenkfi %s/*.hdfenkfo"%(workpath,workpath)
  os.system(M)

print time.strftime('%X %x %Z')
print "end of program"
