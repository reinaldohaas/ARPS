#!/usr/bin/python
## Script file for ncrad2arps to generate radar observations
## in a format that can be used in the EnKF analysis.
##
## Youngsun Jung                                 08/24/2009
##
## Import utilities
import sys,string,time,os,dircache,calendar

## Set parameters (User define variables)
datapath = '.'
workpath = '.'
binpath = '/data3/yjung/arps5.3_EnKF_v1.1/bin'
reftime_start = "20040530-000500"    # First analysis time
reftime_2nd   = "20040530-001000"    # Second analysis time
reftime_end   = "20040530-010000"    # The last analysis time

############################################
## Start excution
############################################
os.chdir(datapath)

## Initialization of arrays
fnstr = []; dtstr = []; dirlist =[]
list1_fn = []; list2_fn = []
dt_1 = []; dt_2 = []; dt_3 = []

for i in range(0,15):
  dt_1.append(10000)
  dt_2.append(10000)
  dt_3.append(10000)
for i in range(0,71):
  list1_fn.append("")
  list2_fn.append("")

## Get reference time and time loop starting and ending time
timeref = time.strptime(reftime_start,"%Y%m%d-%H%M%S")
time_s = calendar.timegm(timeref)
print "Start time: ", time_s
time_temp = time.strptime(reftime_2nd,"%Y%m%d-%H%M%S")
time_interval = calendar.timegm(time_temp)
time_interval = time_interval - time_s
print "Time interval: ",time_interval
time_temp = time.strptime(reftime_end,"%Y%m%d-%H%M%S")
time_e = calendar.timegm(time_temp)+time_interval
print "End   time: ", time_e

## Get file names in the data directory
def lister(dirlist,dirname,files):
    dirlist.append(os.path.abspath(dirname))

os.path.walk('.',lister,dirlist)

##print len(dirlist)

for dirPath in dirlist:
  a = dircache.listdir(dirPath)
  a = a[:]

## Now move to the working directory
os.chdir(workpath)

## time loop 
for time_r in range(time_s,time_e,time_interval):
  print "Reference time:", time_r
  time_temp = time.gmtime(time_r)
  ncrad_time = "%04d%02d%02d-%02d%02d"%(time_temp.tm_year,\
     time_temp.tm_mon,time_temp.tm_mday,time_temp.tm_hour,time_temp.tm_min)
  dt_1[:] = []; dt_2[:] = []; dt_3[:] = []
  list1_fn[:] = []; list2_fn[:] = []

  for i in range(0,15):
    dt_1.append(10000)
    dt_2.append(10000)
    dt_3.append(10000)
  for i in range(0,71):
    list1_fn.append("")
    list2_fn.append("")

  for fileName in a:
    if (fileName[0:3] == "Ref"):
      fnstr[:] = []
      dtstr[:] = []
##      print "     found ",fileName
      fnstr = string.split(fileName,"_")

      ## Elevation
      if(fnstr[1] == "0.00"):
         ilev = 0
      elif (fnstr[1] == "0.50"):
         ilev = 1
      elif (fnstr[1] == "1.50"):
         ilev = 2
      elif (fnstr[1] == "2.50"):
         ilev = 3
      elif (fnstr[1] == "3.50"):
         ilev = 4
      elif (fnstr[1] == "4.50"):
         ilev = 5
      elif (fnstr[1] == "5.50"):
         ilev = 6
      elif (fnstr[1] == "6.50"):
         ilev = 7
      elif (fnstr[1] == "7.50"):
         ilev = 8
      elif (fnstr[1] == "8.70"):
         ilev = 9
      elif (fnstr[1] == "10.00"):
         ilev = 10
      elif (fnstr[1] == "12.00"):
         ilev = 11
      elif (fnstr[1] == "14.00"):
         ilev = 12
      elif (fnstr[1] == "16.70"):
         ilev = 13
      elif (fnstr[1] == "19.50"):
         ilev = 14

      ## Get two closest time data to the reference time
      dtstr = string.split(fnstr[2],".")
      timebf = time.strptime(dtstr[0],"%Y%m%d-%H%M%S")
      time_p = calendar.timegm(timebf)
##      print "     Data time: ", time_p

      dt_1[ilev] = time_r - time_p
      if (dt_1[ilev] > 0 and dt_1[ilev] < dt_2[ilev]):
        dt_2[ilev] = dt_1[ilev]
        list1_fn[ilev] = datapath+'/'+fileName+'\n'

      elif (dt_1[ilev] < 0 and -dt_1[ilev] < dt_3[ilev]):
        dt_3[ilev] = -dt_1[ilev]
        list2_fn[ilev] = datapath+'/'+fileName+'\n'

  ## file lists for Vr, Zdr, Kdp, and RhoHV
  for i in range(1,15):
    list1_fn[i+14] = string.replace(list1_fn[i],"Reflectivity","Velocity")
    list1_fn[i+28] = string.replace(list1_fn[i],"Reflectivity","Zdr")
    list1_fn[i+42] = string.replace(list1_fn[i],"Reflectivity","Kdp")
    list1_fn[i+56] = string.replace(list1_fn[i],"Reflectivity","RhoHV")
    list2_fn[i+14] = string.replace(list2_fn[i],"Reflectivity","Velocity")
    list2_fn[i+28] = string.replace(list2_fn[i],"Reflectivity","Zdr")
    list2_fn[i+42] = string.replace(list2_fn[i],"Reflectivity","Kdp")
    list2_fn[i+56] = string.replace(list2_fn[i],"Reflectivity","RhoHV")

  ## Creat listl and list2 which contain lists of input file names
  fnlist1 = open('list1','w')
  for i in range(1,71):
    fnlist1.write(list1_fn[i])
##      print list1_fn[i]
  fnlist1.close()

  fnlist2 = open('list2','w')
  for i in range(1,71):
    fnlist2.write(list2_fn[i])
  fnlist2.close()
##      print list2_fn[i]

  ## Run ncrad2arps
  ncrad2arps = binpath+'/ncrad2arps KOUN list1 list2 -gridtilt -reftime %s < arps.input'\
               %(ncrad_time)
  print ncrad2arps
  os.system(ncrad2arps)
