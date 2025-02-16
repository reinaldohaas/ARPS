#
#     ##################################################################
#     ##################################################################
#     ######                                                      ######
#     ######      Advanced Regional Prediction System (ARPS)      ######
#     ######                   Version 5.2                        ######
#     ######                                                      ######
#     ######                     Developed by                     ######
#     ######     Center for Analysis and Prediction of Storms     ######
#     ######                University of Oklahoma                ######
#     ######                                                      ######
#     ##################################################################
#     ##################################################################
#
#=======================================================================
#
#  PURPOSE: This makefile generates the ARPS executables. It can be run
#           directly or through the makearps script.
#
#  AUTHOR:  Yuhe Liu
#
#  08/15/1997 (Yuhe Liu)
#  A complete new makefile in responding to the re-structure of ARPS
#  source code.
#
#  13 March 2002 (Eric Kemp)
#  Added WRF BMJ source code files.
#
#  05/01/2002 (Jason Levit)
#  Added the ARPS verification source code files.
#
#  6 June 2002 (Eric Kemp)
#  Added intrpsoil3d.f90.
#
#  OTHER INFORMATION:
#       See the makearps command.
#
#=======================================================================

#
#-----------------------------------------------------------------------
#
# Version number
#
#-----------------------------------------------------------------------
#

VERSION = 5.2.0

#
#-----------------------------------------------------------------------
#
# Files to be included
#
#-----------------------------------------------------------------------
#

INCFILE = /dev/null

include $(INCFILE)

GEMINC = GEMPRM.AIX

MACHDEPLIB = GENLIBOBJ
ASSIMOPT = ASSIMDUMMY
NCARGDEP = NCARGOBJ_OFF
ZXPLOTDEP = ZXPLOTOBJ_POST

MPIDEP  = MPIOBJ_OFF

HDFDEP  = HDFOBJ_OFF
PAKDEP  = PAKOBJ_OFF
SVIDEP  = SVIOBJ_OFF
NETDEP  = NETOBJ_OFF
V5DDEP  = V5DOBJ_OFF
GEMDEP  = GEMOBJ_OFF
PHDF5DEP = PHDF5OBJ_OFF
GRIB2DEP = GRIB2_OFF
CRTMDEP = CRTM_OFF
CITMDEP = CITM_OFF

#
#-----------------------------------------------------------------------
#
# Default shell
#
#-----------------------------------------------------------------------
#

# SHELL=/bin/csh

#
#-----------------------------------------------------------------------
#
# Commands to be used
#
#-----------------------------------------------------------------------
#

FTN = f90
LDR = f90
CC  = cc
CPP = cpp

AWK = awk
CP  = cp
LN  = ln
RM  = rm
AR  = ar
TAR = tar
ZIP = gzip
ARFLAG =

ARPS_LD = $(LDR)

PHDF5PATH =

TOPDIR = `pwd`
WRKDIR = $(TOPDIR)

BINDIR = $(WRKDIR)/bin
LIBDIR = $(WRKDIR)/lib
MODDIR = $(WRKDIR)/modules

#-----------------------------------------------------------------------
#
# Directories
#
#-----------------------------------------------------------------------

INCL_DIR   = include
DATA_DIR   = data
DOC_DIR    = docs
INPUT_DIR  = input
SCRPT_DIR  = scripts
SND_DIR    = sounding
SRC_DIR    = src

ATAB_DIR   = $(DATA_DIR)/adas
PLTD_DIR   = $(DATA_DIR)/arpsplt
VERD_DIR   = $(DATA_DIR)/arpsverif
PLTR_DIR   = $(DATA_DIR)/pltradcol

ARPS_DIR   = $(SRC_DIR)/arps
ADAS_DIR   = $(SRC_DIR)/adas
3DVAR_DIR  = $(SRC_DIR)/arps3dvar
ENKF_DIR   = $(SRC_DIR)/arpsenkf
AGR_DIR    = $(SRC_DIR)/arpsagr
ASSIM_DIR  = $(SRC_DIR)/arpsassim
CVT_DIR    = $(SRC_DIR)/arpscvt
CVTOBS_DIR  = $(SRC_DIR)/arpscvtobs
DIF_DIR    = $(SRC_DIR)/arpsdiff
ENS_DIR    = $(SRC_DIR)/arpsens
E2A_DIR    = $(SRC_DIR)/ext2arps
A2W_DIR    = $(SRC_DIR)/arps2wrf
A4W_DIR    = $(SRC_DIR)/arps4wrf
A2C_DIR    = $(SRC_DIR)/arps2coamps
C2A_DIR    = $(SRC_DIR)/coamps2arps
N2A_DIR    = $(SRC_DIR)/nmm2arps
W2A_DIR    = $(SRC_DIR)/wrf2arps
WRF_EXT_DIR   = $(SRC_DIR)/external
WRF_PHDF5_DIR = $(SRC_DIR)/external/io_phdf5
WRF_INT_DIR   = $(SRC_DIR)/external/io_int
WRF_NET_DIR   = $(SRC_DIR)/external/io_netcdf
G2_EXT_DIR = $(SRC_DIR)/external/g2lib
AVN_DIR    = $(SRC_DIR)/ext2arps
EXTSND_DIR = $(SRC_DIR)/arpsextsnd
SKEWT_DIR  = $(SRC_DIR)/skewt
INTRP_DIR  = $(SRC_DIR)/arpsintrp
TINTRP_DIR  = $(SRC_DIR)/arpstintrp
TRAJC_DIR   = $(SRC_DIR)/arpstrajc
WTRETCOL_DIR = $(SRC_DIR)/wtretcol
MP_DIR     = $(SRC_DIR)/arps_mp
PLT_DIR    = $(SRC_DIR)/arpsplt
PRT_DIR    = $(SRC_DIR)/arpsprt
SFC_DIR    = $(SRC_DIR)/arpssfc
SOIL_DIR   = $(SRC_DIR)/arpssoil
TERN_DIR   = $(SRC_DIR)/arpstern
TRN_DIR    = $(SRC_DIR)/arpstrn
88D_DIR    = $(SRC_DIR)/88d2arps
RDR_EMUL   = $(SRC_DIR)/radaremul
MCI_DIR    = $(SRC_DIR)/mci2arps
MOSAIC_DIR = $(SRC_DIR)/zmosaic2arps
MLSO_DIR   = $(SRC_DIR)/meso2lso
VERIF_DIR  = $(SRC_DIR)/arpsverif
IPLIB_DIR  = $(VERIF_DIR)/iplib
ZXPLOT_DIR = $(SRC_DIR)/zxplot

INCLDIR  = $(TOPDIR)/$(INCL_DIR)
DATADIR  = $(TOPDIR)/$(DATA_DIR)
DOCDIR   = $(TOPDIR)/$(DOC_DIR)
INPUTDIR = $(TOPDIR)/$(INPUT_DIR)
SCRPTDIR = $(TOPDIR)/$(SCRPT_DIR)
SNDDIR   = $(TOPDIR)/$(SND_DIR)
ATABDIR  = $(TOPDIR)/$(ATAB_DIR)
VERDDAT  = $(TOPDIR)/$(VERF_DIR)

ADASDIR  = $(TOPDIR)/$(ADAS_DIR)
3DVARDIR = $(TOPDIR)/$(3DVAR_DIR)
ENKFDIR  = $(TOPDIR)/$(ENKF_DIR)
AGRDIR   = $(TOPDIR)/$(AGR_DIR)
ASSIMDIR = $(TOPDIR)/$(ASSIM_DIR)
ARPSDIR  = $(TOPDIR)/$(ARPS_DIR)
CVTDIR   = $(TOPDIR)/$(CVT_DIR)
CVTOBSDIR = $(TOPDIR)/$(CVTOBS_DIR)
BUDGDIR  = $(TOPDIR)/$(BUDG_DIR)
DIFDIR   = $(TOPDIR)/$(DIF_DIR)
ENSDIR   = $(TOPDIR)/$(ENS_DIR)
E2ADIR   = $(TOPDIR)/$(E2A_DIR)
A2WDIR   = $(TOPDIR)/$(A2W_DIR)
A4WDIR   = $(TOPDIR)/$(A4W_DIR)
A2CDIR   = $(TOPDIR)/$(A2C_DIR)
C2ADIR   = $(TOPDIR)/$(C2A_DIR)
N2ADIR   = $(TOPDIR)/$(N2A_DIR)
W2ADIR   = $(TOPDIR)/$(W2A_DIR)
PHDF5DIR = $(TOPDIR)/$(WRF_PHDF5_DIR)
WRFINTDIR= $(TOPDIR)/$(WRF_INT_DIR)
WRFNETDIR= $(TOPDIR)/$(WRF_NET_DIR)
G2DIR    = $(TOPDIR)/$(G2_EXT_DIR)
AVNDIR   = $(TOPDIR)/$(AVN_DIR)
EXTSNDIR = $(TOPDIR)/$(EXTSND_DIR)
SKEWTDIR = $(TOPDIR)/$(SKEWT_DIR)
INTRPDIR = $(TOPDIR)/$(INTRP_DIR)
TINTRPDIR= $(TOPDIR)/$(TINTRP_DIR)
TRAJCDIR = $(TOPDIR)/$(TRAJC_DIR)
WTRETCOLDIR = $(TOPDIR)/$(WTRETCOL_DIR)
MPDIR    = $(TOPDIR)/$(MP_DIR)
PLTDIR   = $(TOPDIR)/$(PLT_DIR)
PRTDIR   = $(TOPDIR)/$(PRT_DIR)
SFCDIR   = $(TOPDIR)/$(SFC_DIR)
SOILDIR  = $(TOPDIR)/$(SOIL_DIR)
TERNDIR  = $(TOPDIR)/$(TERN_DIR)
TRNDIR   = $(TOPDIR)/$(TRN_DIR)
88DDIR   = $(TOPDIR)/$(88D_DIR)
RDREMULDIR  = $(TOPDIR)/$(RDR_EMUL)
MCIDIR   = $(TOPDIR)/$(MCI_DIR)
MOSAICDIR  = $(TOPDIR)/$(MOSAIC_DIR)
VERIFDIR = $(TOPDIR)/$(VERIF_DIR)
IPLIBDIR = $(TOPDIR)/$(IPLIB_DIR)
ZXPLOTDIR = $(TOPDIR)/$(ZXPLOT_DIR)
MLSODIR  = $(TOPDIR)/$(MLSO_DIR)

F2F90DIR = $(SRC_DIR)/f77tof90

#-----------------------------------------------------------------------
#
# Compiler Flag of Options
#
#-----------------------------------------------------------------------

CPPFLAGS = -C -P
FFLAGS =
FFLAGS_main =
FIXFLAGS =
FREEFLAGS =
CFLAGS =
LDFLAGS =
ICEFLAG =
RADFLAG =
THMFLAG =
MAKEOPT =
TAROPT  = -cf
TARADD  = -uf

#-----------------------------------------------------------------------
#
# Define HDF, Savi3D, NetCDF, PVM, and MPI libraries
#
#-----------------------------------------------------------------------

SVILIBS = -munch /usr/local/MeRAF/bin            \
          -start /usr/local/MeRAF/lib            \
          -L/usr/lib -L/usr/local/MeRAF/lib      \
          /usr/local/MeRAF/lib/libMERAF.a        \
          -lssC /usr/local/MeRAF/lib/libnetcdf.a \
          /usr/local/MeRAF/lib/libCns.a -lxlf

88D2A_RTLIBS = -la2rt -lnexrad -lnssl
88D2A_A2LIBS = /usr1/a2io/liba2tp.a \
               /usr1/tpio/libtpio.a

HDFLIBS = -L/usr/local/hdf/lib -ldf -lmfhdf -lz -ljpeg
NETLIBS = -lnetcdf

PVMLIBS = -L/opt/ctl/mpt/1.1.0.2/pvm3/lib/CRAY/libfpvm3 \
          -lfpvm3 -lpvm3 -lgpvm3

MPILIBS = -lmpi

LIBS = $(SVILIBS) $(HDFLIBS) $(NETLIBS) $(MPILIBS) \
       $(88D2A_RTLIBS) $(88D2A_A2LIBS)

#
#-----------------------------------------------------------------------
#
# Executable commands to be generated by this make file:
#
# ADASEXE      = adas          ADAS - ARPS Data Analysis System
# ADASMPEXE    = adas_mpi      MP version of ADAS
# 3DVAREXE     = arps3dvar     arps3DVAR- ARPS Data Analysis System
# ENKFEXE      = arpsenkf      arpsENKF- ARPS EnKF Data Assimilation System
# SHIFTEXE     = arpsshift     Phase correction analysis
# ARPSEXE      = arps          ARPS model executable
# ARPSMPEXE    = arps_mpi      MP version of ARPS model executable
# AGREXE       = arpsagr       AGR - Adaptive Grid Refinement
# ASSIMEXE     = arpsassim     Assimilation version of the ARPS
# ARPSCVTEXE   = arpscvt       History data dump format convesion exec.
# ARPSCVTOBSEXE  = arpscvtobs  Convert LSO surface observations to HDF
#                              format files.
# ARPSMPBUDGET = arpsmpbudget  Calculates microphysical budget statistics
# INTRPEXE     = arpsintrp     Program to interpolate between ARPS gridded data
# INTRP_LSEXE  = arpsintrp_ls  Program to interpolate large size ARPS data
# TINTRPEXE    = arpstintrp    Program to interpolate between ARPS gridded data
#                              at two times
# RADARINTRPEXE = radardtaintrp Program to interpolate radar data from one grid
#                               to another
# TRAJCEXE    = arpstrajc
# PLTNCAREXE   = arpspltncar   ZXPLOT based graphic analysis program exec.
#                              GKS metafile will be generated.
# PLTPOSTEXE   = arpspltpost   ZXPLOT based graphic analysis program exec.
#                              Postscript file will be generated.
# PLTMAXEXE    = arpspltmax    ZXPLOT based graphic analysis program exec.
# PLTGRIDEXE   = pltgrid       ZXPLOT based graphic analysis program exec.
# PRTEXE       = arpsprt       Formatted table printing prog. exec.
# DIFFEXE      = arpsdiff      Difference of ARPS history files
# ENSCVEXE     = arpsenscv     ARPS ensemble program
# ENSBCEXE     = arpsensbc     ARPS ensemble program
# ENSICEXE     = arpsensic     ARPS ensemble program
# ENSICMPEXE   = arpsensic_mpi MP version of arpsensic
# ENKFICEXE    = arpsenkfic    ARPS ENKF program
# RNDPRTEXE    = rndprt        ARPS ENKF program
# EXTSNDEXE    = extsnd        Extracts a sounding from ARPS history file
# EXTSNDMPEXE  = extsnd_mpi    MP version of "arpsextsnd"
# ARPSSFCEXE   = arpssfc       ARPS surface data pre-processor
# ARPSSOILEXE  = arpssoil      ARPS soil data pre-processor
# TERNEXE      = arpstern      ARPS terrain data pre-processor
# TRNEXE       = arpstrn       ARPS new terrain data pre-processor
# MERGETRNEXE  = mergetrn      ARPS terrain data merger/smoother
# DIR1DEGEXE   = dir1deg       Convertor of terrain data to direct
#                              access file for 1 degree resolution
# DIR5MINEXE   = dir5min       Convertor for 5 minutes resolution
# DIR30SECEXE  = dir30sec      Convertor for 30 seconds resolution
# E2AEXE       = ext2arps        Convertor of external data to ARPS grid
# E2AMPEXE     = ext2arps_mpi    MP version of ext2arps executable
# L2AEXE       = ext2arps.laps   Convertor of OLAPS data to ARPS grid
# G2AEXE       = ext2arps.gempak Convertor of GEMPAK data to ARPS grid
# DIFOBSEXE    = difobs        Difference of ARPS history file from ADAS obs
# 88D2ARPS     = 88d2arps      Remap WSR-88D radar data to ARPS grid
# 88D2ARPS_FAKE= 88d2arps_fake Remap fake radar data to ARPS grid
# NIDS2ARPS    = nids2arps     Remap NIDS format 88D radar data to ARPS grid
# RADMOSAIC    = radmosaic     Mosaic 88d2arps or nids2arps output for display
# ARPS2RAD     = arps2rad      Create simulated radar observations from ARPS data
# NCRAD2AEXE   = ncrad2arps    Remap NetCDF format 88D radar data to ARPS grid
# RADEMULEXE   = radaremul     Create simulated radar az-ran data from ARPS data
# RADB2CDFEXE  = radbin2cdf    Convert simulated radar data from binary to netcdf
# RADSECTEXE   = radsector     Convert simulated radar data from binary to netcdf
# ATTENEXE     = atten         Compute attenuation to tornado locations
# TMATRIXEXE   = tmatrix       Compute scattering amplitudes using T-matrix method
#
# MCI2ARPS     = mci2arps
# SAT2ARPS     = sat2arps
# MERGESAT     = mergesat
# PLTRADCOLEXE = pltradcol
# RADARPLTNCAREXE = radarpltncar
# RADARPLTPOSTEXE = radarpltpost
# PLTSATFLDEXE = pltsatfld
# WTRETCOLEXE  = wtretcol
# VERIFEXE     = arpsverif     ARPS verification package
# VERIFMPEXE   = arpsverif_mpi MP version of arpsverif executable
# A2GEXE       = arps2gem      Convert ARPS history file to GEMPAK format
# A2NEXE       = arps2ncdf     Convert ARPS history file to netCDF format
# A2EEXE       = arps2eta212   Convert ARPS history file to ETA #212 format
# A2WEXE       = arps2wrf      Convert ARPS history file to WRF format
# W2AEXE       = wrf2arps      Convert WRF output to ARPS format
# WEXTSNDEXT   = wrfextsnd     Extract sounding data from WRF file
#
# ARPSREADEXE  = arpsread      Sample program to read ARPS history file
# RAINDIFFEXE  = arpsraindiff
#
# INITENSMBLEXE = initensmbl   Calculate initial ensemble statistics
# OSSEDATAEXE  = ossedata      Generate OSSE data
# OBSSTDEXE    = obsstd        Generate STD of OSSE data
# POSTINNOVEXE = postinnov     Calculate innovation statistics
# ENDGNSEXE    = endgns        Post processing of ensemble
# TESTEXE      = test          test program executable
# ENSSCORESEXE = ensscores     Calculate verification scores
#-----------------------------------------------------------------------
#

ADASEXE      = adas
ADASMPEXE    = adas_mpi
SHIFTEXE     = arpsshift
3DVAREXE     = arps3dvar
3DVARMPEXE   = arps3dvar_mpi
ENKFEXE      = arpsenkf
ENKFMPEXE    = arpsenkf_mpi
ARPSEXE      = arps
ARPSMPEXE    = arps_mpi
AGREXE       = arpsagr
ASSIMEXE     = arpsassim
ARPSCVTEXE   = arpscvt
ARPSCVTOBSEXE   = arpscvtobs
ARPSMPBUDGETEXE = arpsmpbudget
INTRPEXE     = arpsintrp
RADARINTRPEXE= radardtaintrp
INTRP_LSEXE  = arpsintrp_ls
TINTRPEXE    = arpstintrp
TINTRPMPEXE  = arpstintrp_mpi
SUBDMNEXE    = arpssubdomain
HDFSUBDMNEXE = hdfsubdomain
HDFSUBDMNMPEXE = hdfsubdomain_mpi
#TRAJCEXE     = arpstrajc
CALCTEXE     = arpscalctrajc
PLT3DEXE     = plt3dtrajc
PLT1DEXE     = plt1dtrajc
PLTMAXEXE    = arpspltmax
PLTNCAREXE   = arpspltncar
PLTPOSTEXE   = arpspltpost
PLTPOSTMPEXE = arpspltpost_mpi
PLTNCARMPEXE = arpspltncar_mpi
PLTGRIDEXE   = pltgrid
PRTEXE       = arpsprt
A2GEXE       = arps2gem
A2NEXE       = arps2ncdf
A2WIIEXE       = arps2wdssii
A2EEXE       = arps2eta212
A2WEXE       = arps2wrf
A4WEXE       = arps4wrf
A4WMPEXE     = arps4wrf_mpi
A2CEXE       = arps2coamps
A2CMPEXE     = arps2coamps_mpi
C2AEXE       = coamps2arps
C2AMPEXE     = coamps2arps_mpi
N2AEXE       = nmm2arps
N2AMPEXE     = nmm2arps_mpi
WSEXE        = wrfstatic
A2WMPEXE     = arps2wrf_mpi
W2AEXE       = wrf2arps
W2AMPEXE     = wrf2arps_mpi
WEXTSNDEXE   = wrfextsnd
JOINWRFEXE   = joinwrfh
DIFFEXE      = arpsdiff
DIFOBSEXE    = difobs
#ENSCVEXE     = arpsenscv
#ENSBCEXE     = arpsensbc
#ENSICEXE     = arpsensic
ENKFICEXE    = arpsenkfic
RNDPRTEXE   = rndprt
ENSICMPEXE   = arpsensic_mpi
#EPOSTEXE     = arpspost
EPOSTMPEXE   = arpspost_mpi
B2GEXE       = bin2gem
#ENSANAEXE    = ens_ana
#ENSCALEXE    = ens_cal
EXTSNDEXE    = arpsextsnd
EXTSNDMPEXE  = arpsextsnd_mpi
SKEWTNEXE    = skewtncar
SKEWTPEXE    = skewtpost
ARPSSFCEXE   = arpssfc
ARPSSOILEXE  = arpssoil
88D2AEXE     = 88d2arps
88D2ASEXE    = solo2arps
88D2AMPEXE   = 88d2arps_mpi
88D2AEXE_FAKE= 88d2arps_fake
NIDS2AEXE    = nids2arps
NIDS2AMPEXE  = nids2arps_mpi
RADMOSEXE    = radmosaic
ARPS2RADEXE  = arps2rad
NCRAD2AEXE   = ncrad2arps
NCRAD2AMPEXE = ncrad2arps_mpi
CASA2AEXE    = casa2arps
RADEMULEXE   = radaremul
ATTENEXE     = atten
RADB2CDFEXE  = radbin2cdf
RADSECTEXE   = radsector
TMATRIXEXE   = tmatrix

MOSAICEXE    = Zmosaic2arps
MCI2AEXE     = mci2arps
SAT2AEXE     = sat2arps
SATHDF5AEXE  = sathdf52arps
MERGESATEXE  = mergesat
PLTRADCOLEXE = pltradcol
RADARPLTNCAREXE = radarpltncar
RADARPLTPOSTEXE = radarpltpost
PLTSATFLDEXE = pltsatfld
TERNEXE      = arpstern
TRNEXE       = arpstrn
MERGETRNEXE  = mergetrn
DIR1DEGEXE   = dir1deg
DIR5MINEXE   = dir5min
DIR30SECEXE  = dir30sec
E2AEXE       = ext2arps
E2AMPEXE     = ext2arps_mpi
AVNEXE       = extract_avn
L2AEXE       = ext2arps.laps
G2AEXE       = ext2arps.gempak
WTRETCOLEXE  = wtretcol
VERIFEXE     = arpsverif
VERIFMPEXE   = arpsverif_mpi
TYPHOONEXE   = typhoontrack

ARPSREADEXE  = arpsread
RAINDIFFEXE  = arpsraindiff

SPLITEXE     = splitfiles
SPLITMPEXE   = splitfiles_mpi
JOINSEXE     = joinfiles
JOINEXE      = joinfile
SPLITANYEXE  = splitany

H2GEXE       = hdf2grads

INITENSMBLEXE = initensmbl
OSSEDATAEXE  = ossedata
OBSSTDEXE    = obsstd
POSTINNOVEXE = postinnov
ENDGNSEXE    = endgns
TESTEXE      = test

ENSSCORESEXE = ensscores

#-----------------------------------------------------------------------
#
# Object library to be generated for ARPS solver:
#
# ARPSSOLVER  = arpssolver    ARPS solver library
# LIBARPS     = libarps       ARPS shared library
# LIBADAS     = libadas       ADAS shared library
# LIBENKF     = libenkf       EnKF shared library
# PRDLIB      = prdlib	      PRD shared library
# IPLIB       = iplib         NCEP library for verification system
# LIBA2IO     = liba2io       A2IO shared library
#
#-----------------------------------------------------------------------

ARPSSOLVER = arpssolver
LIBARPS    = libarps
LIBRADTN   = libradtn
LIBADAS    = libadas
LIBENKF    = libenkf
PRDLIB     = prdlib
IPLIB      = iplib
LIBA2IO    = liba2io
LIBZXPOST  = libzxpost
LIBZXNCAR  = libzxncar
LIBPHDF5   = libwrfio_phdf5
LIBWRFINT  = libwrfio_int
LIBWRFNET  = libwrfio_net
LIBG2      = libg2

LIBG2_DEP  = /dev/null
LIBPHDF5_OFF =
LIBPHDF5_ON  = $(LIBPHF5)
LIBPHDF5_DEP = $(LIBPHDF5_OFF)

#-----------------------------------------------------------------------
#
# Tar archive files to generated by this make file:
#
# ALLFILESTAR = Allfiles.tar    Archive all files of ARPS model
# ALLFILESGZ  = Allfiles.tar.gz Archive all files of ARPS model
# ADASTAR     = adas.tar        ADAS analysis
# 3DVARTAR    = arps3dvar.tar   3DVAR analysis
# ENKFTAR     = arpsenkf.tar    EnKF analysis
# ARPSTAR     = arps.tar        ARPS model tar file
# AGRTAR      = arpsagr.tar     ARPSAGR tar file
# ASSIMTAR    = arpsassim.tar   Assimulation version of the ARPS model
# CVTTAR      = arpscvt.tar     History data dump format convesion
# CVTOBSTAR   = arpscvtobs.tar  Archive ARPSCVTOBS program.
# BUDGTAR     = arpsmpbudget.tar Archive ARPSMPBUDGET program
# INTRPTAR    = arpsintrp.tar   Gridded data interpolation
# TINTRPTAR   = arpstintrp.tar  Gridded data interpolation
# PLTTAR      = arpsplt.tar     ZXPLOT based graphic analysis program
# PRTTAR      = arpsprt.tar     Formatted table printing program
# A2GTAR      = arps2gem.tar    Convert history files to GEMPAK format
# A2NTAR      = arps2ncdf.tar   Convert history files to netCDF format
# A2ETAR      = arps2eta212.tar Convert history files to Eta 212 format
# DIFFTAR     = arpsdiff.tar    Difference of ARPS history files
# ENSTAR      = arpsens.tar     ARPS ensamble programs
# EXTSNDTAR   = arpsextsnd.tar  External data to ARPS sounding file
# SFCTAR      = arpssfc.tar     ARPS surface data pre-processor
# SOILTAR     = arpssoil.tar    ARPS soil data pre-processor
# TERNTAR     = arpstern.tar    ARPS terrain data pre-processor
# TRNTAR      = arpstrn.tar     ARPS new terrain data pre-processor
# E2ATAR      = ext2arps.tar    External data to ARPS grid
# A2WTAR      = arps2wrf.tar    Convert ARPS files to WRF startup files
# W2ATAR      = wrf2arps.tar    Convert WRF outputs to ARPS format
# 88DTAR      = 88d2arps.tar    Tar file for 88d2arps
# MCITAR      = mci2arps.tar    Tar file for mci2arps
# ARPSVERTAR  = arpsverif.tar   Tar file for arpsverif
# ARPSMPTAR   = arps_mp.tar     Parallel ARPS executable
#
# READTAR     = arpsread.tar    ARPSREAD
# RAINDIFFTAR = arpsraindiff.tar ARPSRAINDIFF
# MERGETRNTAR = mergetrn.tar    ARPS terrain data merger/smoother
#
# OSSEDATATAR = ossedata.tar    OSSEDATA
# OBSSTDTAR   = obsstd.tar      OBSSTD
# POSTINNOVTAR = postinnov.tar  POSTINNOV
# ENDGNSTAR   = endgns.tar      ENDGNS
#-----------------------------------------------------------------------

ALLTAR      = all-tar
ALLFILESTAR = Allfiles.tar
ALLFILESGZ  = Allfiles.tar.gz
ARPSTAR     = arps.tar
ADASTAR     = adas.tar
3DVARTAR    = arps3dvar.tar
ENKFTAR     = arpsenkf.tar
AGRTAR      = arpsagr.tar
ASSIMTAR    = arpsassim.tar
CVTTAR      = arpscvt.tar
CVTOBSTAR   = arpscvtobs.tar
BUDGTAR     = arpsmpbuget.tar
INTRPTAR    = arpsintrp.tar
TINTRPTAR   = arpstintrp.tar
TRAJCTAR    = arpstrajc.tar
PLTTAR      = arpsplt.tar
PRTTAR      = arpsprt.tar
A2GTAR      = arps2gem.tar
A2NTAR      = arps2ncdf.tar
A2ETAR      = arps2eta212.tar
DIFFTAR     = arpsdiff.tar
DIFOBSTAR   = difobs.tar
ARPS2RADTAR = arps2rad.tar
RADMOSTAR   = radmostar.tar
RADEMULTAR  = rademultar.tar
RADB2CDFTAR = radb2cdf.tar
RADSECTTAR  = radsector.tar
ATTENTAR    = atten.tar
ENSTAR      = arpsens.tar
EXTSNDTAR   = arpsextsnd.tar
SKEWTTAR    = skewt.tar
SFCTAR      = arpssfc.tar
SOILTAR     = arpssoil.tar
TERNTAR     = arpstern.tar
TRNTAR      = arpstrn.tar
E2ATAR      = ext2arps.tar
A2WTAR      = arps2wrf.tar
W2ATAR      = wrf2arps.tar
AVNTAR      = extract_avn.tar
88DTAR      = 88d2arps.tar
MCITAR      = mci2arps.tar
MOSAICTAR   = mosaic.tar
ARPSVERTAR  = arpsverif.tar
ARPSMPTAR   = arps_mp.tar

READTAR     = arpsread.tar
RAINDIFFTAR = arpsraindiff.tar
MERGETRNTAR = mergetrn.tar
WTRETCOLTAR = wtretcol.tar
ZXPLOTTAR   = zxplot.tar

OSSEDATATAR = ossedata.tar
OBSSTDTAR   = obsstd.tar
POSTINNOVTAR = postinnov.tar
ENDGNSTAR = endgns.tar

MLSOTAR   = meso2lso.tar

#-----------------------------------------------------------------------
#
# List of general Documentation/Information Files:
#
#-----------------------------------------------------------------------

RELEASE    = RELEASE.NOTES
README     = README $(DATA_DIR)/README.data
HISTORY    = HISTORY
MANIFESTS  = MANIFESTS
UPDATE     = UPDATE TODO
BUGS       = BUGS

ROOTFILES1 = $(RELEASE) $(README) $(HISTORY) $(MANIFESTS) \
             $(UPDATE) $(BUGS)

#-----------------------------------------------------------------------
#
# List of root Makefiles
#
#-----------------------------------------------------------------------

MAKEFILE  = Makefile Makefile.wrkdir
MAKESCRPT = makearps

#-----------------------------------------------------------------------
#
# List of files in data subdirectory
#
#-----------------------------------------------------------------------

DATAREADME = $(DATA_DIR)/README.data

#-----------------------------------------------------------------------
#
# List of sounding files
#
#-----------------------------------------------------------------------

SNDFILES = $(SND_DIR)/may20.snd      $(SND_DIR)/may20_calm.snd \
           $(SND_DIR)/neutral.snd    $(SND_DIR)/stable.snd     \
           $(SND_DIR)/july21_ice.snd $(SND_DIR)/fife.snd       \
           $(SND_DIR)/wangara.snd

#-----------------------------------------------------------------------
#
# List of script files
#
#-----------------------------------------------------------------------

SCRIPTS1 = $(SCRPT_DIR)/linkxlC       $(SCRPT_DIR)/update_arps \
           $(SCRPT_DIR)/nqs.arps      $(SCRPT_DIR)/getdiff    \
           $(SCRPT_DIR)/runflint      $(SCRPT_DIR)/nqs.t3d    \
           $(SCRPT_DIR)/mkarpswd

GRDSSCRPT = $(SCRPT_DIR)/arps.gs              \
            $(SCRPT_DIR)/sfcflx.gs            \
            $(SCRPT_DIR)/density_grads.r      \
            $(SCRPT_DIR)/beltrami_grads.r     \
            $(SCRPT_DIR)/surface_grads.r      \
            $(SCRPT_DIR)/template_sfcdata.ctl \
            $(SCRPT_DIR)/template_trndata.ctl \
            $(SCRPT_DIR)/template_soilvar.ctl \
            $(SCRPT_DIR)/openfn.gs            \
            $(SCRPT_DIR)/openhdf.gs           \
            $(SCRPT_DIR)/useg

TSTSCRPT = $(SCRPT_DIR)/testmkarps         \
           $(SCRPT_DIR)/getincpar.pm       \
           $(SCRPT_DIR)/mkinput.pm         \
           $(SCRPT_DIR)/README_plscripts   \
           $(SCRPT_DIR)/link_data          \
           $(SCRPT_DIR)/testarps.pl        \
           $(SCRPT_DIR)/test_adas.pl       \
           $(SCRPT_DIR)/test_agri.pl       \
           $(SCRPT_DIR)/test_arps2wrf.pl   \
           $(SCRPT_DIR)/test_arpsagr.pl    \
           $(SCRPT_DIR)/test_arpsmpi.pl    \
           $(SCRPT_DIR)/test_arpssfc.pl    \
           $(SCRPT_DIR)/test_arpstern.pl   \
           $(SCRPT_DIR)/test_arpstrn.pl    \
           $(SCRPT_DIR)/test_beltrami.pl   \
           $(SCRPT_DIR)/test_density.pl    \
           $(SCRPT_DIR)/test_ext2arps.pl   \
           $(SCRPT_DIR)/test_nids2arps.pl  \
           $(SCRPT_DIR)/test_may20delcitystorm.pl    \
           $(SCRPT_DIR)/test_mci2arps.pl   \
           $(SCRPT_DIR)/test_realcase.pl   \
           $(SCRPT_DIR)/test_restart.pl    \
           $(SCRPT_DIR)/test_surface.pl    \
           $(SCRPT_DIR)/test_symmetry.pl   \
           $(SCRPT_DIR)/test_warmbubble.pl \
           $(SCRPT_DIR)/test_arpsintrp.pl  \
           $(SCRPT_DIR)/test_intrp.pl

SCRIPTS = $(SCRIPTS1) $(TSTSCRPT) $(GRDSSCRPT)

#-----------------------------------------------------------------------
#
# List of ARPS documentation files in sub-directory docs
#
#-----------------------------------------------------------------------

DOCREADME = $(DOC_DIR)/README.docs      \
            $(DOC_DIR)/README.grib2arps \
            $(DOC_DIR)/README.GrADS     \
            $(DOC_DIR)/README.plscript  \
            $(DOC_DIR)/ARPSQuickGuide.pdf   \
            $(DOC_DIR)/88d2arps50.pdf $(DOC_DIR)/a2io.txt            \
            $(DOC_DIR)/ADAS500.pdf    $(DOC_DIR)/ADASDataFormats.pdf \
            $(DOC_DIR)/ADASErrorFormats.pdf $(DOC_DIR)/ADASNudging.pdf \
            $(DOC_DIR)/ADASPhaseAdjust.pdf  $(DOC_DIR)/ARPSNetCDF.pdf  \
            $(DOC_DIR)/arps2wrf.pdf   $(DOC_DIR)/arps50.ch3.pdf        \
            $(DOC_DIR)/arps50.ch4.pdf $(DOC_DIR)/arps50.ch8.pdf        \
            $(DOC_DIR)/mci2arps.pdf   $(DOC_DIR)/ARPSShading.pdf       \
            $(DOC_DIR)/arpsshift.pdf  $(DOC_DIR)/newprogs.txt         \
            $(DOC_DIR)/nids2arps.pdf  $(DOC_DIR)/TerrainTerm.pdf   \
            $(DOC_DIR)/casa2arpsppi.pdf

TREEDOC = $(DOC_DIR)/all.trees \
          $(DOC_DIR)/arps.tree \
          $(DOC_DIR)/adas.tree \
          $(DOC_DIR)/ext2arps.tree

DOCFILES = $(DOCREADME) $(TREEDOC)

#-----------------------------------------------------------------------
#
# List of all files in ARPS root directory
#
#-----------------------------------------------------------------------

ROOTFILES2 = $(MAKEFILE) $(MAKESCRPT) $(README_DATA)

ROOTFILES = $(ROOTFILES1) $(ROOTFILES2)

ROOTALL = $(ROOTFILES) $(SNDFILES)  $(SCRIPTS) $(DOCFILES)

#-----------------------------------------------------------------------
#
# List of ARPS files in sub-directory src for each program individually
#
#-----------------------------------------------------------------------

ARPSREADME = $(ARPS_DIR)/README.arps

ARPSMKR  = $(ARPS_DIR)/Makefile
ARPSINPUT = $(INPUT_DIR)/arps.input

ARPSSRC = $(ARPS_DIR)/arps.f90        $(ARPS_DIR)/advct3d.f90     \
          $(ARPS_DIR)/advfct3d.f90    $(ARPS_DIR)/bc3d.f90        \
          $(ARPS_DIR)/bcdif3d.f90     $(ARPS_DIR)/celtrk3d.f90    \
          $(ARPS_DIR)/chksym3d.f90    $(ARPS_DIR)/cmix3d.f90      \
          $(ARPS_DIR)/cumulus3d.f90   $(ARPS_DIR)/energy3d.f90    \
          $(ARPS_DIR)/exbc3d.f90      $(ARPS_DIR)/fft99f.f90      \
          $(ARPS_DIR)/force3d.f90     $(ARPS_DIR)/grdtrns3d.f90   \
          $(ARPS_DIR)/img3d.f90       $(ARPS_DIR)/initpara3d.f90  \
          $(ARPS_DIR)/init3d.f90      $(ARPS_DIR)/inibase3d.f90   \
          $(ARPS_DIR)/irrad3d.f90     $(ARPS_DIR)/kfinterfc.f90   \
          $(ARPS_DIR)/kfpara.f90      $(ARPS_DIR)/maproj3d.f90    \
          $(ARPS_DIR)/micro3d.f90     $(ARPS_DIR)/micro_ice3d.f90 \
          $(ARPS_DIR)/micro_nem3d.f90 $(ARPS_DIR)/operat3d.f90    \
          $(ARPS_DIR)/micro_lfo_ice.f90 $(ARPS_DIR)/rhouvw.f90    \
          $(ARPS_DIR)/out3d.f90       $(ARPS_DIR)/raydmp3d.f90    \
          $(ARPS_DIR)/raddata3d.f90   $(ARPS_DIR)/radfrc3d.f90    \
          $(ARPS_DIR)/radlib3d.f90    $(ARPS_DIR)/radtrns3d.f90   \
          $(ARPS_DIR)/rst3d.f90       $(ARPS_DIR)/setbdt3d.f90    \
          $(ARPS_DIR)/sfcphy3d.f90    $(ARPS_DIR)/soildiag3d.f90  \
          $(ARPS_DIR)/soilebm3d.f90   $(ARPS_DIR)/solve3d.f90     \
          $(ARPS_DIR)/sorad3d.f90     $(ARPS_DIR)/tinteg3d.f90    \
          $(ARPS_DIR)/tke3d.f90       $(ARPS_DIR)/tmix3d.f90      \
          $(ARPS_DIR)/vfftpack.f90    $(ARPS_DIR)/wwcont3d.f90    \
          $(ARPS_DIR)/nudge.f90       $(ARPS_DIR)/smooth3d.f90    \
          $(ARPS_DIR)/assimdummy.f90  $(ARPS_DIR)/module_cu_bmj.f90 \
          $(ARPS_DIR)/module_mp_wsm6.F $(ARPS_DIR)/micro_wsm6.f90 \
          $(ARPS_DIR)/module_cu_kfeta.f90 $(ARPS_DIR)/initlib3d.f90 \
          $(ARPS_DIR)/interface_wrf_bmj.f90                         \
          $(ARPS_DIR)/interface_wrf_kfeta.f90 $(ARPS_DIR)/micro_MY.f90

IOSRC = $(ARPS_DIR)/gradsio3d.f90     $(ARPS_DIR)/gribio3d.f90  \
            $(ARPS_DIR)/gribdec3d.f90 $(ARPS_DIR)/gribenc3d.f90 \
            $(ARPS_DIR)/hdfio3d.f90   $(ARPS_DIR)/netio3d.f90   \
            $(ARPS_DIR)/pakio3d.f90   $(ARPS_DIR)/sviio3d.f90   \
            $(ARPS_DIR)/binio.c       $(ARPS_DIR)/v5dio3d.f90   \
            $(ARPS_DIR)/v5d.c         $(ARPS_DIR)/nohdfio3d.f90 \
            $(ARPS_DIR)/nopakio3d.f90 $(ARPS_DIR)/nonetio3d.f90 \
            $(ARPS_DIR)/nosviio3d.f90 $(ARPS_DIR)/nov5dio3d.f90 \
            $(ARPS_DIR)/iolib3d.f90                             \
            $(ARPS_DIR)/exbcio3d.f90  $(ARPS_DIR)/arpsio3d.f90  \
            $(ARPS_DIR)/binio3d.f90   $(ARPS_DIR)/ascio3d.f90   \
            $(ARPS_DIR)/module_arps_netio_metadata.f90

LIBSRC = $(ARPS_DIR)/craylib3d.f90 $(ARPS_DIR)/htchlib3d.f90 \
         $(ARPS_DIR)/genlib3d.f90  $(ARPS_DIR)/ibmlib3d.f90  \
         $(ARPS_DIR)/arpslib3d.f90 $(ARPS_DIR)/outlib3d.f90  \
         $(ARPS_DIR)/sunlib3d.f90  $(ARPS_DIR)/timelib3d.f90 \
         $(ARPS_DIR)/lnxlib3d.F    $(ARPS_DIR)/irixlib3d.f90 \
         $(ARPS_DIR)/t3dlib3d.f90  $(ARPS_DIR)/thermolib3d.f90 \
         $(ARPS_DIR)/nompsubs.f90  $(ARPS_DIR)/mpisubs.f90     \
         $(ARPS_DIR)/tru64lib3d.f90 $(ARPS_DIR)/alloclib.f90 \
         $(ARPS_DIR)/module_arbitrary_vario.f90              \
         $(ARPS_DIR)/module_precision.F $(ARPS_DIR)/mkdir.c  \
         $(ARPS_DIR)/my3mom_main_mod.f90 \
         $(ARPS_DIR)/my3mom_fncs_mod.f90 $(ARPS_DIR)/my3mom_sedi_mod.f90 \
         $(ARPS_DIR)/my2mom_fncs_mod.f90

NCARGSRC = $(ARPS_DIR)/ncarg3d.f90 $(ARPS_DIR)/noncarg3d.f90

ARPSINC = $(INCL_DIR)/bndry.inc    $(INCL_DIR)/cumucst.inc   \
          $(INCL_DIR)/exbc.inc     $(INCL_DIR)/globcst.inc   \
          $(INCL_DIR)/indtflg.inc  $(INCL_DIR)/meraf.inc     \
          $(INCL_DIR)/nemcst.inc   $(INCL_DIR)/phycst.inc    \
          $(INCL_DIR)/radcst.inc   $(INCL_DIR)/sfcphycst.inc \
          $(INCL_DIR)/soilcst.inc  $(INCL_DIR)/nudging.inc   \
          $(INCL_DIR)/binio.h      $(INCL_DIR)/v5d.h         \
          $(INCL_DIR)/v5df.inc     $(INCL_DIR)/vis5d.h       \
          $(INCL_DIR)/mp.inc       $(INCL_DIR)/alloc.inc     \
          $(INCL_DIR)/timelvls.inc $(INCL_DIR)/grid.inc      \
          $(INCL_DIR)/vericst.inc  $(INCL_DIR)/dfilter.inc   \
          $(INCL_DIR)/GEMPRM.AIX.inc     $(INCL_DIR)/GEMPRM.HPUX.inc    \
          $(INCL_DIR)/GEMPRM.IRIX.inc    $(INCL_DIR)/GEMPRM.LINUX.inc   \
          $(INCL_DIR)/GEMPRM.OSF1.inc    $(INCL_DIR)/GEMPRM.SunOS.inc   \
          $(INCL_DIR)/GEMPRM.ULTRIX.inc

#EMK BMJ
LIBARPSSRC = $(ARPS_DIR)/bc3d.f90       $(ARPS_DIR)/celtrk3d.f90 \
             $(ARPS_DIR)/fft99f.f90     $(ARPS_DIR)/rhouvw.f90   \
             $(ARPS_DIR)/grdtrns3d.f90  $(ARPS_DIR)/img3d.f90    \
             $(ARPS_DIR)/inibase3d.f90  $(ARPS_DIR)/initlib3d.f90 \
             $(ARPS_DIR)/initpara3d.f90 $(ARPS_DIR)/maproj3d.f90 \
             $(ARPS_DIR)/operat3d.f90   $(ARPS_DIR)/rst3d.f90    \
             $(ARPS_DIR)/setbdt3d.f90   $(ARPS_DIR)/vfftpack.f90 \
             $(ARPS_DIR)/wwcont3d.f90   $(ARPS_DIR)/nudge.f90    \
             $(ARPS_DIR)/smooth3d.f90   $(ARPS_DIR)/assimdummy.f90 \
             $(ARPS_DIR)/module_cu_bmj.f90     \
             $(ARPS_DIR)/interface_wrf_bmj.f90 \
             $(ARPS_DIR)/module_cu_kfeta.f90   \
             $(ARPS_DIR)/interface_wrf_kfeta.f90

READSRC = $(ARPS_DIR)/arpsread.f90

MPBUDGETSRC = $(ARPS_DIR)/arpsmpbudget.f90 $(ARPS_DIR)/arpsmpbudgetlib.f90

LIBARPSALL = $(ARPSMKR)  $(LIBARPSSRC) $(IOSRC) $(LIBSRC) \
             $(NCARGSRC) $(ARPSINC)

ARPSALL = $(ARPSREADME) $(ARPSMKR)  $(ARPSSRC) $(IOSRC)   \
          $(LIBSRC)     $(NCARGSRC) $(READSRC) $(ARPSINC) \
          $(ASSIMINC)   $(ARPSINPUT) $(MPBUDGETSRC)

#-----------------------------------------------------------------------
#
# List of ARPSAGR files in sub-directory arpsagr
#
#-----------------------------------------------------------------------

AGRREADME = $(AGR_DIR)/README.arpsagr

AGRMKR    = $(AGR_DIR)/Makefile
AGRNPUT   = $(ARPSINPUT)

AGRSRC     = $(AGR_DIR)/arpsagr.f90    $(AGR_DIR)/arpscnst.f90    \
             $(AGR_DIR)/arpsinit.f90   $(AGR_DIR)/arpsolve.f90    \
             $(AGR_DIR)/arpsout.f90    $(AGR_DIR)/arpsplot.f90    \
             $(AGR_DIR)/arpstgrid.f90  $(AGR_DIR)/change_cnst.f90 \
             $(AGR_DIR)/cpucray.f90    $(AGR_DIR)/cpugen.f90      \
             $(AGR_DIR)/cpuibm.f90     $(AGR_DIR)/cpusun.f90      \
             $(AGR_DIR)/filgrd.f90     $(AGR_DIR)/fill.f90        \
             $(AGR_DIR)/initcgrd.f90   $(AGR_DIR)/inivar.f90      \
             $(AGR_DIR)/manage.f90     $(AGR_DIR)/packcray.f90    \
             $(AGR_DIR)/packgen.f90    $(AGR_DIR)/pointf.f90      \
             $(AGR_DIR)/regrid.f90     $(AGR_DIR)/setstr.f90      \
             $(AGR_DIR)/tick.f90       $(AGR_DIR)/updbc.f90       \
             $(AGR_DIR)/updexbc.f90    $(AGR_DIR)/updgrd.f90      \
             $(AGR_DIR)/usrout.f90     $(AGR_DIR)/util1.f90

AGRINC     = $(INCL_DIR)/agricst.inc  $(INCL_DIR)/agrialloc.inc \
             $(INCL_DIR)/agricpu.inc  $(INCL_DIR)/grddsc.inc    \
             $(INCL_DIR)/agrigrid.inc $(INCL_DIR)/manage.inc    \
             $(INCL_DIR)/nodal.inc

AGRALL = $(AGRREADME) $(AGRMKR) $(AGRSRC) $(AGRINC) $(AGRINPUT)

#-----------------------------------------------------------------------
#
# List of ASSIM files in sub-directory assim
#
#-----------------------------------------------------------------------

ASSIMREADME = $(ASSIM_DIR)/README.assim

ASSIMMKR    = $(ASSIM_DIR)/Makefile
ASSIMINPUT  = $(ARPSINPUT)

ASSIMSRC    = $(ASSIM_DIR)/arpsassim.f90  $(ASSIM_DIR)/assimadasread.f90 \
              $(ASSIM_DIR)/assimchk.f90   $(ASSIM_DIR)/assimcon.f90      \
              $(ASSIM_DIR)/assimdriv.f90  $(ASSIM_DIR)/assimfil.f90      \
              $(ASSIM_DIR)/assimfz.f90    $(ASSIM_DIR)/assimout.f90      \
              $(ASSIM_DIR)/assimpois.f90  $(ASSIM_DIR)/assimpsolver.f90  \
              $(ASSIM_DIR)/assimptpr.f90  $(ASSIM_DIR)/assimradread.f90  \
              $(ASSIM_DIR)/assimrd.f90    $(ASSIM_DIR)/assimthermo.f90   \
              $(ASSIM_DIR)/assimtinpl.f90 $(ASSIM_DIR)/assimvel.f90      \
              $(ASSIM_DIR)/assimvrfil.f90 $(ASSIM_DIR)/head3d.f90        \
              $(ASSIM_DIR)/assimblnd.f90  $(ASSIM_DIR)/initassim.f90

ASSIMINC    = $(INCL_DIR)/assim.inc

ASSIMALL = $(ASSIMREADME) $(ASSIMMKR) $(ASSIMSRC) $(ASSIMINC) \
           $(ASSIMINPUT)  $(INPUT_DIR)/arpsassim.input

#-----------------------------------------------------------------------
#
# List of ADAS files in sub-directory adas
#
#-----------------------------------------------------------------------

ADASREADME = $(ADAS_DIR)/README.adas

ADASMKR    = $(ADAS_DIR)/Makefile

ADASSRC    = $(ADAS_DIR)/adas_3dvar_driver.F  $(ADAS_DIR)/initadas.f90   \
             $(ADAS_DIR)/anxiter.f90     $(ADAS_DIR)/anxlib3d.f90   \
             $(ADAS_DIR)/adjuvw3d.f90    $(ADAS_DIR)/grd2obsth.F    \
             $(ADAS_DIR)/prepradar.F     $(ADAS_DIR)/prepretr.f90   \
             $(ADAS_DIR)/prepsfc.f90     $(ADAS_DIR)/prepsng.f90    \
             $(ADAS_DIR)/prepua.f90      $(ADAS_DIR)/incrdump.f90   \
             $(ADAS_DIR)/rdacars.f90     $(ADAS_DIR)/rdprofiles.f90 \
             $(ADAS_DIR)/rdrdapstern.f90 $(ADAS_DIR)/rdradcol.f90   \
             $(ADAS_DIR)/rdretcol.f90    $(ADAS_DIR)/rdsfcobs.f90   \
             $(ADAS_DIR)/rdructern.f90   $(ADAS_DIR)/adjtsfc.f90    \
             $(ADAS_DIR)/fsl2snd.f90                                \
             $(ADAS_DIR)/module_analysisArrays.f90 \
             $(ADAS_DIR)/module_anaIncArray.f90    \
             $(ADAS_DIR)/module_arpsArray.f90      \
             $(ADAS_DIR)/module_soilArray.f90

ADASCLDSRC = $(ADAS_DIR)/cloud_cv.f90   $(ADAS_DIR)/cloud_lwc.f90   \
             $(ADAS_DIR)/cmpclddrv.f90  $(ADAS_DIR)/cmpcldlib.f90   \
             $(ADAS_DIR)/cldinsert.f90  $(ADAS_DIR)/rdsatfld.f90    \
             $(ADAS_DIR)/sunfuncs.f90

ADASINC    = $(INCL_DIR)/adas.inc     $(INCL_DIR)/remap.inc        \
             $(INCL_DIR)/remapsum.inc $(INCL_DIR)/remaptab.inc     \
             $(INCL_DIR)/nudging.inc  $(INCL_DIR)/adjust.inc       \
             $(INCL_DIR)/adassat.inc

LIBADASSRC = $(ADAS_DIR)/anxlib3d.f90   $(ADAS_DIR)/thermo3d.f90   \
             $(ADAS_DIR)/adjuvw3d.f90   $(ADAS_DIR)/intfield.f90   \
             $(ADAS_DIR)/radarlib3d.f90 $(ADAS_DIR)/mthermo.f90    \
             $(ADAS_DIR)/intrpsoil3d.f90 $(ADAS_DIR)/pltradlib.f90

LIBADASALL = $(ADASMKR) $(LIBADASSRC) $(ADASINC)

ADASTABS   = $(ATAB_DIR)/blacklist.sfc         \
             $(ATAB_DIR)/armmnerr.adastab      \
             $(ATAB_DIR)/buoyerr.adastab       \
             $(ATAB_DIR)/fclasserr.adastab     \
             $(ATAB_DIR)/mclasserr.adastab     \
             $(ATAB_DIR)/mesoerr.adastab       \
             $(ATAB_DIR)/moblmerr.adastab      \
             $(ATAB_DIR)/profilerr.adastab     \
             $(ATAB_DIR)/rad88Derr.adastab     \
             $(ATAB_DIR)/radnidserr.adastab    \
             $(ATAB_DIR)/raoberr.adastab       \
             $(ATAB_DIR)/ret88Derr.adastab     \
             $(ATAB_DIR)/ruc3herr.adastab      \
             $(ATAB_DIR)/saoerr.adastab        \
             $(ATAB_DIR)/shiperr.adastab       \
             $(ATAB_DIR)/synoperr.adastab      \
             $(ATAB_DIR)/testerr.adastab       \
             $(ATAB_DIR)/ircalib.adastab       \
             $(ATAB_DIR)/asoserr.adastab       \
             $(ATAB_DIR)/asoserr.adastab       \
             $(ATAB_DIR)/awoserr.adastab       \
             $(ATAB_DIR)/coagmeterr.adastab    \
             $(ATAB_DIR)/hplainserr.adastab    \
             $(ATAB_DIR)/iswserr.adastab       \
             $(ATAB_DIR)/radarinfo.dat         \
             $(ATAB_DIR)/goes08ch4.adastab  $(ATAB_DIR)/goes08vis.adastab  \
             $(ATAB_DIR)/goes09ch4.adastab  $(ATAB_DIR)/goes09vis.adastab  \
             $(ATAB_DIR)/goes10ch4.adastab  $(ATAB_DIR)/goes10vis.adastab  \
             $(ATAB_DIR)/goes11ch4.adastab  $(ATAB_DIR)/goes11vis.adastab  \
             $(ATAB_DIR)/goes12ch4.adastab  $(ATAB_DIR)/goes12vis.adastab  \
             $(ATAB_DIR)/mdcrserr.adastab   $(ATAB_DIR)/wtxmnerr.adastab

ADASALL = $(ADASREADME) $(ADASMKR) $(ADASSRC)   $(LIBADASSRC) \
          $(ADASCLDSRC) $(ADASINC) $(ADASTABS)

#-----------------------------------------------------------------------
#
# List of SHIFTEXE files in sub-directory adas
#
#-----------------------------------------------------------------------

MLSOALL = $(MLSO_DIR)/README.meso2lso $(MLSO_DIR)/Makefile \
          $(MLSO_DIR)/citime.f $(MLSO_DIR)/cvtitoi.c $(MLSO_DIR)/getime.f \
          $(MLSO_DIR)/meso2lso.f $(MLSO_DIR)/mkitime.c $(MLSO_DIR)/readmdf.f \
          $(MLSO_DIR)/readmeso.f $(MLSO_DIR)/readwtx.f $(MLSO_DIR)/strlnth.f \
          $(MLSO_DIR)/upcase.c

#-----------------------------------------------------------------------
#
# List of SHIFTEXE files in sub-directory adas
#
#-----------------------------------------------------------------------

SHIFTREADME = $(ADAS_DIR)/README.arpsshift

SHIFTINPUT  = $(INPUT_DIR)/arpsshift.input

SHIFTSRC    = $(ADAS_DIR)/arpsshift.f90 \
              $(ADAS_DIR)/initshift.f90 \
              $(ADAS_DIR)/prepradsh.f90 \
              $(ADAS_DIR)/grd2obssh.f90 \
              $(ADAS_DIR)/rshift3d.f90  \
              $(ADAS_DIR)/shftsub.f90   \
              $(ADAS_DIR)/shiftlib.f90  \
              $(ADAS_DIR)/dummyshift.f90

SHIFTALL = $(SHIFTREADME) $(SHIFTSRC) $(SHIFTINPUT)

#-----------------------------------------------------------------------
#
# List of 3DVAR files in sub-directory src/arps3dvar
#
#-----------------------------------------------------------------------

3DVARREADME = $(3DVAR_DIR)/README.3dvar $(3DVAR_DIR)/README.history

3DVARMKR    = $(3DVAR_DIR)/Makefile

3DVARINC    = $(INCL_DIR)/3dvarcputime.inc $(INCL_DIR)/varpara.inc

3DVARINPUT  = $(INPUT_DIR)/arps.input

3DVARSRC    = $(3DVAR_DIR)/Vbck_scale.f90  $(3DVAR_DIR)/Vcostf3d.f90        \
              $(3DVAR_DIR)/Vctr_to_vbl.f90 $(3DVAR_DIR)/Vderac_scale.f90    \
              $(3DVAR_DIR)/Vgradc3d.f90    $(3DVAR_DIR)/Vinit3dvar.f90      \
              $(3DVAR_DIR)/Vlinear_int2d.f90 $(3DVAR_DIR)/Vlinear_int3d.f90 \
              $(3DVAR_DIR)/Vminimization.f90 $(3DVAR_DIR)/Vminim_sub.f90    \
              $(3DVAR_DIR)/Vmodule_grd.f90   $(3DVAR_DIR)/Vrecur_filt1d.f90 \
              $(3DVAR_DIR)/Vrecur_filt3d.f90 $(3DVAR_DIR)/Vrecur_filtmpi.f90 \
              $(3DVAR_DIR)/Vrecur_filtnompi.f90 $(3DVAR_DIR)/Vtrans.f90      \
              $(3DVAR_DIR)/Vvbl_to_ctr.f90 $(3DVAR_DIR)/Vmodule_3dvar.f90

3DVARALL = $(3DVARREADME) $(3DVARMKR) $(3DVARINC) $(3DVARINPUT) $(3DVARSRC)

#-----------------------------------------------------------------------
#
# List of ENKF files in sub-directory EnKF (mtong)
#
#-----------------------------------------------------------------------

ENKFREADME = $(ENKF_DIR)/README.enkf

ENKFMKR    = $(ENKF_DIR)/Makefile

ENKFINPUT = $(INPUT_DIR)/arpsenkf.input

ENKFSRC    = $(ENKF_DIR)/arpsenkf.f90        $(ENKF_DIR)/initenkf.f90  \
             $(ENKF_DIR)/enkfio3d.f90        $(ENKF_DIR)/enkfda.f90    \
             $(ENKF_DIR)/enkfhelp.f90        $(ENKF_DIR)/obsio.f90     \
             $(ENKF_DIR)/obsoperator.f90     $(ENKF_DIR)/radarlib.f90  \
             $(ENKF_DIR)/stateest.f90        $(ENKF_DIR)/paraest.f90   \
             $(ENKF_DIR)/errscore.f90        $(ENKF_DIR)/rmserror.f90  \
             $(ENKF_DIR)/innovation.f90      $(ENKF_DIR)/covinflt.f90  \
             $(ENKF_DIR)/covvar.f90          $(ENKF_DIR)/attenuate.f90 \
             $(ENKF_DIR)/report.f90          $(ENKF_DIR)/module_mpi_enkf.F

LIBENKFSRC = $(ENKF_DIR)/enkflib.f90

INITENSMBLSRC = $(ENKF_DIR)/initensmbl.f90

OSSEDATASRC = $(ENKF_DIR)/ossedata.f90

OBSSTDSRC = $(ENKF_DIR)/obsstd.f90

POSTINNOVSRC = $(ENKF_DIR)/postinnov.f90

ENDGNS = $(ENKF_DIR)/endgns.f90

TESTSRC = = $(ENKF_DIR)/test.f90

ENKFINC    = $(INCL_DIR)/globcst.inc    $(INCL_DIR)/enkf.inc          \
             $(INCL_DIR)/phycst.inc     $(INCL_DIR)/indtflg.inc       \
             $(INCL_DIR)/grid.inc       $(INCL_DIR)/mp.inc

ENKFALL = $(ENKFREADME) $(ENKFMKR) $(ENKFSRC) $(LIBENKFSRC) $(POSTINNOVSRC)  \
          $(INITENSMBLSRC) $(OSSEDATASRC) $(ENDGNSSRC) $(ENKFINC)      \
          $(ENKFINPUT)

#-----------------------------------------------------------------------
#
# List of ARPSDIFF files in sub-directory arpsdiff
#
#-----------------------------------------------------------------------

DIFFREADME = $(DIF_DIR)/README.arpsdiff

DIFFMKR    = $(DIF_DIR)/Makefile
DIFFINPUT  = $(INPUT_DIR)/arpsdiff.input $(INPUT_DIR)/arpsraindiff.input

DIFFSRC    = $(DIF_DIR)/arpsdiff.f90  $(DIF_DIR)/arpsraindiff.f90

DIFFINC    =

DIFFALL = $(DIFFREADME) $(DIFFMKR) $(DIFFSRC) $(DIFFINC) $(DIFFINPUT)

#-----------------------------------------------------------------------
#
# List of DIFOBS files in sub-directory adas
#
#-----------------------------------------------------------------------

#DIFOBSREADME = $(ADAS_DIR)/README.difobs
DIFOBSSRC    = $(ADAS_DIR)/difobs.f90      $(ADAS_DIR)/difstats.f90

DIFOBSALL    = $(ADASMKR) $(DIFOBSSRC) $(ADASSRC)                       \
               $(ADASINC) $(ADASINPUT)

#-----------------------------------------------------------------------
#
# List of ARPSENS files in sub-directory arpsens
#
#-----------------------------------------------------------------------

ENSREADME = $(ENS_DIR)/README.arpsens

ENSMKR    = $(ENS_DIR)/Makefile
ENSINPUT  = $(INPUT_DIR)/arpsenscv.input  $(INPUT_DIR)/arpsensic.input \
            $(INPUT_DIR)/arpsensbc.input  $(INPUT_DIR)/arpspost.input  \
            $(INPUT_DIR)/bin2gem.input    $(INPUT_DIR)/ens_ana.input   \
            $(INPUT_DIR)/arpsenkfic.input $(INPUT_DIR)/ensscores.input

ENSSRC    = $(ENS_DIR)/arpsenscv.f90 $(ENS_DIR)/arpsensic.f90  \
            $(ENS_DIR)/arpsensbc.f90 $(ENS_DIR)/arpspost.f90   \
            $(ENS_DIR)/arpsenslib.f90 $(ENS_DIR)/postcore.f90  \
            $(ENS_DIR)/bin2gem.f90    $(ENS_DIR)/ens_ana.f90   \
            $(ENS_DIR)/gemio.f90      $(ENS_DIR)/nogemio.f90   \
            $(ENS_DIR)/arpspostlib.f90 $(ENS_DIR)/mrgrnk.f90   \
            $(ENS_DIR)/extrefleclib.f90 $(ENS_DIR)/ensscores.f90 \
            $(ENS_DIR)/ens_cal.f90      $(ENS_DIR)/arpsenkfic.f90 \
            $(ENS_DIR)/crtmpost.f90     $(ENS_DIR)/crtm_no.f90    \
            $(ENS_DIR)/lsq.f90          $(ENS_DIR)/logistic.f90

ENSINC    = $(INCL_DIR)/enscv.inc $(INCL_DIR)/arpsenkfic.inc \
            $(INCL_DIR)/enkf.inc

ENSALL = $(ENSREADME) $(ENSMKR) $(ENSSRC) $(ENSINC) $(ENSINPUT)

#-----------------------------------------------------------------------
#
# List of ARPSPLT files in sub-directory arpsplt
#
#-----------------------------------------------------------------------

PLTREADME = $(PLT_DIR)/README.arpsplt $(PLT_DIR)/HISTORY

PLTMKR    = $(PLT_DIR)/Makefile
PLTINPUT  = $(INPUT_DIR)/arpsplt.input $(INPUT_DIR)/arpspltmax.input \
            $(INPUT_DIR)/pltgrid.input

PLTINC    = $(INCL_DIR)/arpsplt.inc

PLTSRC    = $(PLT_DIR)/arpsplt.f90    $(PLT_DIR)/arpspltlib.f90       \
            $(PLT_DIR)/arpspltmax.f90 $(PLT_DIR)/arpsplt_cpu.f90      \
            $(PLT_DIR)/wirfrmstub.f90 $(PLT_DIR)/read_surface_obs.f90 \
            $(PLT_DIR)/wirfrm.f90     $(PLT_DIR)/arpspltderive.f90    \
            $(PLT_DIR)/pltgrid.f90    $(PLT_DIR)/plot_coltab.f90      \
            $(PLT_DIR)/mpisubs.f90    $(PLT_DIR)/nompisubs.f90

PLTDATA = $(PLTD_DIR)/us_spcounty.mapdata      \
          $(PLTD_DIR)/us_state.mapdata         \
          $(PLTD_DIR)/world_coast.mapdata      \
          $(PLTD_DIR)/world_country.mapdata    \
          $(PLTD_DIR)/world_us_country.mapdata \
          $(PLTD_DIR)/zx_color.tbl             \
          $(PLTD_DIR)/aa-dfw.meta

PLTALL = $(PLTREADME) $(PLTMKR)   $(PLTSRC) $(PLTINC) \
         $(PLTDATA)   $(PLTINPUT)

#-----------------------------------------------------------------------
#
# List of ARPSCVT files in sub-directory cvt
#
#-----------------------------------------------------------------------

CVTREADME = $(CVT_DIR)/README.arpscvt

CVTMKR    = $(CVT_DIR)/Makefile
CVTINPUT  = $(INPUT_DIR)/arpscvt.input

CVTSRC    = $(CVT_DIR)/arpscvt.f90 $(CVT_DIR)/cvtiolib3d.f90

H2GSRC    = $(CVT_DIR)/hdf2grads.f90

CVTINC    =

CVTALL = $(CVTREADME) $(CVTMKR)  $(CVTSRC) $(CVTINPUT) $(CVTINC)        \
         $(H2GSRC)

#-----------------------------------------------------------------------
#
# List of ARPSCVTOBS files in sub-directory cvt
#
#-----------------------------------------------------------------------

CVTOBSREADME = $(CVTOBS_DIR)/README.arpscvtobs

CVTOBSMKR    = $(CVTOBS_DIR)/Makefile
CVTOBSINPUT  = $(INPUT_DIR)/arpscvtobs.input

CVTOBSSRC    = $(CVTOBS_DIR)/arpscvtobs.f90 $(CVTOBS_DIR)/wrtlsohdf.f90  \
               $(CVTOBS_DIR)/module_verifhdf.f90

CVTOBSINC    =

CVTOBSALL = $(CVTOBSREADME) $(CVTOBSMKR)  $(CVTOBSSRC) $(CVTOBSINPUT) \
        $(CVTOBSINC)

#-----------------------------------------------------------------------
#
# List of ARPSINTRP and ARPSSUBDOMAIN files in sub-directory intrp
#
#-----------------------------------------------------------------------

INTRPREADME = $(INTRP_DIR)/README.intrp

INTRPMKR    = $(INTRP_DIR)/Makefile
INTRPINPUT  = $(INPUT_DIR)/arpsintrp.input $(INPUT_DIR)/arpsintrp_ls.input \
              $(INPUT_DIR)/arpssubdomain.input $(INPUT_DIR)/radardtaintrp.input

INTRPSRC    = $(INTRP_DIR)/arpsintrp.f90     \
              $(INTRP_DIR)/arpsintrp_ls.f90  \
              $(INTRP_DIR)/arpssubdomain.f90 \
              $(INTRP_DIR)/intrplib.f90      \
              $(INTRP_DIR)/check_file.f90    \
              $(INTRP_DIR)/radardtaintrp.f90 \
              $(INTRP_DIR)/hdfsubdomain.f90 $(INTRP_DIR)/hdfio.o \
              $(INTRP_DIR)/module_hdfsubdomain_namelist.f90

INTRPINC    =

INTRPALL = $(INTRPREADME) $(INTRPMKR) $(INTRPSRC) $(INTRPINC) \
           $(INTRPINPUT)

#-----------------------------------------------------------------------
# List of ARPSTRAJC files
#-----------------------------------------------------------------------

TRAJCMKR    = $(TRAJC_DIR)/Makefile $(TRAJC_DIR)/README
TRAJCINPUT  = $(INPUT_DIR)/arpstrajc.input
TRAJCSRC    = $(TRAJC_DIR)/arpstrajc.f90  $(TRAJC_DIR)/arpscalctrajc.f90 \
              $(TRAJC_DIR)/plt3dtrajc.f90 $(TRAJC_DIR)/plt1dtrajc.f90

TRAJCALL = $(TRAJCMKR) $(TRAJCINPUT) $(TRAJCSRC)


#-----------------------------------------------------------------------
#
# List of ARPSTINTRP files in sub-directory tintrp
#
#-----------------------------------------------------------------------

TINTRPREADME= $(TINTRP_DIR)/README.tintrp

TINTRPMKR   = $(TINTRP_DIR)/Makefile
TINTRPINPUT = $(INPUT_DIR)/arpstintrp.input

TINTRPSRC   = $(TINTRP_DIR)/arpstintrp.f90

TINTRPINC   =

TINTRPALL   = $(TINTRPREADME) $(TINTRPMKR) $(TINTRPSRC) $(TINTRPINC) \
              $(TINTRPINPUT)

#-----------------------------------------------------------------------
#
# List of ARPSPRT files in sub-directory prt
#
#-----------------------------------------------------------------------

PRTREADME = $(PRT_DIR)/README.arpsprt

PRTMKR    = $(PRT_DIR)/Makefile
PRTINPUT  = $(INPUT_DIR)/arpsprt.input

PRTSRC    = $(PRT_DIR)/arpsprt.f90 $(PRT_DIR)/arpsprtlib.f90

PRTALL = $(PRTREADME) $(PRTMKR) $(PRTSRC) $(PRTINPUT)

#-----------------------------------------------------------------------
#
# List of ARPSSFC files in sub-directory sfc
#
#-----------------------------------------------------------------------

SFCREADME = $(SFC_DIR)/README.arpssfc

SFCMKR    = $(SFC_DIR)/Makefile
SFCINPUT  = $(ARPSINPUT)

SFCSRC    = $(SFC_DIR)/arpssfc.f90 $(SFC_DIR)/arpssfclib.f90

SFCINC    = $(INCL_DIR)/arpssfc.inc

SFCMISC   = $(SFC_DIR)/sfc_winter.tbl

SFCALL = $(SFCREADME) $(SFCMKR) $(SFCSRC) $(SFCINC) $(SFCMISC) \
         $(SFCINPUT)

#-----------------------------------------------------------------------
#
# List of ARPSVERIF files in sub-directory arpsverif
#
#-----------------------------------------------------------------------

VERIFREADME = $(VERIF_DIR)/README.arpsverif

VERIFMKR    = $(VERIF_DIR)/Makefile
VERIFINPUT  = $(INPUT_DIR)/arpsverif.input

VERIFSRC = $(VERIF_DIR)/arpshis2ver.f90  $(VERIF_DIR)/arpsverif.f90 \
           $(VERIF_DIR)/avgverifgrid.f90 $(VERIF_DIR)/cvt2verif.f90 \
           $(VERIF_DIR)/cvt2verif.f90    $(VERIF_DIR)/extdims2.f90 \
           $(VERIF_DIR)/getgrib.f90      $(VERIF_DIR)/gribcst2.f90 \
           $(VERIF_DIR)/indtflg.f90      $(VERIF_DIR)/intfield2.f90 \
           $(VERIF_DIR)/intqpf.f90       $(VERIF_DIR)/intverifh.f90 \
           $(VERIF_DIR)/module_verif.f90 $(VERIF_DIR)/module_verifhdf.f90 \
           $(VERIF_DIR)/qpfmask.f90      $(VERIF_DIR)/qpfstats.f90 \
           $(VERIF_DIR)/rdarpsqpf.f90    $(VERIF_DIR)/rdnmcgrb3d2.f90 \
           $(VERIF_DIR)/rdverifstats.f90 $(VERIF_DIR)/setijloc2.f90 \
           $(VERIF_DIR)/verifgrid.f90    $(VERIF_DIR)/verif_hdf.f90 \
           $(VERIF_DIR)/veriflib.f90     $(VERIF_DIR)/verifstatslib.f90 \
           $(VERIF_DIR)/vgridconst.f90   $(VERIF_DIR)/arps_nn.c     \
           $(VERIF_DIR)/arpshis2verGrid.f90 $(VERIF_DIR)/typhoontrack.f90

IPLIBSRC = $(VERIF_DIR)/iplib/decode.f  $(VERIF_DIR)/iplib/encode.f \
           $(VERIF_DIR)/iplib/gausslat.f $(VERIF_DIR)/iplib/gcdist.f \
           $(VERIF_DIR)/iplib/gdsawt.f $(VERIF_DIR)/iplib/gdswiz00.f \
           $(VERIF_DIR)/iplib/gdswiz01.f $(VERIF_DIR)/iplib/gdswiz03.f \
           $(VERIF_DIR)/iplib/gdswiz04.f $(VERIF_DIR)/iplib/gdswiz05.f \
           $(VERIF_DIR)/iplib/gdswizc9.f $(VERIF_DIR)/iplib/gdswizca.f \
           $(VERIF_DIR)/iplib/gdswzd00.f $(VERIF_DIR)/iplib/gdswzd01.f \
           $(VERIF_DIR)/iplib/gdswzd03.f $(VERIF_DIR)/iplib/gdswzd04.f \
           $(VERIF_DIR)/iplib/gdswzd05.f $(VERIF_DIR)/iplib/gdswzdc9.f \
           $(VERIF_DIR)/iplib/gdswzdca.f $(VERIF_DIR)/iplib/gdswzd.f \
           $(VERIF_DIR)/iplib/ijkgds.f $(VERIF_DIR)/iplib/ipolates.f \
           $(VERIF_DIR)/iplib/ipolatev.f $(VERIF_DIR)/iplib/ipxetas.f \
           $(VERIF_DIR)/iplib/ipxwafs2.f $(VERIF_DIR)/iplib/ipxwafs.f \
           $(VERIF_DIR)/iplib/makgds.f $(VERIF_DIR)/iplib/movect.f \
           $(VERIF_DIR)/iplib/polateg0.f $(VERIF_DIR)/iplib/polateg1.f \
           $(VERIF_DIR)/iplib/polateg4.f \
           $(VERIF_DIR)/iplib/polates0.f $(VERIF_DIR)/iplib/polates1.f \
           $(VERIF_DIR)/iplib/polates2.f $(VERIF_DIR)/iplib/polates3.f \
           $(VERIF_DIR)/iplib/polates4.f $(VERIF_DIR)/iplib/polates5.f \
           $(VERIF_DIR)/iplib/polates6.f $(VERIF_DIR)/iplib/polatev0.f \
           $(VERIF_DIR)/iplib/polatev1.f $(VERIF_DIR)/iplib/polatev2.f \
           $(VERIF_DIR)/iplib/polatev3.f $(VERIF_DIR)/iplib/polatev4.f \
           $(VERIF_DIR)/iplib/polatev6.f $(VERIF_DIR)/iplib/polfixs.f \
           $(VERIF_DIR)/iplib/polfixv.f $(VERIF_DIR)/iplib/r63w72.f \
           $(VERIF_DIR)/iplib/w3fi63.f $(VERIF_DIR)/iplib/w3fi68.f \
           $(VERIF_DIR)/iplib/w3fi71.f $(VERIF_DIR)/iplib/w3fi74.f

IPLIBMKR = $(VERIF_DIR)/iplib/Makefile

VERIFDAT = $(VERD_DIR)/sfcstns.tbl

VERIFALL = $(VERIFREADME) $(VERIFMKR) $(VERIFSRC) $(IPLIBSRC) $(VERIFINPUT) \
           $(VERIFDAT) $(IPLIBMKR)

#-----------------------------------------------------------------------
#
# List of ARPSSOIL files in sub-directory soil
#
#-----------------------------------------------------------------------

SOILREADME = $(SOIL_DIR)/README.arpssoil

SOILMKR    = $(SOIL_DIR)/Makefile
SOILINPUT  = $(INPUT_DIR)/arpssoil.input

SOILSRC    = $(SOIL_DIR)/arpssoil.f90    $(SOIL_DIR)/barnes3d.f90  \
             $(SOIL_DIR)/mksoilvar.f90   $(SOIL_DIR)/mnet2arps.f90 \
             $(SOIL_DIR)/read_obs.f90

SOILALL = $(SOILREADME) $(SOILMKR) $(SOILSRC) $(SOILINPUT)

#-----------------------------------------------------------------------
#
# List of ARPSMP files in sub-directory arps_mp
#
#-----------------------------------------------------------------------

MPREADME = $(MP_DIR)/README.arps_mp

MPMKR    = $(MP_DIR)/Makefile

MPINPUT  = $(INPUT_DIR)/splitany.input

MPSRC    = $(MP_DIR)/fjoindumps.f90    $(MP_DIR)/fsplitdump.f90    \
           $(MP_DIR)/fsplitexbc.f90    $(MP_DIR)/fsplitsoil.f90    \
           $(MP_DIR)/fsplitsoilini.f90 $(MP_DIR)/fsplitterrain.f90 \
           $(MP_DIR)/joinfile.f90      $(MP_DIR)/joinfiles.f90     \
           $(MP_DIR)/splitfiles.f90    $(MP_DIR)/fsplitrestart.f90 \
           $(MP_DIR)/fjoinbin2hdf.f90  $(MP_DIR)/joinbin2hdf.f90

MPHDFSRC = $(MP_DIR)/fjoin_hdf.f90     $(MP_DIR)/fsplit_hdf.f90 \
           $(MP_DIR)/joinhdf.f90       $(MP_DIR)/splithdf.f90   \
           $(MP_DIR)/fsplithdf.f90     $(MP_DIR)/fsplitnet.f90  \
           $(MP_DIR)/splitany.f90                               \
           $(MP_DIR)/nohdf.f90         $(MP_DIR)/nonet.f90

MPALL = $(MPREADME) $(MPMKR) $(MPSRC) $(MPHDFSRC) $(MPINPUT)

#-----------------------------------------------------------------------
#
# List of ARPSTERN files in sub-directory terrain
#
#-----------------------------------------------------------------------

TERNREADME = $(TERN_DIR)/README.terrain

TERNMKR    = $(TERN_DIR)/Makefile
TERNINPUT  = $(INPUT_DIR)/arpstern.input

TERNSRC    = $(TERN_DIR)/arpstern.f90 $(TERN_DIR)/dir1deg.f90    \
             $(TERN_DIR)/dir5min.f90  $(TERN_DIR)/dir30sec.f90

TERNINC    =

TERNMISC   = $(TERN_DIR)/load1deg.file  \
             $(TERN_DIR)/load5min.file  $(TERN_DIR)/load30sec.file

TERNALL = $(TERNREADME) $(TERNMKR) $(TERNSRC) $(TERNINC) $(TERNMISC) \
          $(TERNINPUT)

#-----------------------------------------------------------------------
#
# List of ARPSTRN files in sub-directory terrain
#
#-----------------------------------------------------------------------

TRNREADME = $(TRN_DIR)/README.arpstrn

TRNMKR    = $(TRN_DIR)/Makefile
TRNINPUT  = $(INPUT_DIR)/arpstrn.input

TRNSRC    = $(TRN_DIR)/arpstrn.f90 $(TRN_DIR)/nozxplot.f90

TRNMISC   = $(TRN_DIR)/usgs_dem.index

TRNALL = $(TRNREADME) $(TRNMKR) $(TRNSRC) $(TRNINC) $(TRNMISC) \
         $(TRNINPUT) \
         $(MERGETRNSRC) $(MERGETRNINPUT)

MERGETRNSRC = $(TRN_DIR)/mergetrn.f90
MERGETRNINPUT = $(INPUT_DIR)/mergetrn.input

#-----------------------------------------------------------------------
#
# List of EXT2ARPS files in sub-directory e2a
#
#-----------------------------------------------------------------------

E2AREADME = $(E2A_DIR)/README.ext2arps

E2AMKR    = $(E2A_DIR)/Makefile
E2AINPUT  = $(ARPSINPUT)

E2ASRC    = $(E2A_DIR)/ext2arps.f90   $(E2A_DIR)/extlib.f90      \
            $(E2A_DIR)/getextd3d.f90  $(E2A_DIR)/getgempak3d.f90 \
            $(E2A_DIR)/getlaps3d.f90  $(E2A_DIR)/gribio_c.c      \
            $(E2A_DIR)/griblib3d.f90  $(E2A_DIR)/nmcdecode.f90   \
            $(E2A_DIR)/nogempak3d.f90 $(E2A_DIR)/nolaps3d.f90    \
            $(E2A_DIR)/rdnmcgrb3d.f90 $(E2A_DIR)/getcoamps.f90   \
            $(E2A_DIR)/getwrfdata.f90 $(E2A_DIR)/nowrfnet.f90

E2AINC    = $(INCL_DIR)/gribcst.inc $(INCL_DIR)/lapsparms.cmn \
            $(INCL_DIR)/ext2arps.inc

E2AALL    = $(E2AREADME) $(E2AMKR) $(E2ASRC) $(E2AINC) $(E2AINPUT)

A2GINPUT  = $(INPUT_DIR)/arps2gem.input
A2GSRC    = $(E2A_DIR)/arps2gem.f90 $(E2A_DIR)/v2dint.f90
A2GINC    = $(INCL_DIR)/arps2gem.inc
A2GALL    = $(A2GINPUT) $(A2GSRC) $(A2GINC)

A2NINPUT  = $(INPUT_DIR)/arps2ncdf.input
A2NSRC    = $(E2A_DIR)/arps2ncdf.f90 $(E2A_DIR)/v2dint.f90 $(E2A_DIR)/arps2wdssii.f90
A2NALL    = $(A2NINPUT) $(A2NSRC)

A2EINPUT  = $(INPUT_DIR)/arps2eta212.input
A2ESRC    = $(E2A_DIR)/arps2eta212.f90  $(E2A_DIR)/v2dint.f90     \
            $(E2A_DIR)/extlib.f90
A2EALL    = $(A2EINPUT) $(A2ESRC)

A2WINPUT  = $(INPUT_DIR)/arps2wrf.input
A2WSRC    = $(A2W_DIR)/arps2wrf.f90  $(A2W_DIR)/interplib.f90     \
            $(A2W_DIR)/wrf_iolib.f90 $(A2W_DIR)/wrf_mpsubs.f90    \
            $(A2W_DIR)/wrf_ioncd.f90 $(A2W_DIR)/wrf_nompsubs.f90  \
            $(A2W_DIR)/wrf_iophdf5.f90 $(A2W_DIR)/wrf_ionophdf5.f90 \
            $(A2W_DIR)/dump_wrf_input.f90 $(A2W_DIR)/wrf_subs.F   \
            $(A2W_DIR)/dump_wrf_bdy.f90                           \
            $(A2W_DIR)/dump_wrf_static.f90                        \
            $(A2W_DIR)/module_wrf_metadata.f90                    \
            $(A2W_DIR)/readnamelist.f90

WSSRC     = $(A2W_DIR)/wrfstatic.f90 $(A2W_DIR)/wrfstaticlib.f90  \
            $(A2W_DIR)/wrfsi_subs.f  $(A2W_DIR)/wrfsi_cio.c       \
            $(A2W_DIR)/dump_wrf_static.f90

A2WALL    = $(A2WINPUT) $(A2WSRC) $(A2W_DIR)/Makefile $(WSSRC)

AVNINPUT  = $(INPUT_DIR)/extract_avn.input
AVNSRC    = $(AVN_DIR)/extract_avn.f90
AVNINC    =
AVNALL    = $(AVNINPUT) $(AVNSRC) $(AVNINC)

W2AINPUT  = $(INPUT_DIR)/wrf2arps.input  $(INPUT_DIR)/joinwrf.input
W2ASRC    = $(W2A_DIR)/wrf2arps.f90      $(W2A_DIR)/getwrfd3d.f90  \
            $(W2A_DIR)/wrfintrplib.f90   $(W2A_DIR)/wrf_subs.F     \
            $(W2A_DIR)/wrf_nompsubs.f90  $(W2A_DIR)/wrf_mpsubs.f90 \
            $(W2A_DIR)/wrf_ionophdf5.f90 $(W2A_DIR)/wrf_iophdf5.f90 \
            $(W2A_DIR)/wrf_ioncd.f90     $(W2A_DIR)/joinwrf.f90    \
            $(W2A_DIR)/get_wrf_one_file.f90 $(W2A_DIR)/fjoinwrf.f90 \
            $(W2A_DIR)/get_wrf_multi_files.f90 $(W2A_DIR)/module_wrf2arps_post.f90

WEXTSNDSRC = $(W2A_DIR)/wrfextsnd.f90

W2AALL    = $(W2AINPUT) $(W2ASRC) $(WEXTSNDSRC) $(W2A_DIR)/Makefile

#----------------------------------------------------------------------
#
# List of files of arps4wrf and nmm2arps
#
#-----------------------------------------------------------------------

A4WINPUT  = $(INPUT_DIR)/arps4wrf.input
A4WSRC    = $(A4W_DIR)/Makefile            $(A4W_DIR)/arps4wrf.f90  \
            $(A4W_DIR)/process_domain.f90  $(A4W_DIR)/output_domain.f90 \
            $(A4W_DIR)/interplib.f90       $(A4W_DIR)/interpmpi.F   \
            $(A4W_DIR)/module_arpsgrid.f90 $(A4W_DIR)/module_metgrid.f90 \
            $(A4W_DIR)/module_geogrid.f90  $(A4W_DIR)/module_commontypes.f90 \
            $(A4W_DIR)/module_constants.f90 $(A4W_DIR)/module_namelist.f90   \
            $(A4W_DIR)/wrf_maputils_module.f90 $(A4W_DIR)/wrf_mapprojs_module.f90 \
            $(A4W_DIR)/wrf_parallel_module.F $(A4W_DIR)/queue_module.f90 \
            $(A4W_DIR)/read_arpsmp.f90     $(A4W_DIR)/static_input_module.f90

A4WALL    = $(A4WINPUT) $(A4WSRC)

N2AINPUT  = $(INPUT_DIR)/nmm2arps.input
N2ASRC    = $(N2A_DIR)/Makefile                 $(N2A_DIR)/nmm2arps.f90   \
            $(N2A_DIR)/process_file.f90         $(N2A_DIR)/interplib.f90  \
            $(N2A_DIR)/module_arpsgrid.f90      $(N2A_DIR)/module_nmmgrid.f90  \
            $(N2A_DIR)/module_interpolation.f90 $(N2A_DIR)/module_namelist.f90 \
            $(N2A_DIR)/module_nmm_input.f90

N2AALL    = $(N2AINPUT) $(N2ASRC)

#----------------------------------------------------------------------
#
# List of files of arps2coamps
#
#-----------------------------------------------------------------------

A2CINPUT  = $(INPUT_DIR)/arps2coamps.input
A2CSRC    = $(A2C_DIR)/Makefile            $(A2C_DIR)/arps2coamps.f90   \
            $(A2C_DIR)/process_domain.f90  $(A2C_DIR)/output_domain.f90 \
            $(A2C_DIR)/interplib.f90       $(A2C_DIR)/module_interpolation.f90 \
            $(A2C_DIR)/module_arpsgrid.f90 $(A2C_DIR)/module_coampsgrid.f90 \
            $(A2C_DIR)/module_constants.f90 $(A2C_DIR)/module_namelist.f90  \
            $(A2C_DIR)/coamps_mapprojs_module.f90  $(A2C_DIR)/xgetgg.c      \
            $(A2C_DIR)/coamps_parallel_module.F $(A2C_DIR)/module_topo.F      \
            $(A2C_DIR)/read_arpsmp.f90     $(A2C_DIR)/static_input_module.f90 \

A2CALL    = $(A2CINPUT) $(A2CSRC)

#----------------------------------------------------------------------
#
# List of files of coamps2arps
#
#-----------------------------------------------------------------------

C2AINPUT  = $(INPUT_DIR)/coamps2arps.input
C2ASRC    = $(C2A_DIR)/Makefile            $(C2A_DIR)/coamps2arps.f90   \
            $(C2A_DIR)/process_file.f90    $(C2A_DIR)/interplib.f90 \
            $(C2A_DIR)/module_coampsgrid.f90 $(A2C_DIR)/module_namelist.f90 \
            $(C2A_DIR)/module_interpolation.f90
C2AIMPORT = $(N2A_DIR)/module_arpsgrid.f90 $(A2C_DIR)/coamps_constants.f90 \
            $(A2C_DIR)/coamps_mapprojs_module.f90                       \
            $(A2C_DIR)/coamps_parallel_module.F

C2AALL    = $(C2AINPUT) $(C2ASRC) $(C2AIMPORT)

#-----------------------------------------------------------------------
#
# List of 88D2ARPS files in sub-directory 88d2arps
#
#-----------------------------------------------------------------------

88DREADME = $(88D_DIR)/README.88d2arps

88DMKR    = $(88D_DIR)/Makefile

88DSRC    = $(88D_DIR)/f88d2arps.f90  $(88D_DIR)/fakeio.c               \
            $(88D_DIR)/fakerad.f90    $(88D_DIR)/remaptilt.f90          \
            $(88D_DIR)/remaplib.f90                                     \
            $(88D_DIR)/pltradcol.f90  $(88D_DIR)/pltradarscan.f90       \
            $(88D_DIR)/ncrad2arps.f90 $(88D_DIR)/rdtiltcdf.f90          \
            $(88D_DIR)/nids2arps.f90  $(88D_DIR)/rasterlib.f90          \
            $(88D_DIR)/nidsdecoders.f90                                 \
            $(88D_DIR)/radarinfo.c    $(88D_DIR)/adjustreflib.f90       \
            $(88D_DIR)/apdetect.f90   $(88D_DIR)/unfldlib.f90           \
            $(88D_DIR)/unfold.f90     $(88D_DIR)/wrttilts.f90           \
            $(88D_DIR)/casa2arpsppi.f90 $(88D_DIR)/rdtilt88d.f90

88DINC    = $(INCL_DIR)/remap.inc    $(INCL_DIR)/remapsum.inc           \
            $(INCL_DIR)/remaptab.inc $(INCL_DIR)/remapcst.inc           \
            $(INCL_DIR)/ulog.h       $(INCL_DIR)/alarm.h                \
            $(INCL_DIR)/cfortran.h   $(INCL_DIR)/udposix.h              \
            $(INCL_DIR)/nidscst.inc

88DINPUT  = $(ARPSINPUT) $(INPUT_DIR)/fakerad.input                     \
            $(INPUT_DIR)/pltradcol.input                                \
            $(INPUT_DIR)/radremap_NIDS.input                            \
            $(INPUT_DIR)/radremap_KTLX.input                            \
            $(INPUT_DIR)/radremap_KSAO.input                            \
            $(INPUT_DIR)/radarplt.input                                 \
            $(INPUT_DIR)/casa2arps.input  $(INPUT_DIR)/ppiplt.input

LIBA2IOSRC= $(88D_DIR)/a2io/*.c $(88D_DIR)/a2io/*.h $(88D_DIR)/a2io/Makefile

SOLOSRC = $(88D_DIR)/soloio/*.c $(88D_DIR)/soloio/*.h $(88D_DIR)/soloio/Makefile

88DDAT    = $(PLTR_DIR)/spcounty.mapdata

88DALL    = $(88DREADME)  $(88DMKR) $(88DSRC) $(88DINC) $(88DINPUT)     \
            $(LIBA2IOSRC) $(SOLOSRC) $(88DDAT)

#-----------------------------------------------------------------------
#
# List of RDREMUL files in sub-directory radaremul
#
#-----------------------------------------------------------------------

RDREMULSRC = $(RDR_EMUL)/convert2radar.f90 $(RDR_EMUL)/dualpara.f90   \
             $(RDR_EMUL)/global_module.f90 $(RDR_EMUL)/tmatrix.f90

#-----------------------------------------------------------------------
#
# List of RADMOSAIC files in sub-directory adas
#
#-----------------------------------------------------------------------

RADMOSSRC   = $(ADAS_DIR)/radmosaic.f90
RADMOSINPUT = $(INPUT_DIR)/radmosaic.input
RADMOSALL   = $(ADASMKR) $(RADMOSSRC) $(ARPSINC) $(RADMOSINPUT)

#-----------------------------------------------------------------------
#
# List of ARPS2RAD files in sub-directory 88d2arps
#
#-----------------------------------------------------------------------

88DALL    = $(88DREADME)  $(88DMKR) $(88DSRC) $(88DINC) $(88DINPUT)     \
            $(LIBA2IOSRC) $(88DDAT)
ARPS2RADSRC   = $(88D_DIR)/arps2rad.f90
ARPS2RADINPUT = $(INPUT_DIR)/arps2rad.input
ARPS2RADALL = $(88DMKR) $(ARPS2RADSRC) $(ARPSINC) $(ARPS2RADINPUT)

#-----------------------------------------------------------------------
#
# List of RADAREMUL  files in sub-directory 88d2arps
#
#-----------------------------------------------------------------------

RADEMULSRC   = $(88D_DIR)/radaremul.f90 $(88D_DIR)/dualpara.f90  \
               $(88D_DIR)/dualtmat.f90  $(88D_DIR)/duallib.f90   \
               $(88D_DIR)/wttiltcdf.f90 $(88D_DIR)/dummycdf.f90

#RADEMULINPUT = $(INPUT_DIR)/radaremul.input

RADEMULALL = $(88DMKR) $(RADEMULSRC) $(ARPSINC) $(RADEMULINPUT)

#-----------------------------------------------------------------------
#
# List of RADAREMUL  files in sub-directory radaremul
#
#-----------------------------------------------------------------------

PRDLIBSRC  = $(RDREMULDIR)/convert2radar.f90      \
             $(RDREMULDIR)/dualpara.f90           \
             $(RDREMULDIR)/global_module.f90      \

PRDLIBINC  = $(INCL_DIR)/globcst.inc    $(INCL_DIR)/phycst.inc

PRDLIBALL  = $(PRDLIBSRC) $(PRDLIBINC)

#-----------------------------------------------------------------------
#
# List of ATTEN files in sub-directory 88d2arps
#
#-----------------------------------------------------------------------

ATTENSRC   = $(88D_DIR)/atten.f90     \
             $(88D_DIR)/wttiltcdf.f90 $(88D_DIR)/dummycdf.f90

#ATTENINPUT = $(INPUT_DIR)/atten.input

ATTENALL = $(88DMKR) $(ATTENSRC) $(ARPSINC) $(ATTENINPUT)

#-----------------------------------------------------------------------
#
# List of RADBIN2CDF  files in sub-directory 88d2arps
#
#-----------------------------------------------------------------------

RADB2CDFSRC  = $(88D_DIR)/radbin2cdf.f90 $(88D_DIR)/wttiltcdf.f90

RADB2CDFALL = $(88DMKR) $(RADB2CDFSRC)

#
#-----------------------------------------------------------------------
#
# List of RADSECTOR  files in sub-directory 88d2arps
#
#-----------------------------------------------------------------------

RADSECTSRC  = $(88D_DIR)/radsector.f90 $(88D_DIR)/rdtiltcdf.f90  \
              $(88D_DIR)/wttiltcdf.f90

RADSECTALL = $(88DMKR) $(RADSECTSRC)

#-----------------------------------------------------------------------
#
# List of MCI2ARPS files in sub-directory mci2arps
#
#-----------------------------------------------------------------------

MCIREADME = $(MCI_DIR)/README.mci2arps

MCIMKR    = $(MCI_DIR)/Makefile

MCISRC    = $(MCI_DIR)/mci2arps.c   $(MCI_DIR)/calib.c                  \
            $(MCI_DIR)/gldtonat.c   $(MCI_DIR)/inisatarps.f90           \
            $(MCI_DIR)/maptran.c    $(MCI_DIR)/sat_read.c               \
            $(MCI_DIR)/solcor.f90   $(MCI_DIR)/sunfuncs.f90             \
            $(MCI_DIR)/wtsatfld.f90 $(MCI_DIR)/pltsatfld.f90            \
            $(MCI_DIR)/pltmap.f90   $(MCI_DIR)/sat_tran.c               \
            $(MCI_DIR)/palgrey.hdf  $(MCI_DIR)/coldfilt.f90             \
            $(MCI_DIR)/sat2arps.c   $(MCI_DIR)/mergesat.c               \
            $(MCI_DIR)/sathdf52arps.f90

MCIINC    = $(INCL_DIR)/mc_area.h  $(INCL_DIR)/wxp.h                    \
            $(INCL_DIR)/wxpresrc.h $(INCL_DIR)/wxpgks.h

MCIINPUT  = $(ARPSINPUT) $(INPUT_DIR)/pltsatfld.input

MCIALL    = $(MCIREADME) $(MCIMKR) $(MCISRC) $(MCIINC) $(MCIINPUT)

#-----------------------------------------------------------------------
#
# List of MOSAIC2ARPS files in sub-directory mosaic2arps
#
#-----------------------------------------------------------------------

MOSAICMKR    = $(MOSAIC_DIR)/Makefile $(MOSAIC_DIR)/README

MOSAICSRC    = $(MOSAIC_DIR)/Zmosaic2arps.f90   $(MOSAIC_DIR)/Zmosaiclib.f90

MOSAICINPUT  = $(INPUT_DIR)/Zmosaic2arps.input

MOSAICALL    = $(MOSAICMKR) $(MOSAICSRC) $(MOSAICINPUT)

#-----------------------------------------------------------------------
#
# List of ARPSEXTSND files in sub-directory extsnd
#
#-----------------------------------------------------------------------

EXTSNDREADME = $(EXTSND_DIR)/README.arpsextsnd
EXTSNDMKR    = $(EXTSND_DIR)/Makefile

EXTSNDINPUT  = $(INPUT_DIR)/arpsextsnd.input
EXTSNDSRC    = $(EXTSND_DIR)/arpsextsnd.f90 $(EXTSND_DIR)/arpsextsndderive.f90

EXTSNDALL = $(EXTSNDREADME) $(EXTSNDMKR) $(EXTSNDSRC) $(EXTSNDINPUT)

#-----------------------------------------------------------------------
#
# List of SKEWT files in sub-directory src/skewt
#
#-----------------------------------------------------------------------

SKEWTREADME = $(SKEWT_DIR)/README.skewt
SKEWTMKR    = $(SKEWT_DIR)/Makefile

SKEWTSRC    = $(SKEWT_DIR)/skewt.F       $(SKEWT_DIR)/barb.f         \
              $(SKEWT_DIR)/c.miscsubs2.f $(SKEWT_DIR)/denmwind.f     \
              $(SKEWT_DIR)/helcont.f     $(SKEWT_DIR)/hodograph.f    \
              $(SKEWT_DIR)/indices.f     $(SKEWT_DIR)/misc_kb.f      \
              $(SKEWT_DIR)/myzxsubs.f    $(SKEWT_DIR)/parseinput.f90 \
              $(SKEWT_DIR)/preciptype.f  $(SKEWT_DIR)/readsound.f    \
              $(SKEWT_DIR)/sfcupd.f      $(SKEWT_DIR)/thermosubs.f

SKEWTINC    = $(INCL_DIR)/thermo.consts   $(INCL_DIR)/thermo.stfunc
SKEWTDAT    = $(PLTD_DIR)/skewt.pltcbar   $(DATA_DIR)/data.test/SHV.txt \
              $(DATA_DIR)/data.test/AMA16may9212z.txt

SKEWTALL    = $(SKEWTREADME) $(SKEWTMKR) $(SKEWTSRC) $(SKEWTINC) $(SKEWTDAT)

#-----------------------------------------------------------------------
#
# List of WTRETCOL files in sub-directory wtretcol
#
#-----------------------------------------------------------------------

WTRETCOLRDME = $(WTRETCOL_DIR)/README.wtretcol

WTRETCOLMKR  = $(WTRETCOL_DIR)/Makefile
WTRETCOLINPUT=

WTRETCOLSRC  = $(WTRETCOL_DIR)/wtretcol.f90

WTRETCOLALL = $(WTRETCOLRDME) $(WTRETCOLMKR) $(WTRETCOLSRC) \
              $(WTRETCOLINPUT)

#-----------------------------------------------------------------------
#
# List of ZXPLOT files
#
#-----------------------------------------------------------------------

ZXPLOTMKFILE = $(ZXPLOT_DIR)/Makefile

ZXPLOTRDME   = $(ZXPLOT_DIR)/README

ZXPLOTSRC = $(ZXPLOT_DIR)/zxplot3.f                                     \
            $(ZXPLOT_DIR)/xpost3.f     $(ZXPLOT_DIR)/zxpslib3.f         \
            $(ZXPLOT_DIR)/xncar3.f     $(ZXPLOT_DIR)/zxnglib3.f

ZXPLOTALL = $(ZXPLOTSRC) $(ZXPLOTMKFILE) $(ZXPLOTRDME)

#-----------------------------------------------------------------------
#
# List of external PHDF5 package (From WRFV2.0.3)
#
#-----------------------------------------------------------------------

PHDF5MKFILE = $(WRF_PHDF5_DIR)/Makefile

PHDF5RDME   = $(WRF_PHDF5_DIR)/README.phdf5

PHDF5INC    = $(WRF_EXT_DIR)/inc/wrf_io_flags.h                         \
              $(WRF_EXT_DIR)/inc/wrf_status_codes.h

PHDF5SRC  = $(WRF_PHDF5_DIR)/wrf-phdf5.F90                              \
            $(WRF_PHDF5_DIR)/wrf-phdf5attr.F90                          \
            $(WRF_PHDF5_DIR)/wrf-phdf5bdy.F90                           \
            $(WRF_PHDF5_DIR)/wrf-phdf5support.F90

LIBPHDF5ALL = $(PHDF5INC) $(PHDF5SRC) $(PHDF5MKFILE) $(PHDF5RDME)

#-----------------------------------------------------------------------
#
# List of external WRF binary package (From WRFV2.0.3)
#
#-----------------------------------------------------------------------

WRFINTMKFILE = $(WRF_INT_DIR)/Makefile

WRFINTRDME = $(WRF_INT_DIR)/README.int

WRFINTINC  = $(WRF_EXT_DIR)/inc/intio_tags.h                        \
             $(WRF_EXT_DIR)/inc/wrf_io_flags.h                      \
             $(WRF_EXT_DIR)/inc/wrf_status_codes.h

WRFINTSRC  = $(WRF_INT_DIR)/io_int.F                              \
             $(WRF_INT_DIR)/module_internal_header_util.F

LIBWRFINTALL = $(WRFINTINC) $(WRFINTSRC) $(WRFINTMKFILE) $(WRFINTRDME)

#-----------------------------------------------------------------------
#
# List of external WRF netcdf package (From WRFV3.1.1)
#
#-----------------------------------------------------------------------

WRFNETMKFILE = $(WRF_NET_DIR)/Makefile

WRFNETRDME = $(WRF_NET_DIR)/README

WRFNETSRC  = $(WRF_NET_DIR)/wrf_io.F $(WRF_NET_DIR)/module_wrfsi_static.F \
             $(WRF_NET_DIR)/field_routines.F $(WRF_NET_DIR)/bitwise_operators.c \
             $(WRF_NET_DIR)/ext_ncd_get_dom_ti.code  \
             $(WRF_NET_DIR)/ext_ncd_get_var_td.code  \
             $(WRF_NET_DIR)/ext_ncd_get_var_ti.code  \
             $(WRF_NET_DIR)/ext_ncd_put_dom_ti.code  \
             $(WRF_NET_DIR)/ext_ncd_put_var_td.code  \
             $(WRF_NET_DIR)/ext_ncd_put_var_ti.code  \
             $(WRF_NET_DIR)/transpose.code

LIBWRFNETALL = $(WRFNETSRC) $(WRFNETMKFILE) $(WRFNETRDME)

#-----------------------------------------------------------------------
#
# List of external GRIB 2 package (Version 1.0.7)
#
#-----------------------------------------------------------------------

LIBG2RDME = $(G2_EXT_DIR)/Makefile $(G2_EXT_DIR)/README $(G2_EXT_DIR)/makefile.orig

LIBG2SRC = $(G2_EXT_DIR)/addfield.F  $(G2_EXT_DIR)/gf_unpack7.F \
           $(G2_EXT_DIR)/addgrid.f   $(G2_EXT_DIR)/addlocal.f   \
           $(G2_EXT_DIR)/cmplxpack.f $(G2_EXT_DIR)/compack.f \
           $(G2_EXT_DIR)/comunpack.f $(G2_EXT_DIR)/dec_jpeg2000.c \
           $(G2_EXT_DIR)/dec_png.c   $(G2_EXT_DIR)/drstemplates.f \
           $(G2_EXT_DIR)/enc_jpeg2000.c $(G2_EXT_DIR)/enc_png.c \
           $(G2_EXT_DIR)/g2grids.f   $(G2_EXT_DIR)/gb_info.f \
           $(G2_EXT_DIR)/gbytesc.f   $(G2_EXT_DIR)/gdt2gds.f \
           $(G2_EXT_DIR)/getdim.f    $(G2_EXT_DIR)/getfield.f \
           $(G2_EXT_DIR)/getg2i.f    $(G2_EXT_DIR)/getg2ir.f \
           $(G2_EXT_DIR)/getgb2.f    $(G2_EXT_DIR)/getgb2l.f \
           $(G2_EXT_DIR)/getgb2p.f   $(G2_EXT_DIR)/getgb2r.f \
           $(G2_EXT_DIR)/getgb2rp.f  $(G2_EXT_DIR)/getgb2s.f \
           $(G2_EXT_DIR)/getidx.f    $(G2_EXT_DIR)/getlocal.f \
           $(G2_EXT_DIR)/getpoly.f   $(G2_EXT_DIR)/gettemplates.f \
           $(G2_EXT_DIR)/gf_free.f   $(G2_EXT_DIR)/gf_getfld.f \
           $(G2_EXT_DIR)/gf_unpack1.f $(G2_EXT_DIR)/gf_unpack2.f \
           $(G2_EXT_DIR)/gf_unpack3.f $(G2_EXT_DIR)/gf_unpack4.f \
           $(G2_EXT_DIR)/gf_unpack5.f $(G2_EXT_DIR)/gf_unpack6.f \
           $(G2_EXT_DIR)/gf_unpack7.F $(G2_EXT_DIR)/grib2.doc \
           $(G2_EXT_DIR)/gribcreate.f $(G2_EXT_DIR)/gribend.f \
           $(G2_EXT_DIR)/gribinfo.f   $(G2_EXT_DIR)/gribmod.f \
           $(G2_EXT_DIR)/gridtemplates.f $(G2_EXT_DIR)/ixgb2.f \
           $(G2_EXT_DIR)/jpcpack.f    $(G2_EXT_DIR)/jpcunpack.f \
           $(G2_EXT_DIR)/misspack.f   $(G2_EXT_DIR)/mkieee.f \
           $(G2_EXT_DIR)/mova2i.c     $(G2_EXT_DIR)/pack_gp.f \
           $(G2_EXT_DIR)/params.f     $(G2_EXT_DIR)/pdstemplates.f \
           $(G2_EXT_DIR)/pngpack.f $(G2_EXT_DIR)/pngunpack.f \
           $(G2_EXT_DIR)/putgb2.f  $(G2_EXT_DIR)/rdieee.f \
           $(G2_EXT_DIR)/realloc.f $(G2_EXT_DIR)/reduce.f \
           $(G2_EXT_DIR)/simpack.f $(G2_EXT_DIR)/simunpack.f \
           $(G2_EXT_DIR)/skgb.f    $(G2_EXT_DIR)/specpack.f \
           $(G2_EXT_DIR)/specunpack.f

LIBG2ALL = $(LIBG2RDME) $(LIBG2SRC)

#-----------------------------------------------------------------------
#
# List of f77 to f90 converter files
#
#-----------------------------------------------------------------------

F2F90SRC = $(F2F90DIR)/convert_src_to_f90.f90 \
           $(F2F90DIR)/convert_inc_to_mod.f90 \
           $(F2F90DIR)/NOTES $(F2F90DIR)/F90ARPS.html

#-----------------------------------------------------------------------
#
# List of machine-dependent object codes
#
#-----------------------------------------------------------------------

MACHLIBOBJ = $(MACHDEPLIB)
NCARGOBJ   = $(NCARGDEP)
ZXPLOTOBJ  = $(ZXPLOTDEP)
ASSIMOBJ   = $(ASSIMOPT)

MPIOBJ  = $(MPIDEP)

HDFOBJ  = $(HDFDEP)
PAKOBJ  = $(PAKDEP)
SVIOBJ  = $(SVIDEP)
NETOBJ  = $(NETDEP)
V5DOBJ  = $(V5DDEP)
GEMOBJ  = $(GEMDEP)
PHDF5OBJ = $(PHDF5DEP)
GRIB2OBJ = $(GRIB2DEP)
CRTMOBJ = $(CRTMDEP)
CITMOBJ = $(CITMDEP)

#-----------------------------------------------------------------------
#
# Set Default
#
#-----------------------------------------------------------------------

default:
	echo; echo Please use makearps script to do make.; echo

#-----------------------------------------------------------------------
#
# Make all executables except arps_mp
#
#-----------------------------------------------------------------------

#
# Not working executables:
#
#    $(AGREXE)
#
all: $(ADASEXE)     $(ARPSEXE)     $(ARPSCVTEXE)   $(ARPSREADEXE)      \
     $(ARPSSOILEXE) $(ARPSSFCEXE)  $(INTRPEXE)     $(INTRP_LSEXE)      \
     $(PRTEXE)      $(DIFFEXE)     $(EXTSNDEXE)    $(E2AEXE)           \
     $(TERNEXE)     $(DIR1DEGEXE)   $(DIR5MINEXE)                      \
     $(DIR30SECEXE) $(WTRETCOLEXE) $(SPLITEXE)     $(JOINSEXE)         \
     $(JOINEXE)     $(ENSCVEXE)    $(ENSBCEXE)     $(ENSICEXE)         \
     $(EPOSTEXE)    $(ENSANAEXE)   ${ENSCALEXE}    $(RAINDIFFEXE)  $(ASSIMEXE) \
     $(ARPSVERIF)   $(TRNEXE)      $(VERIFEXE)     $(TRAJCEXE)
#     $(ENKFEXE)     $(ENKFICEXE)   $(RNDPRTEXE)   $(RADARINTRPEXE)     \
#     $(INITENSMBLEXE)

#
# All mpi
#
all.mpi: $(ARPSMPEXE) $(E2AMPEXE) $(VERIFMPEXE) $(A2WMPEXE) $(W2AMPEXE) \
         $(PLTNCARMPEXE) $(PLTPOSTMPEXE) $(EXTSNDMPEXE) $(SPLITMPEXE)  \
         $(ENSICMPEXE) $(EPOSTMPEXE) $(ADASMPEXE) $(3DVARMPEXT) \
         $(88D2AMPEXE)

#
# Depends on ZXPLOT library.
#
all.plot: $(PLTMAXEXE) $(PLTNCAREXE) $(PLTPOSTEXE)   $(PLTGRIDEXE)     \
          $(SKEWNEXE)  $(SKEWPEXE)   $(PLTPOSTMPEXE) $(PLTNCARMPEXE)

#
# Depends on NetCDF 3.0.
#
all.netcdf: $(A2NEXE) $(A2WIIEXE) $(A2WEXE) $(W2AEXE) $(WSEXE) \
            $(A2WMPEXE) $(W2AMPEXE)

#
# Depends on Gempak.
#
all.gempak: $(ARPSCVTOBSEXE) $(A2GEXE) $(G2AEXE) $(EPOSTEXE) \
            $(EPOSTMPEXE) $(B2GEXE)

#-----------------------------------------------------------------------
#
# Compile and link ARPS model executable, arps
#
#-----------------------------------------------------------------------

$(ARPSEXE):
	cd $(ARPSDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'CPP=$(CPP)'                   \
                            'CPPFLAGS=$(CPPFLAGS)'         \
                            'CFLAGS=$(CFLAGS)'             \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'ICEFLAG=$(ICEFLAG)'           \
                            'RADFLAG=$(RADFLAG)'           \
                            'THMFLAG=$(THMFLAG)'           \
                            'LIBS=$(LIBS)'                 \
                            'MACHLIBOBJ=$$($(MACHLIBOBJ))' \
                            'MPIOBJ=$$($(MPIOBJ))'         \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'PAKOBJ=$$($(PAKOBJ))'         \
                            'SVIOBJ=$$($(SVIOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            'V5DOBJ=$$($(V5DOBJ))'         \
                            $(ARPSEXE)

#-----------------------------------------------------------------------
#
# Compile and link ARPS MP model executable, arps_mpi
#
#-----------------------------------------------------------------------

$(ARPSMPEXE):
	cd $(ARPSDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'CPP=$(CPP)'                   \
                            'CPPFLAGS=$(CPPFLAGS)'         \
                            'CFLAGS=$(CFLAGS)'             \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'ICEFLAG=$(ICEFLAG)'           \
                            'RADFLAG=$(RADFLAG)'           \
                            'THMFLAG=$(THMFLAG)'           \
                            'LIBS=$(LIBS)'                 \
                            'MACHLIBOBJ=$$($(MACHLIBOBJ))' \
                            'MPIOBJ=$$($(MPIOBJ))'         \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'PAKOBJ=$$($(PAKOBJ))'         \
                            'SVIOBJ=$$($(SVIOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            'V5DOBJ=$$($(V5DOBJ))'         \
                            $(ARPSMPEXE)

#-----------------------------------------------------------------------
#
# Make arpssolver.a, libarps.a, and shared.a
#
#-----------------------------------------------------------------------

$(ARPSSOLVER):
	cd $(ARPSDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'CFLAGS=$(CFLAGS)'             \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'ICEFLAG=$(ICEFLAG)'           \
                            'RADFLAG=$(RADFLAG)'           \
                            'THMFLAG=$(THMFLAG)'           \
                            'LIBS=$(LIBS)'                 \
                            'MACHLIBOBJ=$$($(MACHLIBOBJ))' \
                            'ASSIMOPT=$$($(ASSIMOBJ))'     \
                            'MPIOBJ=$$($(MPIOBJ))'         \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'PAKOBJ=$$($(PAKOBJ))'         \
                            'SVIOBJ=$$($(SVIOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            'V5DOBJ=$$($(V5DOBJ))'         \
                            $(ARPSSOLVER)

$(LIBARPS):
	cd $(ARPSDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'CPP=$(CPP)'                   \
                            'CPPFLAGS=$(CPPFLAGS)'         \
                            'CFLAGS=$(CFLAGS)'             \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'LIBS=$(LIBS)'                 \
                            'MACHLIBOBJ=$$($(MACHLIBOBJ))' \
                            'NCARGOBJ=$$($(NCARGOBJ))'     \
                            'MPIOBJ=$$($(MPIOBJ))'         \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'PAKOBJ=$$($(PAKOBJ))'         \
                            'SVIOBJ=$$($(SVIOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            'V5DOBJ=$$($(V5DOBJ))'         \
                            'ARFLAG=$(ARFLAG)'             \
                            $(LIBARPS)

$(LIBZXPOST):
	cd $(ZXPLOTDIR); make $(MAKEOPT)                   \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FIXFLAGS=$(FIXFLAGS)'         \
                            'CFLAGS=$(CFLAGS)'             \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'LIBS=$(LIBS)'                 \
                            'ARFLAG=$(ARFLAG)'             \
                            $(LIBZXPOST)

$(LIBZXNCAR):
	cd $(ZXPLOTDIR); make $(MAKEOPT)                   \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FIXFLAGS=$(FIXFLAGS)'         \
                            'CFLAGS=$(CFLAGS)'             \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'LIBS=$(LIBS)'                 \
                            'ARFLAG=$(ARFLAG)'             \
                            $(LIBZXNCAR)

$(LIBRADTN):
	cd $(ARPSDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'CFLAGS=$(CFLAGS)'             \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'LIBS=$(LIBS)'                 \
                            'MACHLIBOBJ=$$($(MACHLIBOBJ))' \
                            'NCARGOBJ=$$($(NCARGOBJ))'     \
                            'MPIOBJ=$$($(MPIOBJ))'         \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'PAKOBJ=$$($(PAKOBJ))'         \
                            'SVIOBJ=$$($(SVIOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            'V5DOBJ=$$($(V5DOBJ))'         \
                            'ARFLAG=$(ARFLAG)'             \
                            $(LIBRADTN)

#-----------------------------------------------------------------------
#
# Make iplib.a for the verification subsystem
#
#-----------------------------------------------------------------------

$(IPLIB):
	cd $(IPLIBDIR); make $(MAKEOPT)                    \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FIXFLAGS=$(FIXFLAGS)'         \
                            'CFLAGS=$(CFLAGS)'             \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'LIBS=$(LIBS)'                 \
                            'MACHLIBOBJ=$$($(MACHLIBOBJ))' \
                            'NCARGOBJ=$$($(NCARGOBJ))'     \
                            'MPIOBJ=$$($(MPIOBJ))'         \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'PAKOBJ=$$($(PAKOBJ))'         \
                            'SVIOBJ=$$($(SVIOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            'V5DOBJ=$$($(V5DOBJ))'         \
                            $(IPLIB)

#-----------------------------------------------------------------------
#
# Compile the ARPS Data Analysis System.
#
# The executable is  adas
#
#-----------------------------------------------------------------------

$(ADASEXE): $(LIBARPS) $(LIBRADTN)
	cd $(ADASDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'ARPSDIR=$(ARPSDIR)'           \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            $(ADASEXE)

$(ADASMPEXE): $(LIBARPS) $(LIBRADTN)
	cd $(ADASDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'ARPSDIR=$(ARPSDIR)'           \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            $(ADASMPEXE)

$(LIBADAS):
	cd $(ADASDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'ARFLAG=$(ARFLAG)'             \
                            $(LIBADAS)


fsl2snd:
	cd $(ADASDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'ARPSDIR=$(ARPSDIR)'           \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            fsl2snd

#-----------------------------------------------------------------------
#
# Compile the ARPS 3D variational Data Analysis System.
#
# The executable is  arps3dvar
#
#-----------------------------------------------------------------------

$(3DVAREXE): $(LIBARPS) $(LIBADAS) $(LIBRADTN)
	cd $(3DVARDIR); make $(MAKEOPT)                    \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'MPIOBJ=$$($(MPIOBJ))'         \
                            $(3DVAREXE)

$(3DVARMPEXE): $(LIBARPS) $(LIBADAS) $(LIBRADTN)
	cd $(3DVARDIR); make $(MAKEOPT)                    \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'MPIOBJ=$$($(MPIOBJ))'         \
                            $(3DVARMPEXE)

meso2lso:
	cd $(MLSODIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'CC=$(CC)'                     \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'FIXFLAGS=$(FIXFLAGS)'         \
                            'CFLAGS=$(CFLAGS)'             \
                            'LDFLAGS=$(LDFLAGS)'           \
                            meso2lso

#-----------------------------------------------------------------------
#
# Compile the Ensemble Kalman Filter Data Assimilation System.
#
# The executable is arpsenkf (mtong)
#
#-----------------------------------------------------------------------

$(ENKFEXE): $(LIBARPS) $(LIBADAS) $(LIBRADTN) $(PRDLIB)
	cd $(ENKFDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            $(ENKFEXE)

$(ENKFMPEXE): $(LIBARPS) $(LIBADAS) $(LIBRADTN)
	cd $(ENKFDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            $(ENKFMPEXE)


$(LIBENKF):
	cd $(ENKFDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'ARFLAG=$(ARFLAG)'             \
                            $(LIBENKF)

$(PRDLIB):
	cd $(RDREMULDIR); make $(MAKEOPT)                  \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'ARFLAG=$(ARFLAG)'             \
                            $(PRDLIB)

#-----------------------------------------------------------------------
#
# Compile and link ARPS AGR executable, arpsagr
#
#-----------------------------------------------------------------------

$(AGREXE): $(ARPSSOLVER)
	cd $(AGRDIR);  make $(MAKEOPT)                 \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'MACHLIBOBJ=$$($(MACHLIBOBJ))' \
                            $(AGREXE)

#-----------------------------------------------------------------------
#
# Compile and link ARPS ASSIM executable, arpsassim
#
#-----------------------------------------------------------------------

$(ASSIMEXE): $(ARPSSOLVER)
	cd $(ASSIMDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'MACHLIBOBJ=$$($(MACHLIBOBJ))' \
                            $(ASSIMEXE)

#-----------------------------------------------------------------------
#
# Graphic plotting:
#
#-----------------------------------------------------------------------

$(PLTMAXEXE): $(LIBARPS) $(LIBADAS) $(LIBZXPOST) $(LIBZXNCAR)
	cd $(PLTDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'ZXPLOTOBJ=$$($(ZXPLOTOBJ))'   \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           $(PLTMAXEXE)

$(PLTNCAREXE): $(LIBARPS) $(LIBADAS)  $(LIBZXNCAR)
	cd $(PLTDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'ADASDIR=$(ADASDIR)'           \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           $(PLTNCAREXE)

$(PLTPOSTEXE): $(LIBARPS) $(LIBADAS)  $(LIBZXPOST)
	cd $(PLTDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'ADASDIR=$(ADASDIR)'           \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           $(PLTPOSTEXE)

$(PLTPOSTMPEXE): $(LIBARPS) $(LIBADAS) $(PRDLIB) $(LIBZXPOST)
	cd $(PLTDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'ADASDIR=$(ADASDIR)'           \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           $(PLTPOSTMPEXE)

$(PLTNCARMPEXE): $(LIBARPS) $(LIBADAS) $(PRDLIB) $(LIBZXNCAR)
	cd $(PLTDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'ADASDIR=$(ADASDIR)'           \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           $(PLTNCARMPEXE)

$(PLTGRIDEXE): $(LIBZXPOST) $(LIBZXNCAR) $(LIBADAS)
	cd $(PLTDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'ZXPLOTOBJ=$$($(ZXPLOTOBJ))'   \
                           $(PLTGRIDEXE)

#-----------------------------------------------------------------------
#
# Compile and link data conversion program, arpscvt
#
#-----------------------------------------------------------------------

$(ARPSCVTEXE): $(LIBARPS)
	cd $(CVTDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(ARPSCVTEXE)

#-----------------------------------------------------------------------
#
# Compile and link data conversion program, arpscvtobs
#
#-----------------------------------------------------------------------

$(ARPSCVTOBSEXE):
	cd $(CVTOBSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'INCFILE=$(INCFILE)'           \
                           'GEMINC=$(GEMINC)'             \
                           $(ARPSCVTOBSEXE)

#-----------------------------------------------------------------------
#
# Compile and link microphysics budget program, arpsmpbudget
#
#-----------------------------------------------------------------------

$(ARPSMPBUDGETEXE): $(LIBARPS)
	cd $(ARPSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(ARPSMPBUDGETEXE)

#-----------------------------------------------------------------------
#
# Compile and link trajectory program
#
#-----------------------------------------------------------------------

$(TRAJCEXE): $(LIBARPS)
	cd $(TRAJCDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             $(TRAJCEXE)

$(CALCTEXE): $(LIBARPS)
	cd $(TRAJCDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             $(CALCTEXE)

$(PLT3DEXE): $(LIBZXNCAR)
	cd $(TRAJCDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             plt3dtrajc

$(PLT1DEXE): $(LIBZXNCAR)
	cd $(TRAJCDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             plt1dtrajc

#-----------------------------------------------------------------------
#
# Compile and link gridded data interpolation program - arpsintrp.
#
#-----------------------------------------------------------------------

$(INTRPEXE): $(LIBARPS) $(LIBADAS)
	cd $(INTRPDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             $(INTRPEXE)


$(INTRP_LSEXE): $(LIBARPS)
	cd $(INTRPDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             $(INTRP_LSEXE)

$(RADARINTRPEXE): $(LIBARPS) $(LIBADAS) $(LIBENKF)
	cd $(INTRPDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                              $(RADARINTRPEXE)

$(SUBDMNEXE): $(LIBARPS) $(LIBADAS)
	cd $(INTRPDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             $(SUBDMNEXE)

$(HDFSUBDMNEXE): $(LIBARPS)
	cd $(INTRPDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             $(HDFSUBDMNEXE)

$(HDFSUBDMNMPEXE): $(LIBARPS)
	cd $(INTRPDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             $(HDFSUBDMNMPEXE)

#-----------------------------------------------------------------------
#
# Compile and link gridded data interpolation program - arpstintrp.
#
#-----------------------------------------------------------------------

$(TINTRPEXE): $(LIBARPS)
	cd $(TINTRPDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             $(TINTRPEXE)

$(TINTRPMPEXE): $(LIBARPS)
	cd $(TINTRPDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'LIBS=$(LIBS)'                 \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             $(TINTRPMPEXE)

#-----------------------------------------------------------------------
#
# Compile and link soil conversion program, arpssoil
#
#-----------------------------------------------------------------------

$(ARPSSOILEXE): $(LIBARPS) $(LIBADAS)
	cd $(SOILDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            $(ARPSSOILEXE)

#-----------------------------------------------------------------------
#
# Compile and link data conversion program for systems
# and generate executable arpsprt
#
#-----------------------------------------------------------------------

$(PRTEXE): $(LIBARPS)
	cd $(PRTDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(PRTEXE)

#-----------------------------------------------------------------------
#
# Compile and link ARPS to GEMPAK data conversion program
# and generate executable arps2gem
#
#-----------------------------------------------------------------------
#

$(A2GEXE): $(LIBARPS) $(LIBADAS)
	cd $(E2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'INCFILE=$(INCFILE)'           \
                           'GEMINC=$(GEMINC)'             \
                           $(A2GEXE)

#-----------------------------------------------------------------------
#
# Compile and link ARPS to netCDF data conversion program
# and generate executable arps2ncdf
#
#-----------------------------------------------------------------------
#

$(A2NEXE): $(LIBARPS) $(LIBADAS)
	cd $(E2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(A2NEXE)

$(A2WIIEXE): $(LIBARPS)
	cd $(E2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(A2WIIEXE)

$(A2EEXE): $(LIBARPS) $(LIBADAS)
	cd $(E2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(A2EEXE)

$(A2WEXE): $(LIBARPS) $(LIBADAS) $(LIBWRFINT)
	cd $(A2WDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'FTN=$(FTN)'                   \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(A2WEXE)

$(A2WMPEXE): $(LIBARPS) $(LIBADAS) $(LIBWRFINT) $(LIBPHDF5_DEP)
	cd $(A2WDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'PHDF5OBJ=$$($(PHDF5OBJ))'     \
                           'FTN=$(FTN)'                   \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(A2WMPEXE)

$(A4WEXE): $(LIBARPS) $(LIBWRFINT) $(LIBWRFNET)
	cd $(A4WDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'FTN=$(FTN)'                   \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(A4WEXE)

$(A4WMPEXE): $(LIBARPS) $(LIBWRFINT) $(LIBWRFNET)
	cd $(A4WDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'PHDF5OBJ=$$($(PHDF5OBJ))'     \
                           'FTN=$(FTN)'                   \
                           'FTNMP=$(FTNMP)'               \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(A4WMPEXE)

$(N2AEXE): $(LIBARPS) $(LIBADAS) $(LIBWRFINT) $(LIBWRFNET)
	cd $(N2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'FTN=$(FTN)'                   \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           'GEMOBJ=$$($(GEMOBJ))'         \
                           'CRTMOBJ=$$($(CRTMOBJ))'       \
                           'CITMOBJ=$$($(CITMOBJ))'       \
                           $(N2AEXE)

$(N2AMPEXE): $(LIBARPS) $(LIBADAS) $(LIBWRFINT) $(LIBWRFNET)
	cd $(N2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'PHDF5OBJ=$$($(PHDF5OBJ))'     \
                           'FTN=$(FTN)'                   \
                           'FTNMP=$(FTNMP)'               \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           'GEMOBJ=$$($(GEMOBJ))'         \
                           'CRTMOBJ=$$($(CRTMOBJ))'       \
                           'CITMOBJ=$$($(CITMOBJ))'       \
                           $(N2AMPEXE)

$(WSEXE): $(LIBARPS) $(LIBWRFINT)
	cd $(A2WDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'FTN=$(FTN)'                   \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(WSEXE)

$(W2AEXE): $(LIBARPS) $(LIBADAS) $(LIBWRFINT)
	cd $(W2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'FTN=$(FTN)'                   \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           'GEMOBJ=$$($(GEMOBJ))'         \
                           'CRTMOBJ=$$($(CRTMOBJ))'       \
                           'CITMOBJ=$$($(CITMOBJ))'       \
                           $(W2AEXE)

$(W2AMPEXE): $(LIBARPS) $(LIBADAS) $(LIBWRFINT) $(LIBPHDF5_DEP)
	cd $(W2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'PHDF5OBJ=$$($(PHDF5OBJ))'     \
                           'FTN=$(FTN)'                   \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           'GEMOBJ=$$($(GEMOBJ))'         \
                           'CRTMOBJ=$$($(CRTMOBJ))'       \
                           'CITMOBJ=$$($(CITMOBJ))'       \
                           $(W2AMPEXE)

$(WEXTSNDEXE): $(LIBARPS) $(LIBADAS) $(LIBWRFINT)
	cd $(W2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'FTN=$(FTN)'                   \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(WEXTSNDEXE)

$(JOINWRFEXE): $(LIBARPS)
	cd $(W2ADIR); make $(MAKEOPT)                           \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(JOINWRFEXE)

$(A2CEXE): $(LIBARPS)
	cd $(A2CDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'FTN=$(FTN)'                   \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(A2CEXE)

$(A2CMPEXE): $(LIBARPS)
	cd $(A2CDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'PHDF5OBJ=$$($(PHDF5OBJ))'     \
                           'FTN=$(FTN)'                   \
                           'FTNMP=$(FTNMP)'               \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(A2CMPEXE)

$(C2AEXE): $(LIBARPS) $(LIBADAS)
	cd $(C2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'FTN=$(FTN)'                   \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           'GEMOBJ=$$($(GEMOBJ))'         \
                           'CRTMOBJ=$$($(CRTMOBJ))'       \
                           'CITMOBJ=$$($(CITMOBJ))'       \
                           $(C2AEXE)

$(C2AMPEXE): $(LIBARPS) $(LIBADAS)
	cd $(C2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'PHDF5OBJ=$$($(PHDF5OBJ))'     \
                           'FTN=$(FTN)'                   \
                           'FTNMP=$(FTNMP)'               \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           'GEMOBJ=$$($(GEMOBJ))'         \
                           'CRTMOBJ=$$($(CRTMOBJ))'       \
                           'CITMOBJ=$$($(CITMOBJ))'       \
                           $(C2AMPEXE)

#----------------------------------------------------------
#
# External libraries
#
#----------------------------------------------------------

$(LIBPHDF5):
	cd $(PHDF5DIR); make $(MAKEOPT)                   \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS) -w'      \
                           'PHDF5PATH=$(PHDF5PATH)'       \
                           'FC=$(FTN)'                    \
                           'FFLAGS=$(FFLAGS) -I$(PHDF5PATH)/lib -w' \
                           $(LIBPHDF5)

$(LIBWRFINT):
	cd $(WRFINTDIR); make $(MAKEOPT)                  \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'FC=$(FTN)'                    \
                           'CC=$(CC)'                     \
                           'FFLAGS=$(FFLAGS)'             \
                           'FREEFLAGS=$(FREEFLAGS) -w'    \
                           'CFLAGS=$(CFLAGS) -w'          \
                           $(LIBWRFINT)

$(LIBWRFNET):
	cd $(WRFNETDIR); make $(MAKEOPT)                  \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'FC=$(FTN)'                    \
                           'CC=$(CC)'                     \
                           'FFLAGS=$(FFLAGS)'             \
                           'FREEFLAGS=$(FREEFLAGS) -w'    \
                           'CFLAGS=$(CFLAGS) -w'          \
                           $(LIBWRFNET)

$(LIBG2):
	cd $(G2DIR); make $(MAKEOPT)                      \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           'FC=$(FTN)'                    \
                           'CC=$(CC)'                     \
                           'FFLAGS=$(FFLAGS)'             \
                           'FIXFLAGS=$(FIXFLAGS)'         \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'CFLAGS=$(CFLAGS)'             \
                           $(LIBG2)

#-----------------------------------------------------------------------
#
# Generate verification statistics comparing analysis or forecast
# to observations.
#
#-----------------------------------------------------------------------

$(DIFOBSEXE): $(LIBARPS) $(LIBADAS)
	cd $(ADASDIR); make $(MAKEOPT)                    \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(DIFOBSEXE)

#-----------------------------------------------------------------------
#
# Compare two data dumps and write out the difference in the
# history dump format.
#
#-----------------------------------------------------------------------

$(DIFFEXE): $(LIBARPS) $(LIBADAS)
	cd $(DIFDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(DIFFEXE)

#-----------------------------------------------------------------------
#
# Write out differenec between accumulated rainfall at two times.
#
#-----------------------------------------------------------------------

$(RAINDIFFEXE): $(LIBARPS)
	cd $(DIFDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(RAINDIFFEXE)

#-----------------------------------------------------------------------
#
# ARPS ensemble programs
#
#-----------------------------------------------------------------------

$(ENSCVEXE): $(LIBARPS) $(LIBADAS) $(LIBZXPOST) $(LIBZXNCAR)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'ZXPLOTOBJ=$$($(ZXPLOTOBJ))'   \
                           $(ENSCVEXE)

$(ENSBCEXE): $(LIBARPS)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(ENSBCEXE)

$(ENSICEXE): $(LIBARPS)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(ENSICEXE)
$(ENSICMPEXE): $(LIBARPS)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           $(ENSICMPEXE)

$(EPOSTEXE): $(LIBARPS) $(LIBADAS)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'INCFILE=$(INCFILE)'           \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'GEMOBJ=$$($(GEMOBJ))'         \
                           'CRTMOBJ=$$($(CRTMOBJ))'       \
                           'CITMOBJ=$$($(CITMOBJ))'       \
                           $(EPOSTEXE)

$(EPOSTMPEXE): $(LIBARPS) $(LIBADAS)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'INCFILE=$(INCFILE)'           \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'GEMOBJ=$$($(GEMOBJ))'         \
                           'CRTMOBJ=$$($(CRTMOBJ))'       \
                           'CITMOBJ=$$($(CITMOBJ))'       \
                           $(EPOSTMPEXE)

$(B2GEXE): $(LIBARPS) $(LIBADAS)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'INCFILE=$(INCFILE)'           \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(B2GEXE)

$(ENSANAEXE): $(LIBARPS) $(LIBADAS)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'INCFILE=$(INCFILE)'           \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'GEMOBJ=$$($(GEMOBJ))'         \
                           $(ENSANAEXE)

$(ENSCALEXE): $(LIBARPS) $(LIBADAS)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'INCFILE=$(INCFILE)'           \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'GEMOBJ=$$($(GEMOBJ))'         \
                           $(ENSCALEXE)

$(ENKFICEXE): $(LIBARPS) $(LIBENKF) $(LIBADAS)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(ENKFICEXE)


$(RNDPRTEXE): $(LIBENKF) $(LIBARPS)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
			   $(RNDPRTEXE)

$(ENSSCORESEXE): $(LIBARPS)
	cd $(ENSDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(ENSSCORESEXE)
#-----------------------------------------------------------------------
#
# Sample program to read ARPS history dump.
#
#-----------------------------------------------------------------------

$(ARPSREADEXE): $(LIBARPS)
	cd $(ARPSDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            $(ARPSREADEXE)

#-----------------------------------------------------------------------
#
# Compile and link the program that extracts a sounding profile at
# a given column in the model history data.
#
# The executable is arpsextsnd
#
#-----------------------------------------------------------------------
#

$(EXTSNDEXE): $(LIBARPS) $(LIBADAS)
	cd $(EXTSNDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             $(EXTSNDEXE)

$(EXTSNDMPEXE): $(LIBARPS) $(LIBADAS)
	cd $(EXTSNDIR); make $(MAKEOPT)                     \
                             'TOPDIR=$(TOPDIR)'             \
                             'BINDIR=$(BINDIR)'             \
                             'INCLDIR=$(INCLDIR)'           \
                             'LIBDIR=$(LIBDIR)'             \
                             'FTN=$(FTN)'                   \
                             'ARPS_LD=$(ARPS_LD)'           \
                             'FFLAGS=$(FFLAGS)'             \
                             'FFLAGS_main=$(FFLAGS_main)'   \
                             'FREEFLAGS=$(FREEFLAGS)'       \
                             'LDFLAGS=$(LDFLAGS)'           \
                             $(EXTSNDMPEXE)

$(SKEWTNEXE): $(LIBZXNCAR)
	cd $(SKEWTDIR); make $(MAKEOPT)                   \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'FIXFLAGS=$(FIXFLAGS)'         \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           $(SKEWTNEXE)

$(SKEWTPEXE): $(LIBZXPOST)
	cd $(SKEWTDIR); make $(MAKEOPT)                   \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'FIXFLAGS=$(FIXFLAGS)'         \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'CPP=$(CPP)'                   \
                           'CPPFLAGS=$(CPPFLAGS)'         \
                           $(SKEWTPEXE)

#-----------------------------------------------------------------------
#
# Compile surface data pre-processing program arpssfcxx
#
# The executable is  arpssfc
#
#-----------------------------------------------------------------------

$(ARPSSFCEXE): $(LIBARPS) $(LIBADAS)
	cd $(SFCDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(ARPSSFCEXE)

#-----------------------------------------------------------------------
#
# Compile and link radar remapping program, 88d2arps.
#
# The executable is 88d2arps.
#
#-----------------------------------------------------------------------

$(88D2AEXE): $(LIBARPS) $(LIBADAS) $(LIBA2IO)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(88D2AEXE)

$(88D2ASEXE): $(LIBARPS) $(LIBADAS)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(88D2ASEXE)

$(88D2AMPEXE): $(LIBARPS) $(LIBADAS) $(LIBA2IO)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(88D2AMPEXE)

$(88D2AEXE_FAKE): $(LIBARPS) $(LIBADAS)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(88D2AEXE_FAKE)

$(LIBA2IO):
	cd $(88DDIR)/a2io; make $(MAKEOPT)                \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           $(LIBA2IO)

#-----------------------------------------------------------------------
#
# Compile and link radar column plotting program, nids2arps
#
# The executable is nids2arps
#
#-----------------------------------------------------------------------

$(NIDS2AEXE): $(LIBARPS) $(LIBADAS)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(NIDS2AEXE)

$(NIDS2AMPEXE): $(LIBARPS) $(LIBADAS)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(NIDS2AMPEXE)

#-----------------------------------------------------------------------
#
# Compile and link NetCDF radar data remapping program, ncrad2arps
#
# The executable is ncrad2arps
#
#-----------------------------------------------------------------------

$(NCRAD2AEXE): $(LIBARPS) $(LIBADAS)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(NCRAD2AEXE)

$(NCRAD2AMPEXE): $(LIBARPS) $(LIBADAS)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(NCRAD2AMPEXE)

$(CASA2AEXE): $(LIBARPS) $(LIBADAS)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(CASA2AEXE)

#-----------------------------------------------------------------------
#
# Compile and link the simulated radar data generating program.
#
#-----------------------------------------------------------------------

$(ARPS2RADEXE): $(LIBARPS) $(LIBADAS)
	cd $(88DDIR); make $(MAKEOPT)                    \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(ARPS2RADEXE)

#-----------------------------------------------------------------------
#
# Compile and link the simulated radar data generating program.
#
#-----------------------------------------------------------------------

$(RADEMULEXE): $(LIBARPS) $(LIBADAS)
	cd $(88DDIR); make $(MAKEOPT)                    \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'NETOBJ=$$($(NETOBJ))'         \
                           $(RADEMULEXE)

#-----------------------------------------------------------------------
#
# Compile and link the scattering amplitudes generating program.
#
#-----------------------------------------------------------------------

$(TMATRIXEXE):
	cd $(RDREMULDIR); make $(MAKEOPT)                    \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(TMATRIXEXE)

#-----------------------------------------------------------------------
#
# Compile and link a program to convert binary emulated
# radar data to NetCDF
#
#-----------------------------------------------------------------------

$(RADB2CDFEXE): $(LIBARPS)
	cd $(88DDIR); make $(MAKEOPT)                    \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'NETOBJ=$$($(NETOBJ))'         \
                           $(RADB2CDFEXE)

#-----------------------------------------------------------------------
#
# Compile and link a program to sectorize NetCDF radar files
#
#-----------------------------------------------------------------------


$(RADSECTEXE): $(LIBARPS)
	cd $(88DDIR); make $(MAKEOPT)                    \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'NETOBJ=$$($(NETOBJ))'         \
                           $(RADSECTEXE)

#-----------------------------------------------------------------------
#
# Compile and link a program to compute attenuation along path
# to a specific target location
#
#-----------------------------------------------------------------------

$(ATTENEXE): $(LIBARPS) $(LIBADAS)
	cd $(88DDIR); make $(MAKEOPT)                    \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'NETOBJ=$$($(NETOBJ))'         \
                           $(ATTENEXE)

#-----------------------------------------------------------------------
#
# Compile and link the radar mosaicking program.
#
#-----------------------------------------------------------------------

$(RADMOSEXE): $(LIBARPS) $(LIBADAS)
	cd $(ADASDIR); make $(MAKEOPT)                    \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(RADMOSEXE)

#-----------------------------------------------------------------------
#
# Compile and link radar column plotting program, pltradcol
#
# The executable is pltradcol
#
#-----------------------------------------------------------------------

$(PLTRADCOLEXE): $(LIBARPS) $(LIBZXPOST) $(LIBZXNCAR)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'ZXPLOTOBJ=$$($(ZXPLOTOBJ))'   \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(PLTRADCOLEXE)

ppiplt: $(LIBARPS) $(LIBZXPOST) $(LIBZXNCAR)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'ZXPLOTOBJ=$$($(ZXPLOTOBJ))'   \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           ppiplt

$(RADARPLTNCAREXE): $(LIBARPS) $(LIBZXNCAR)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(RADARPLTNCAREXE)

$(RADARPLTPOSTEXE): $(LIBARPS) $(LIBZXPOST)
	cd $(88DDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(RADARPLTPOSTEXE)

#-----------------------------------------------------------------------
#
# Compile and link satellite data conversion program, mci2arps
#
# The executable is mci2arps
#
#-----------------------------------------------------------------------

$(MOSAICEXE): $(LIBARPS)
	cd $(MOSAICDIR); make $(MAKEOPT)                  \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(MOSAICEXE)

$(MCI2AEXE): $(LIBARPS)
	cd $(MCIDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(MCI2AEXE)

$(SAT2AEXE): $(LIBARPS)
	cd $(MCIDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(SAT2AEXE)
$(MERGESATEXE): $(LIBARPS)
	cd $(MCIDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(MERGESATEXE)
$(SATHDF5AEXE): $(LIBARPS)
	cd $(MCIDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(SATHDF5AEXE)

#-----------------------------------------------------------------------
#
# Compile and link radar column plotting program, pltsatfld
#
# The executable is pltsatfld
#
#-----------------------------------------------------------------------

$(PLTSATFLDEXE): $(LIBARPS)
	cd $(MCIDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSDIR=$(ARPSDIR)'           \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(PLTSATFLDEXE)

#-----------------------------------------------------------------------
#
# Compile and link terrain data preprocessing program arpstern.
#
# The executable is arpstern.
#
#-----------------------------------------------------------------------

$(TERNEXE): $(LIBARPS)
	cd $(TERNDIR); make $(MAKEOPT)                    \
                            'TOPDIR=$(TOPDIR)'            \
                            'BINDIR=$(BINDIR)'            \
                            'INCLDIR=$(INCLDIR)'          \
                            'LIBDIR=$(LIBDIR)'            \
                            'FTN=$(FTN)'                  \
                            'ARPS_LD=$(ARPS_LD)'          \
                            'FFLAGS=$(FFLAGS)'            \
                            'FFLAGS_main=$(FFLAGS_main)'  \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'          \
                            'LIBS=$(LIBS)'                \
                            $(TERNEXE)

#-----------------------------------------------------------------------
#
# Compile and link terrain data preprocessing program arpstrn.
#
# The executable is arpstrn.
#
#-----------------------------------------------------------------------

$(TRNEXE): $(LIBARPS) $(LIBZXPOST) $(LIBZXNCAR)
	cd $(TRNDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'ZXPLOTOBJ=$$($(ZXPLOTOBJ))'   \
                           $(TRNEXE)

$(DIR1DEGEXE):
	cd $(TERNDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(DIR1DEGEXE)

$(DIR5MINEXE):
	cd $(TERNDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(DIR5MINEXE)

$(DIR30SECEXE):
	cd $(TERNDIR); make $(MAKEOPT)                    \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(DIR30SECEXE)

#-----------------------------------------------------------------------
#
# Compile and link terrain data merger and edge smoother.
#
# The executable is mergetrn.
#
#-----------------------------------------------------------------------

$(MERGETRNEXE): $(LIBARPS) $(LIBADAS) $(LIBZXPOST) $(LIBZXNCAR)
	cd $(TRNDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'LIBDIR=$(LIBDIR)'             \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'ZXPLOTOBJ=$$($(ZXPLOTOBJ))'   \
                           $(MERGETRNEXE)

#-----------------------------------------------------------------------
#
# Program ext2arps to convert external data set to the ARPS grid.
#
#-----------------------------------------------------------------------

$(E2AEXE): $(LIBARPS) $(LIBADAS) $(LIBG2_DEP)
	cd $(E2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'GRIB2OBJ=$$($(GRIB2OBJ))'     \
                           'NETOBJ=$$($(NETOBJ))'         \
                           $(E2AEXE)

$(E2AMPEXE): $(LIBARPS) $(LIBADAS) $(LIBG2_DEP)
	cd $(E2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           'GRIB2OBJ=$$($(GRIB2OBJ))'     \
                           'NETOBJ=$$($(NETOBJ))'         \
                           $(E2AMPEXE)

$(L2AEXE): $(LIBARPS) $(LIBADAS)
	cd $(E2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(L2AEXE)

$(G2AEXE): $(LIBARPS) $(LIBADAS)
	cd $(E2ADIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'INCFILE=$(INCFILE)'           \
                           'GEMINC=$(GEMINC)'             \
                           $(G2AEXE)

$(AVNEXE): $(LIBARPS) $(LIBADAS)
	cd $(AVNDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'LIBS=$(LIBS)'                 \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS=$(FFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           $(AVNEXE)

#-----------------------------------------------------------------------
#
# Compile the ARPS verification package.
#
# The executable name is 'arpsverif'.
#
#-----------------------------------------------------------------------

$(VERIFEXE): $(LIBARPS) $(LIBADAS) $(IPLIB)
	cd $(VERIFDIR); make $(MAKEOPT)                    \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'ARPSDIR=$(ARPSDIR)'           \
                            'MPIOBJ=$$($(MPIOBJ))'         \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            $(VERIFEXE)

$(VERIFMPEXE): $(LIBARPS) $(LIBADAS) $(IPLIB)
	cd $(VERIFDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'MPIOBJ=$$($(MPIOBJ))'         \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(VERIFMPEXE)

$(TYPHOONEXE): $(LIBARPS)
	cd $(VERIFDIR); make $(MAKEOPT)                    \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'ARPSDIR=$(ARPSDIR)'           \
                            'MPIOBJ=$$($(MPIOBJ))'         \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'FREEFLAGS=$(FREEFLAGS)'       \
                            'LDFLAGS=$(LDFLAGS)'           \
                            $(TYPHOONEXE)

#-----------------------------------------------------------------------
#
# Compile and link program to write retrieval output
#
#-----------------------------------------------------------------------

$(WTRETCOLEXE): $(LIBARPS)
	cd $(WTRETCOLDIR); make $(MAKEOPT)                     \
                                'TOPDIR=$(TOPDIR)'             \
                                'BINDIR=$(BINDIR)'             \
                                'INCLDIR=$(INCLDIR)'           \
                                'LIBDIR=$(LIBDIR)'             \
                                'FTN=$(FTN)'                   \
                                'ARPS_LD=$(ARPS_LD)'           \
                                'FFLAGS=$(FFLAGS)'             \
                                'FFLAGS_main=$(FFLAGS_main)'   \
                                'FREEFLAGS=$(FREEFLAGS)'       \
                                'LDFLAGS=$(LDFLAGS)'           \
                                $(WTRETCOLEXE)

#-----------------------------------------------------------------------
#
# Create data file splitter and joiner for general arps_mp
#
#-----------------------------------------------------------------------

$(SPLITEXE): $(LIBARPS)
	cd $(MPDIR); make $(MAKEOPT)                     \
                          'TOPDIR=$(TOPDIR)'             \
                          'BINDIR=$(BINDIR)'             \
                          'INCLDIR=$(INCLDIR)'           \
                          'LIBDIR=$(LIBDIR)'             \
                          'FTN=$(FTN)'                   \
                          'ARPS_LD=$(ARPS_LD)'           \
                          'FFLAGS=$(FFLAGS)'             \
                          'FFLAGS_main=$(FFLAGS_main)'   \
                          'FREEFLAGS=$(FREEFLAGS)'       \
                          'LDFLAGS=$(LDFLAGS)'           \
                          'HDFOBJ=$$($(HDFOBJ))'         \
                          $(SPLITEXE)

$(SPLITMPEXE): $(LIBARPS)
	cd $(MPDIR); make $(MAKEOPT)                     \
                          'TOPDIR=$(TOPDIR)'             \
                          'BINDIR=$(BINDIR)'             \
                          'INCLDIR=$(INCLDIR)'           \
                          'LIBDIR=$(LIBDIR)'             \
                          'FTN=$(FTN)'                   \
                          'ARPS_LD=$(ARPS_LD)'           \
                          'FFLAGS=$(FFLAGS)'             \
                          'FFLAGS_main=$(FFLAGS_main)'   \
                          'FREEFLAGS=$(FREEFLAGS)'       \
                          'LDFLAGS=$(LDFLAGS)'           \
                          'HDFOBJ=$$($(HDFOBJ))'         \
                          $(SPLITMPEXE)

$(JOINSEXE): $(LIBARPS)
	cd $(MPDIR); make $(MAKEOPT)                     \
                          'TOPDIR=$(TOPDIR)'             \
                          'BINDIR=$(BINDIR)'             \
                          'INCLDIR=$(INCLDIR)'           \
                          'LIBDIR=$(LIBDIR)'             \
                          'FTN=$(FTN)'                   \
                          'ARPS_LD=$(ARPS_LD)'           \
                          'FFLAGS=$(FFLAGS)'             \
                          'FFLAGS_main=$(FFLAGS_main)'   \
                          'FREEFLAGS=$(FREEFLAGS)'       \
                          'LDFLAGS=$(LDFLAGS)'           \
                          'HDFOBJ=$$($(HDFOBJ))'         \
                          $(JOINSEXE)

$(JOINEXE): $(LIBARPS)
	cd $(MPDIR); make $(MAKEOPT)                     \
                          'TOPDIR=$(TOPDIR)'             \
                          'BINDIR=$(BINDIR)'             \
                          'INCLDIR=$(INCLDIR)'           \
                          'LIBDIR=$(LIBDIR)'             \
                          'FTN=$(FTN)'                   \
                          'ARPS_LD=$(ARPS_LD)'           \
                          'FFLAGS=$(FFLAGS)'             \
                          'FFLAGS_main=$(FFLAGS_main)'   \
                          'FREEFLAGS=$(FREEFLAGS)'       \
                          'LDFLAGS=$(LDFLAGS)'           \
                          'HDFOBJ=$$($(HDFOBJ))'         \
                          $(JOINEXE)
joinbin2hdf: $(LIBARPS)
	cd $(MPDIR); make $(MAKEOPT)                     \
                          'TOPDIR=$(TOPDIR)'             \
                          'BINDIR=$(BINDIR)'             \
                          'INCLDIR=$(INCLDIR)'           \
                          'LIBDIR=$(LIBDIR)'             \
                          'FTN=$(FTN)'                   \
                          'ARPS_LD=$(ARPS_LD)'           \
                          'FFLAGS=$(FFLAGS)'             \
                          'FFLAGS_main=$(FFLAGS_main)'   \
                          'FREEFLAGS=$(FREEFLAGS)'       \
                          'LDFLAGS=$(LDFLAGS)'           \
                          joinbin2hdf

joinhdf: $(LIBARPS)
	cd $(MPDIR); make $(MAKEOPT)                     \
                          'TOPDIR=$(TOPDIR)'             \
                          'BINDIR=$(BINDIR)'             \
                          'INCLDIR=$(INCLDIR)'           \
                          'LIBDIR=$(LIBDIR)'             \
                          'FTN=$(FTN)'                   \
                          'ARPS_LD=$(ARPS_LD)'           \
                          'FFLAGS=$(FFLAGS)'             \
                          'FFLAGS_main=$(FFLAGS_main)'   \
                          'FREEFLAGS=$(FREEFLAGS)'       \
                          'LDFLAGS=$(LDFLAGS)'           \
                          joinhdf
splithdf: $(LIBARPS)
	cd $(MPDIR); make $(MAKEOPT)                     \
                          'TOPDIR=$(TOPDIR)'             \
                          'BINDIR=$(BINDIR)'             \
                          'INCLDIR=$(INCLDIR)'           \
                          'LIBDIR=$(LIBDIR)'             \
                          'FTN=$(FTN)'                   \
                          'ARPS_LD=$(ARPS_LD)'           \
                          'FFLAGS=$(FFLAGS)'             \
                          'FFLAGS_main=$(FFLAGS_main)'   \
                          'FREEFLAGS=$(FREEFLAGS)'       \
                          'LDFLAGS=$(LDFLAGS)'           \
                          splithdf

splithdf_mpi: $(LIBARPS)
	cd $(MPDIR); make $(MAKEOPT)                     \
                          'TOPDIR=$(TOPDIR)'             \
                          'BINDIR=$(BINDIR)'             \
                          'INCLDIR=$(INCLDIR)'           \
                          'LIBDIR=$(LIBDIR)'             \
                          'FTN=$(FTN)'                   \
                          'ARPS_LD=$(ARPS_LD)'           \
                          'FFLAGS=$(FFLAGS)'             \
                          'FFLAGS_main=$(FFLAGS_main)'   \
                          'FREEFLAGS=$(FREEFLAGS)'       \
                          'LDFLAGS=$(LDFLAGS)'           \
                          splithdf_mpi

$(SPLITANYEXE):
	cd $(MPDIR); make $(MAKEOPT)                     \
                          'TOPDIR=$(TOPDIR)'             \
                          'BINDIR=$(BINDIR)'             \
                          'INCLDIR=$(INCLDIR)'           \
                          'LIBDIR=$(LIBDIR)'             \
                          'FTN=$(FTN)'                   \
                          'ARPS_LD=$(ARPS_LD)'           \
                          'FFLAGS=$(FFLAGS)'             \
                          'FFLAGS_main=$(FFLAGS_main)'   \
                          'FREEFLAGS=$(FREEFLAGS)'       \
                          'LDFLAGS=$(LDFLAGS)'           \
                          'HDFOBJ=$$($(HDFOBJ))'         \
                          'NETOBJ=$$($(NETOBJ))'         \
                          $(SPLITANYEXE)

$(H2GEXE):
	cd $(CVTDIR); make $(MAKEOPT)                     \
                           'TOPDIR=$(TOPDIR)'             \
                           'BINDIR=$(BINDIR)'             \
                           'INCLDIR=$(INCLDIR)'           \
                           'ARPSLIBDIR=$(LIBDIR)'         \
                           'FTN=$(FTN)'                   \
                           'ARPS_LD=$(ARPS_LD)'           \
                           'FFLAGS=$(FFLAGS)'             \
                           'CFLAGS=$(CFLAGS)'             \
                           'FFLAGS_main=$(FFLAGS_main)'   \
                           'FREEFLAGS=$(FREEFLAGS)'       \
                           'LDFLAGS=$(LDFLAGS)'           \
                           'LIBS=$(LIBS)'                 \
                           $(H2GEXE)

#-----------------------------------------------------------------------
#
# EnKF associated programs (mtong)
#
#-----------------------------------------------------------------------

$(INITENSMBLEXE): $(LIBARPS) $(LIBENKF) $(LIBADAS)
	cd $(ENKFDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            $(INITENSMBLEXE)

$(OSSEDATAEXE): $(LIBARPS) $(LIBENKF) $(LIBADAS) $(PRDLIB)
	cd $(ENKFDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            $(OSSEDATAEXE)

$(OBSSTDEXE)  : $(LIBARPS) $(LIBENKF) $(LIBADAS) $(PRDLIB)
	cd $(ENKFDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            $(OBSSTDEXE)

$(POSTINNOVEXE): $(LIBARPS) $(LIBENKF) $(LIBADAS) $(PRDLIB)
	cd $(ENKFDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            $(POSTINNOVEXE)

$(ENDGNSEXE): $(LIBARPS) $(LIBENKF) $(LIBADAS) $(PRDLIB)
	cd $(ENKFDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
		            $(ENDGNSEXE)

$(TESTEXE): $(LIBARPS) $(LIBENKF) $(LIBADAS) $(PRDLIB)
	cd $(ENKFDIR); make $(MAKEOPT)                     \
                            'TOPDIR=$(TOPDIR)'             \
                            'BINDIR=$(BINDIR)'             \
                            'INCLDIR=$(INCLDIR)'           \
                            'LIBDIR=$(LIBDIR)'             \
                            'LIBS=$(LIBS)'                 \
                            'FTN=$(FTN)'                   \
                            'ARPS_LD=$(ARPS_LD)'           \
                            'FFLAGS=$(FFLAGS)'             \
                            'FFLAGS_main=$(FFLAGS_main)'   \
                            'LDFLAGS=$(LDFLAGS)'           \
                            'HDFOBJ=$$($(HDFOBJ))'         \
                            'NETOBJ=$$($(NETOBJ))'         \
                            $(TESTEXE)

#-----------------------------------------------------------------------
#
# Archive the source code, data, and input files for each individule
# programs as well as the tar file including everything
#
#-----------------------------------------------------------------------

$(ALLTAR): $(ALLFILESTAR) $(ARPSTAR) $(AGRTAR) $(CVTTAR) $(PLTTAR)  \
           $(PRTTAR) $(DIFFTAR) $(ADASTAR) $(SFCTAR) $(SOILTAR)     \
           $(TRNTAR) $(EXTSNDTAR) $(ARPSMPTAR) $(READTAR) $(A2GTAR) \
           $(MERGETRNTAR) $(ENSTAR) $(A2NTAR) $(ASSIMTAR)           \
           $(ARPSVERTAR) $(CVTOBSTAR) $(A2ETAR) $(3DVARTAR)         \
           $(DIFOBSTAR)  $(ARPS2RADTAR) $(RADMOSTAR)  $(ATTENTAR)   \
           $(RADEMULTAR) $(RADB2CDFTAR) $(RADSECTTAR) $(88DTAR)     \
           $(A2WTAR)     $(W2ATAR)      $(ZXPLOTTAR)  $(MLSOTAR)    \
           $(BUDGTAR)
           #$(ENKFTAR)


#-----------------------------------------------------------------------
#
# Archive all the source code, data, and input files into one tar file
#
#-----------------------------------------------------------------------

$(ALLFILESTAR): $(ROOTALL)  $(ARPSALL) $(AGRALL)   $(CVTALL)    \
                $(INTRPALL) $(PLTALL)  $(PRTALL)    $(DIFFALL)  \
                $(ADASALL)  $(SFCALL)  $(SOILALL)   $(MPALL)    \
                $(TRNALL)   $(E2AALL)  $(EXTSNDALL) $(88DALL)   \
                $(A2GALL)   $(AVNALL)  $(WTRETCOLALL) $(MCIALL) \
                $(ENSALL)   $(A2NALL)   $(TINTRPALL) $(TERNALL) \
                $(ASSIMALL) $(F2F90SRC) $(VERIFALL)  $(A2EALL)  \
                $(SKEWTALL) $(TRAJCALL) $(RADEMULALL) $(MOSAICALL) \
                $(BUDGALL)
	@$(TAR) $(TAROPT) $@ \
                $(ROOTALL)  $(ARPSALL) $(AGRALL)    $(CVTALL)   \
                $(INTRPALL) $(PLTALL)  $(PRTALL)    $(DIFFALL)  \
                $(ADASALL)  $(SFCALL)  $(SOILALL)   $(MPALL)    \
                $(TRNALL)   $(E2AALL)  $(EXTSNDALL) $(88DALL)   \
                $(A2GALL)   $(AVNALL)  $(WTRETCOLALL) $(MCIALL)
	@$(TAR) $(TARADD) $@ \
                $(ENSALL)   $(A2NALL)   $(TINTRPALL) $(TERNALL) \
                $(ASSIMALL) $(F2F90SRC) $(VERIFALL)  $(CVTOBSALL) \
                $(DIFOBSALL) $(ARPS2RADALL) $(RADMOSALL)        \
                $(A2EALL)   $(A2WALL)   $(W2AALL)  $(ZXPLOTALL) \
                $(3DVARALL) $(A4WALL)   $(N2AALL)     \
                $(LIBPHDF5ALL) $(LIBWRFINTALL) $(LIBWRFNETALL)  \
                $(SKEWTALL) $(LIBG2ALL) $(TRAJCALL) $(RADEMULALL) \
                $(RADB2CDFALL) $(RADSECTALL) $(ATTENALL) $(BUDGALL)


$(ALLFILESGZ): $(ALLFILESTAR)
	$(ZIP) -f $(ALLFILESTAR)

#-----------------------------------------------------------------------
#
# Archive all the source code, data, and input files into one tar file
#
#-----------------------------------------------------------------------

$(ARPSTAR): $(ROOTALL) $(ARPSALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(ARPSALL)

$(AGRTAR): $(ROOTALL) $(ARPSALL) $(AGRALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(ARPSALL) $(AGRALL)

$(ASSIMTAR): $(ROOTALL) $(ARPSALL) $(ASSIMALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(ARPSALL) $(ASSIMALL)

$(CVTTAR): $(ROOTALL) $(LIBARPSALL) $(CVTALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(CVTALL)

$(CVTOBSTAR): $(ROOTALL) $(CVTOBSALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(CVTOBSALL)

$(BUDGTAR): $(ROOTALL) $(LIBARPSALL) $(BUDGALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(BUDGALL)

$(INTRPTAR): $(ROOTALL) $(LIBARPSALL) $(INTRPALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(INTRPALL)

$(TINTRPTAR): $(ROOTALL) $(LIBARPSALL) $(TINTRPALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(TINTRPALL)

$(TRAJCTAR): $(ROOTALL) $(LIBARPSALL) $(TRAJCALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(TRAJCALL)

$(PLTTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(PLTALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                       $(PLTALL)

$(PRTTAR): $(ROOTALL) $(LIBARPSALL) $(PRTALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(PRTALL)

$(DIFFTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(DIFFALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                       $(DIFFALL)

$(ENSTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(ENSALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                       $(ENSALL)

$(ADASTAR): $(ROOTALL) $(LIBARPSALL) $(ADASALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(ADASALL)

$(3DVARTAR): $(ROOTALL) $(LIBARPSALL) $(ADASALL) $(3DVARALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(ADASALL) $(3DVARALL)

$(ENKFTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(ENKFALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(ENKFALL)

$(DIFOBSTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(DIFOBSALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                       $(DIFOBSALL)

$(SFCTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(SFCALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                       $(SFCALL)

$(SOILTAR): $(ROOTALL) $(LIBARPSALL) $(SOILALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(SOILALL)

$(TRNTAR): $(ROOTALL) $(LIBARPSALL) $(TRNALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(TRNALL)

$(E2ATAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(E2AALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                             $(E2AALL)

$(A2WTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(A2WALL)         \
           $(LIBWRFINTALL) $(LIBPHDF5ALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                             $(A2WALL)  $(LIBPHDF5ALL) $(LIBWRFINTALL)

$(W2ATAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(W2AALL)         \
           $(LIBWRFINTALL) $(LIBPHDF5ALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                             $(W2AALL)  $(LIBPHDF5ALL) $(LIBWRFINTALL)

$(AVNTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(AVNALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                             $(AVNALL)

$(88DTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(88DALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                             $(88DALL)

$(ARPS2RADTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(ARPS2RADALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                             $(ARPS2RADALL)

$(RADMOSTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(RADMOSALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                             $(RADMOSALL)

$(RADEMULTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(RADEMULALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                             $(RADEMULALL)

$(RADB2CDFTAR): $(ROOTALL) $(LIBARPSALL) $(RADB2CDFALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) \
                             $(RADB2CDFALL)

$(RADSECTTAR): $(ROOTALL) $(LIBARPSALL) $(RADSECTALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) \
                             $(RADSECTALL)

$(ATTENTAR): $(ROOTALL) $(LIBARPSALL) $(ATTENALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) \
                             $(ATTENALL)

$(MOSAICTAR): $(ROOTALL) $(LIBARPSALL)  $(MOSAICALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) \
                             $(MOSAICALL)

$(MCITAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(MCIALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                             $(MCIALL) $(AGRINC) $(SFCINC)

$(EXTSNDTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(EXTSNDALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                             $(EXTSNDALL)

$(SKEWTTAR): $(ROOTALL) $(ZXPLOTALL) $(SKEWTALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(ZXPLOTALL) $(SKEWTALL)

$(ARPSMPTAR): $(ROOTALL) $(ARPSALL) $(MPALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(ARPSALL) $(MPALL)

$(READTAR): $(ROOTALL) $(LIBARPSALL) $(READSRC)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(READSRC)

$(A2GTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(A2GALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                       $(A2GALL)

$(A2NTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(A2NALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                       $(A2NALL)

$(A2ETAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(A2EALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                       $(A2EALL)

$(MERGETRNTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(MERGETRNSRC)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                       $(MERGETRNSRC)

$(ARPSVERTAR): $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) $(VERIFALL)
	@$(TAR) $(TAROPT) $@ $(ROOTALL) $(LIBARPSALL) $(LIBADASALL) \
                       $(VERIFALL)

$(ZXPLOTTAR) : $(ZXPLOTALL)
	@$(TAR) $(TAROPT) $@ $(ZXPLOTALL)

$(MLSOTAR): $(MLSOALL)
	@$(TAR) $(TAROPT) $@ $(MLSOALL)

#-----------------------------------------------------------------------
#
# Remove the object code for individual programs
#
#-----------------------------------------------------------------------

clean: clean.exe clean.arps clean.solver clean.arpsread clean.libarps \
       clean.arpscvt clean.arpsdiff clean.arpsplt clean.arpsintrp     \
       clean.arpstintrp clean.arpstrn clean.arpsprt clean.arpssfc     \
       clean.arpstern clean.dirtern clean.ext2arps clean.extract_avn  \
       clean.arpsextsnd clean.adas clean.arps2gem clean.arps3dvar     \
       clean.mergetrn clean.arpsagr clean.arpssoil clean.wtretcol     \
       clean.split clean.join clean.joins clean.libadas clean.88d2arps \
       clean.arps2rad clean.joinhdf clean.splithdf clean.joinbin2hdf \
       clean.tar clean.mci2arps clean.difobs clean.arpsens clean.arpsraindiff \
       clean.arpsassim clean.arps2ncdf clean.arpsverif clean.libradtn \
       clean.arpscvtobs clean.arps2wrf clean.wrf2arps  clean.mosaic   \
       clean.libzxpost clean.libzxncar clean.skewt clean.arpstrajc    \
       clean.ncrad2arps clean.atten clean.hdf2grads clean.radsector   \
       clean.radbin2cdf clean.radmosaic clean.radaremul clean.arps2eta212 \
       clean.arpsenkf clean.initensmbl clean.ossedata clean.postinnov \
       clean.endgns clean.libenkf clean.meso2lso clean.casa2arps      \
       clean.arps4wrf clean.nmm2arps clean.obsstd clean.arpsmpbudget  \
       clean.test clean.prdlib clean.arps2coamps clean.coamps2arps    \
       clean.tmatrix
	-$(RM) -f $(MODDIR)/*.mod

clean.exe: clean.spexe clean.mpexe clean.accessories

clean.spexe:
	-$(RM) -f $(BINDIR)/$(ARPSEXE)     $(BINDIR)/$(ADASEXE)     \
	          $(BINDIR)/$(E2AEXE)      $(BINDIR)/$(SPLITEXE)    \
	          $(BINDIR)/$(PLTNCAREXE)  $(BINDIR)/$(PLTPOSTEXE)  \
	          $(BINDIR)/$(INTRPEXE)    $(BINDIR)/$(INTRP_LSEXE) \
	          $(BINDIR)/$(TINTRPEXE)   $(BINDIR)/$(ARPSSFCEXE)  \
	          $(BINDIR)/$(ARPSSOILEXE) $(BINDIR)/$(NCRAD2AEXE)  \
	          $(BINDIR)/$(TRNEXE)      $(BINDIR)/$(MERGETRNEXE) \
	          $(BINDIR)/$(DIR1DEGEXE)  $(BINDIR)/$(DIR5MINEXE)  \
	          $(BINDIR)/$(DIR30SECEXE) $(BINDIR)/$(DIR1DEGEXE)  \
	          $(BINDIR)/$(88D2AEXE)    $(BINDIR)/$(88D2AEXE_FAKE) \
	          $(BINDIR)/$(JOINEXE)     $(BINDIR)/$(JOINSEXE)    \
	          $(BINDIR)/$(A2WEXE)      $(BINDIR)/$(W2AEXE)      \
	          $(BINDIR)/$(VERIFEXE)    $(BINDIR)/$(ENSICEXE)    \
                  $(BINDIR)/$(ENKFEXE)     $(BINDIR)/$(INITENSMBLEXE) \
                  $(BINDIR)/$(TINTRPMPEXE)

clean.mpexe:
	-$(RM) -f $(BINDIR)/$(ARPSMPEXE)    $(BINDIR)/$(ADASMPEXE)    \
	          $(BINDIR)/$(E2AMPEXE)     $(BINDIR)/$(SPLITMPEXE)   \
	          $(BINDIR)/$(PLTPOSTMPEXE) $(BINDIR)/$(PLTNCARMPEXE) \
	          $(BINDIR)/$(A2WMPEXE)     $(BINDIR)/$(W2AMPEXE)     \
	          $(BINDIR)/$(VERIFMPEXE)   $(BINDIR)/$(ENSICMPEXE)   \
	          $(BINDIR)/$(88D2AMPEXE)   $(BINDIR)/$(NCRAD2AMPEXE)

clean.accessories:
	-$(RM) -f $(BINDIR)/$(ARPSCVTEXE)  $(BINDIR)/$(ARPSCVTOBSEXE)  \
	          $(BINDIR)/$(A2NEXE)      $(BINDIR)/$(A2EEXE)      \
	          $(BINDIR)/$(PLTGRIDEXE)  $(BINDIR)/$(PLTMAXEXE)   \
	          $(BINDIR)/$(PRTEXE)      $(BINDIR)/$(A2GEXE)      \
	          $(BINDIR)/$(DIFFEXE)     $(BINDIR)/$(EXTSNDEXE)   \
	          $(BINDIR)/$(DIFOBSEXE)   $(BINDIR)/$(ARPS2RADEXE) \
	          $(BINDIR)/$(PLTRADCOLEXE) $(BINDIR)/$(RADMOSEXE)  \
	          $(BINDIR)/$(ARPSREADEXE) $(BINDIR)/$(WTRETCOLEXE) \
	          $(BINDIR)/$(TRAJCEXE)    $(BINDIR)/$(ASSIMEXE)    \
	          $(BINDIR)/$(AVNEXE)      $(BINDIR)/$(H2GEXE)      \
	          $(BINDIR)/arpscalctrajc  $(BINDIR)/joinbin2hdf    \
	          $(BINDIR)/$(RADARPLTNCAREXE) $(BINDIR)/$(RADARPLTPOSTEXE) \
	          $(BINDIR)/$(WSEXE)       $(BINDIR)/$(SUBDMNEXE)   \
	          $(BINDIR)/$(RADEMULEXE)  $(BINDIR)/$(RADB2CDFEXE) \
	          $(BINDIR)/$(RADSECTEXE)  $(BINDIR)/$(SPLITANYEXE) \
	          $(BINDIR)/$(WEXTSNDEXE)  $(BINDIR)/$(ATTENEXE)    \
	          $(BINDIR)/$(OSSEDATAEXE) $(BINDIR)/$(POSTINNOVEXE)\
	          $(BINDIR)/$(RNDPRTEXE)   $(BINDIR)/$(ENKFICEXE)   \
	          $(BINDIR)/$(RADARINTRPEXE) $(BINDIR)/$(A2WIIEXE)

clean.arps: clean.solver clean.libarps
	cd $(ARPSDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.arps

clean.solver:
	cd $(ARPSDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.solver

clean.libarps:
	cd $(ARPSDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.libarps

clean.libenkf:
	cd $(ENKFDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.libenkf

clean.prdlib:
	cd $(RDREMULDIR); make 'TOPDIR=$(TOPDIR)'   \
                               'BINDIR=$(BINDIR)'   \
                               'INCLDIR=$(INCLDIR)' \
                               'LIBDIR=$(LIBDIR)'   \
                               clean.prdlib

clean.arpsmpbudget:
	cd $(ARPSDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.arpsmpbudget

clean.libzxpost:
	cd $(ZXPLOTDIR); make 'TOPDIR=$(TOPDIR)'   \
                              'BINDIR=$(BINDIR)'   \
                              'INCLDIR=$(INCLDIR)' \
                              'LIBDIR=$(LIBDIR)'   \
                              clean.libzxpost
clean.libzxncar:
	cd $(ZXPLOTDIR); make 'TOPDIR=$(TOPDIR)'   \
                              'BINDIR=$(BINDIR)'   \
                              'INCLDIR=$(INCLDIR)' \
                              'LIBDIR=$(LIBDIR)'   \
                              clean.libzxncar

clean.libradtn:
	cd $(ARPSDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.libradtn

clean.adas: clean.libarps clean.libadas clean.libradtn
	cd $(ADASDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.adas

clean.arps3dvar: clean.libarps clean.libadas clean.libradtn
	cd $(3DVARDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.arps3dvar

clean.arpsenkf: clean.libarps clean.libenkf clean.libadas \
                clean.prdlib
	cd $(ENKFDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.arpsenkf

clean.libadas:
	cd $(ADASDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.libadas

clean.arpsagr: clean.solver
	cd $(AGRDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.arpsagr

clean.arpsassim: clean.solver
	cd $(ASSIMDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.arpsassim

clean.arpscvt: clean.libarps
	cd $(CVTDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.arpscvt

clean.hdf2grads:
	cd $(CVTDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.hdf2grads

clean.arpscvtobs: clean.libarps
	cd $(CVTOBS_DIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.arpscvtobs

clean.arpsplt: clean.libarps clean.libadas      \
               clean.libzxpost clean.libzxncar  \
               clean.prdlib
	cd $(PLTDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.arpsplt

clean.arpsprt: clean.libarps
	cd $(PRTDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.arpsprt

clean.arps2gem: clean.libarps clean.libadas
	cd $(E2ADIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'ARPSLIBDIR=$(LIBDIR)'\
                           clean.arps2gem

clean.arpsdiff: clean.libarps clean.libadas
	cd $(DIFDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.arpsdiff

clean.arpsraindiff: clean.libarps
	cd $(DIFDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.arpsraindiff

clean.arpsens: clean.libarps clean.libadas
	cd $(ENSDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.arpsens

clean.arpsextsnd: clean.libarps clean.libadas
	cd $(EXTSNDIR); make 'TOPDIR=$(TOPDIR)'   \
                             'BINDIR=$(BINDIR)'   \
                             'INCLDIR=$(INCLDIR)' \
                             'LIBDIR=$(LIBDIR)'   \
                             clean.arpsextsnd

clean.skewt: clean.libzxncar clean.libzxpost
	cd $(SKEWTDIR); make 'TOPDIR=$(TOPDIR)'   \
                             'BINDIR=$(BINDIR)'   \
                             'INCLDIR=$(INCLDIR)' \
                             'LIBDIR=$(LIBDIR)'   \
                             clean.skewt

clean.arpsintrp: clean.libarps
	cd $(INTRPDIR); make 'TOPDIR=$(TOPDIR)'   \
                             'BINDIR=$(BINDIR)'   \
                             'INCLDIR=$(INCLDIR)' \
                             'LIBDIR=$(LIBDIR)'   \
                             clean.arpsintrp

clean.arpstrajc: clean.libarps
	cd $(TRAJCDIR); make 'TOPDIR=$(TOPDIR)'   \
                             'BINDIR=$(BINDIR)'   \
                             'INCLDIR=$(INCLDIR)' \
                             'LIBDIR=$(LIBDIR)'   \
                             clean.arpstrajc


clean.arpstintrp: clean.libarps
	cd $(TINTRPDIR); make 'TOPDIR=$(TOPDIR)'   \
                             'BINDIR=$(BINDIR)'   \
                             'INCLDIR=$(INCLDIR)' \
                             'LIBDIR=$(LIBDIR)'   \
                             clean.arpstintrp

clean.arpssfc: clean.libarps clean.libadas
	cd $(SFCDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.arpssfc

clean.arpssoil: clean.libarps
	cd $(SOILDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.arpssoil

clean.arpstern: clean.libarps
	cd $(TERNDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.arpstern

clean.arpstrn: clean.libarps
	cd $(TRNDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.arpstrn

clean.mergetrn: clean.libarps clean.libadas
	cd $(TRNDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.mergetrn

clean.arps2ncdf: clean.libarps clean.libadas
	cd $(E2ADIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'ARPSLIBDIR=$(LIBDIR)'\
                           clean.arps2ncdf

clean.arps2eta212: clean.libarps clean.libadas
	cd $(E2ADIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'ARPSLIBDIR=$(LIBDIR)'\
                           clean.arps2eta212

clean.arps4wrf: clean.libarps
	cd $(A4WDIR); make 'TOPDIR=$(TOPDIR)'     \
                           'BINDIR=$(BINDIR)'     \
                           'INCLDIR=$(INCLDIR)'   \
                           'ARPSLIBDIR=$(LIBDIR)' \
                           clean.arps4wrf
	cd $(WRFINTDIR); make $(MAKEOPT)                  \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           superclean
	cd $(WRFNETDIR); make $(MAKEOPT)                   \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           'PHDF5PATH=$(PHDF5PATH)'       \
                           superclean

clean.nmm2arps: clean.libarps
	cd $(N2ADIR); make 'TOPDIR=$(TOPDIR)'     \
                           'BINDIR=$(BINDIR)'     \
                           'INCLDIR=$(INCLDIR)'   \
                           'ARPSLIBDIR=$(LIBDIR)' \
                           clean.nmm2arps
	cd $(WRFINTDIR); make $(MAKEOPT)                  \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           superclean
	cd $(WRFNETDIR); make $(MAKEOPT)                   \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           'PHDF5PATH=$(PHDF5PATH)'       \
                           superclean

clean.arps2coamps: clean.libarps
	cd $(A2CDIR); make 'TOPDIR=$(TOPDIR)'     \
                           'BINDIR=$(BINDIR)'     \
                           'INCLDIR=$(INCLDIR)'   \
                           'ARPSLIBDIR=$(LIBDIR)' \
                           clean.arps2coamps

clean.coamps2arps: clean.libarps
	cd $(C2ADIR); make 'TOPDIR=$(TOPDIR)'     \
                           'BINDIR=$(BINDIR)'     \
                           'INCLDIR=$(INCLDIR)'   \
                           'ARPSLIBDIR=$(LIBDIR)' \
                           clean.coamps2arps

clean.arps2wrf: clean.libarps clean.libadas
	cd $(A2WDIR); make 'TOPDIR=$(TOPDIR)'     \
                           'BINDIR=$(BINDIR)'     \
                           'INCLDIR=$(INCLDIR)'   \
                           'ARPSLIBDIR=$(LIBDIR)' \
                           clean.arps2wrf
	cd $(WRFINTDIR); make $(MAKEOPT)                  \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           superclean
	cd $(PHDF5DIR); make $(MAKEOPT)                   \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           'PHDF5PATH=$(PHDF5PATH)'       \
                           superclean

clean.wrf2arps: clean.libarps clean.libadas
	cd $(W2ADIR); make 'TOPDIR=$(TOPDIR)'     \
                           'BINDIR=$(BINDIR)'     \
                           'INCLDIR=$(INCLDIR)'   \
                           'ARPSLIBDIR=$(LIBDIR)' \
                           clean.wrf2arps
	cd $(WRFINTDIR); make $(MAKEOPT)                  \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           superclean
	cd $(PHDF5DIR); make $(MAKEOPT)                   \
                           'TOPDIR=$(TOPDIR)'             \
                           'LIBDIR=$(LIBDIR)'             \
                           'ARFLAG=$(ARFLAG)'             \
                           'PHDF5PATH=$(PHDF5PATH)'       \
                           superclean

clean.dirtern: clean.libarps
	cd $(TERNDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.dirtern

clean.ext2arps: clean.libarps clean.libadas clean.libg2
	cd $(E2ADIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'ARPSLIBDIR=$(LIBDIR)'\
                           clean.ext2arps

clean.extract_avn: clean.libarps clean.libadas
	cd $(AVNDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'APRSLIBDIR=$(LIBDIR)'\
                           clean.extract_avn

clean.88d2arps: clean.libarps clean.libadas
	cd $(88DDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.88d2arps

clean.ncrad2arps: clean.libarps clean.libadas
	cd $(88DDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.ncrad2arps

clean.casa2arps: clean.libarps clean.libadas
	cd $(88DDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.casa2arps

clean.arps2rad: clean.libarps clean.libadas
	cd $(88DDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.arps2rad

clean.radmosaic: clean.libarps clean.libadas
	cd $(ADASDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.radmosaic

clean.mosaic: clean.libarps
	cd $(MOSAICDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean

clean.tmatrix:
	cd $(RDREMULDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean

clean.mci2arps: clean.libarps clean.libadas
	cd $(MCIDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.mci2arps

clean.arpsread: clean.libarps
	cd $(ARPSDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.arpsread

clean.wtretcol: clean.libarps
	cd $(WTRETCOLDIR); make 'TOPDIR=$(TOPDIR)'   \
                                'BINDIR=$(BINDIR)'   \
                                'INCLDIR=$(INCLDIR)' \
                                'LIBDIR=$(LIBDIR)'   \
                                clean.wtretcol

clean.split: clean.libarps
	cd $(MPDIR); make 'TOPDIR=$(TOPDIR)'   \
                          'BINDIR=$(BINDIR)'   \
                          'INCLDIR=$(INCLDIR)' \
                          'LIBDIR=$(LIBDIR)'   \
                          clean.split

clean.join: clean.libarps
	cd $(MPDIR); make 'TOPDIR=$(TOPDIR)'   \
                          'BINDIR=$(BINDIR)'   \
                          'INCLDIR=$(INCLDIR)' \
                          'LIBDIR=$(LIBDIR)'   \
                          clean.join

clean.joins: clean.libarps
	cd $(MPDIR); make 'TOPDIR=$(TOPDIR)'   \
                          'BINDIR=$(BINDIR)'   \
                          'INCLDIR=$(INCLDIR)' \
                          'LIBDIR=$(LIBDIR)'   \
                          clean.joins

clean.joinbin2hdf: clean.libarps
	cd $(MPDIR); make 'TOPDIR=$(TOPDIR)'   \
                          'BINDIR=$(BINDIR)'   \
                          'INCLDIR=$(INCLDIR)' \
                          'LIBDIR=$(LIBDIR)'   \
                          clean.joinbin2hdf

clean.joinhdf: clean.libarps
	cd $(MPDIR); make 'TOPDIR=$(TOPDIR)'   \
                          'BINDIR=$(BINDIR)'   \
                          'INCLDIR=$(INCLDIR)' \
                          'LIBDIR=$(LIBDIR)'   \
                          clean.joinhdf
clean.splithdf: clean.libarps
	cd $(MPDIR); make 'TOPDIR=$(TOPDIR)'   \
                          'BINDIR=$(BINDIR)'   \
                          'INCLDIR=$(INCLDIR)' \
                          'LIBDIR=$(LIBDIR)'   \
                          clean.splithdf

clean.arpsverif: clean.libarps clean.libadas
	cd $(VERIFDIR); make 'TOPDIR=$(TOPDIR)'  \
                          'BINDIR=$(BINDIR)'     \
                          'INCLDIR=$(INCLDIR)'   \
                          'LIBDIR=$(LIBDIR)'     \
                          clean.arpsverif

clean.difobs: clean.libarps clean.libadas
	cd $(ADASDIR); make 'TOPDIR=$(TOPDIR)'  \
                          'BINDIR=$(BINDIR)'     \
                          'INCLDIR=$(INCLDIR)'   \
                          'LIBDIR=$(LIBDIR)'     \
                          clean.difobs

clean.initensmbl: clean.libarps clean.libadas clean.libenkf
	cd $(ENKFDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.initensmbl

clean.ossedata: clean.libarps clean.libadas clean.libenkf
	cd $(ENKFDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.ossedata

clean.obsstd  : clean.libenkf
	cd $(ENKFDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.obsstd

clean.postinnov:
	cd $(ENKFDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.postinnov

clean.endgns: clean.libarps clean.libadas clean.libenkf
	cd $(ENKFDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.endgns

clean.test:
	cd $(ENKFDIR); make 'TOPDIR=$(TOPDIR)'   \
                            'BINDIR=$(BINDIR)'   \
                            'INCLDIR=$(INCLDIR)' \
                            'LIBDIR=$(LIBDIR)'   \
                            clean.test

clean.libg2:
	cd $(G2DIR); make 'TOPDIR=$(TOPDIR)'   \
                          'BINDIR=$(BINDIR)'   \
                          'INCLDIR=$(INCLDIR)' \
                          'LIBDIR=$(LIBDIR)'   \
                          clean
clean.radaremul: clean.libarps clean.libadas
	cd $(88DDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.radaremul

clean.atten: clean.libarps
	cd $(88DDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.atten

clean.radbin2cdf: clean.libarps
	cd $(88DDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.radbin2cdf

clean.radsector: clean.libarps
	cd $(88DDIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean.radsector

clean.meso2lso:
	cd $(MLSODIR); make 'TOPDIR=$(TOPDIR)'   \
                           'BINDIR=$(BINDIR)'   \
                           'INCLDIR=$(INCLDIR)' \
                           'LIBDIR=$(LIBDIR)'   \
                           clean

clean.tar:
	-$(RM) -f $(ALLFILESTAR) $(ARPSTAR) $(AGRTAR) $(CVTTAR)    \
                  $(PLTTAR) $(PRTTAR) $(DIFFTAR) $(ADASTAR)        \
                  $(SFCTAR) $(SOILTAR) $(TRNTAR) $(EXTSNDTAR)      \
                  $(ARPSMPTAR) $(READTAR) $(A2GTAR) $(MERGETRNTAR) \
                  $(88DTAR) $(MCITAR) $(ENSTAR) $(A2NTAR)          \
                  $(TINTRPTAR) $(INTRPTAR) $(TERNTAR)              \
                  $(ARPSVERTAR) $(TRAJCTAR) $(MOSAICTAR)           \
                  $(ENKFTAR)
