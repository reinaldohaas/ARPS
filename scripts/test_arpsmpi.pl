#!/usr/bin/perl -w
######################################################################
#
# Usage:
#
#    test_arpsmpi.pl [options]
#
#     Options:
#       -v           verbose
#       -n           do not actually build and execute the program
#       -w dir	     working directory
#       -s dir	     ARPS source directory
#       -d dir	     Data directory
#       -io bin/hdf  history dump format
#       -opt 0-4     compiler optimization level
#       -real        option to run with real data
#
#####################################################################

use Cwd;
use File::Copy "cp";
use File::Basename;
use Getopt::Long;

use lib "./scripts";
use mkinput;

$nx = 67;
$ny = 67;
$nz = 43;

$nxcpu = 2;
$nycpu = 2;

$runname = "mpi25may1998";
$adasname = "adas25may1998";
$exbcname = "eta25may1998";
$execmd = "arps_mpi";
$execmd2 = "arps";

$curdir = cwd();
$subdir = "arpsmpi";

$topdir = "$ENV{'ARPSDIR'}";
$arpsversion = "$ENV{'ARPSVER'}";

&GetOptions( "h", => \$help,
             "v", => \$verbose,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir,
             "d=s", => \$datadir,
	     "io=s",  => \$io,
	     "opt=s", => \$opt,
             "real",  => \$real,
             "split",  => \$split );

&Usage if ($help);

( $io  ) or $io  = "bin";
( $opt ) or $opt = "3";

$verb=""; $run = "";
( ! $verbose ) or $verb = "-v";
( ! $norun )   or $run = "-n";

$real = 1;
( $split ) or $split = 0;

( $srcdir  ) or $srcdir  = "$topdir/arps${arpsversion}";
( $wrkdir  ) or $wrkdir  = "$srcdir/$subdir";
( $datadir ) or $datadir = "$srcdir/data";

( -d "$srcdir/src/arps" ) or die "ARPS root directory not found in $srcdir\n";
( -d $wrkdir ) or mkdir $wrkdir, 0755;

chdir $curdir;
chdir $srcdir;
$srcdir = cwd();

chdir $curdir;
chdir $wrkdir;
$wrkdir = cwd();

if ( $verbose ) {
  print "\n";
  print "Current directory: $curdir\n";
  print "Source directory: $srcdir\n";
  print "Working directory: $wrkdir\n";
  print "Data directory: $datadir\n";
  print "\n";
}

#$scrptdir = "$srcdir/scripts";

$tdatdir = "$datadir/data.test";

$arpsinput = "$srcdir/input/arps.input";
#$arpscvttempl = "$srcdir/input/arpscvt.input";

$trndata = "trn32km${nx}x${ny}.trndata";
$sfcdata = "sfc32km${nx}x${ny}.sfcdata";
if ($io eq "hdf") {
  $inifile  = "$adasname.hdf000000";
  $inigbf   = "$adasname.hdfgrdbas";
} else {
  $inifile  = "$adasname.bin000000";
  $inigbf   = "$adasname.bingrdbas";
}
$soilvar = "$exbcname.soilvar.000000";

#-----------------------------------------------------------
#
# set the input parameters for arps_mpi
#
#-----------------------------------------------------------

%newpar = (
  nx => $nx,
  ny => $ny,
  nz => $nz,
 
  nproc_x => $nxcpu,
  nproc_y => $nycpu,
  max_fopen => $nxcpu*$nycpu, 
  
  runmod  => 1,

  initime => qq('1998-05-25.00:00:00'),
  initopt => 3,
  inifmt  => 1,
  inisplit => 0,
  inifile => qq('$inifile'),
  inigbf  => qq('$inigbf'),

  tstop   => 120.0,
  dtbig   => 40.0,
  dtsml   => 40.0,

  dx      => 32000.0,
  dy      => 32000.0,
  dz      => 500.0,
  dzmin   => 20.0,
  strhopt => 2,
  ctrlat  =>    38.0,
  ctrlon  =>   -98.0,

  mapproj => 2,
  trulat1 => 30.0,
  trulat2 => 60.0,
  trulon  => -100.0,

  ternopt => 2,
  terndta => qq('$trndata'),

  lbcopt  => 2,
  wbc     => 5,
  ebc     => 5,
  nbc     => 5,
  sbc     => 5,

  exbcname => qq('$exbcname'),
  brlxhw => 2.8,
  cbcdmp => 0.001,
  cbcmix => 0.0001,

  madvopt => 3,
  sadvopt => 3,

  coriopt => 4,

  tmixopt  => 4,
  trbisotp => 0,
  tkeopt   => 3,

  cmix2nd => 1,
  cfcm2h  => 0.0,
  cfcm2v  => 0.0001,
  cmix4th => 1,
  cfcm4h  => 0.0001,
  cfcm4v  => 0.0,

  divdmp => 0,
  raydmp => 2,
  cfrdmp => 0.001,
  zbrdmp => 14000.0,

  mphyopt => 2,
  cnvctopt => 3,
  confrq   => 300.0,

  radopt  => 2,

  sfcphy  => 4,
  landwtr => 1,
  pbldopt => 2,
  sfcdat  => 2,
  sfcdtfl => qq('$sfcdata'),
  soilinit => 2,
  soilinfl => qq('$exbcname.soilvar.000000'),

  hdmpfmt => 1,
  dmp_out_joined => 1,
  thisdmp => 40.0,
  trstout => 0.0,

  "dirname" => qq('./'),
  exbcdmp => 0,
  varout  => 1,
  mstout  => 1,
  iceout  => 1,
  tkeout  => 1,
  trbout  => 1,
  rainout => 1,
  sfcout  => 1,
  landout => 1,
  prcout  => 1,
  radout  => 1,
  flxout  => 1,
); 

if ($io eq "hdf") {
  $newpar{hdmpfmt} = 3;
  $newpar{inifmt}  = 3;
  $newpar{exbcfmt} = 3;
  $newpar{soilfmt} = 3;
  $newpar{sfcfmt}  = 3;
  $newpar{ternfmt} = 3;
} else {
  $newpar{hdmpfmt} = 1;
  $newpar{inifmt}  = 1;
  $newpar{exbcfmt} = 1;
  $newpar{soilfmt} = 1;
  $newpar{sfcfmt}  = 1;
  $newpar{ternfmt} = 1;
}

if ($split) {
  $newpar{inisplited} =1;
  $newpar{dmp_out_joined} = 0;
}

$newpar{tinitebd} = $newpar{initime};
$newpar{runname} = qq('$runname');

$mpiinput = "$wrkdir/arps.input";
mkinput( \%newpar, $arpsinput, $mpiinput );
system ("cp $mpiinput $wrkdir/arps.input") if ($mpiinput ne "$wrkdir/arps.input");

$newpar{runname} = qq('no$runname');
$nompiinput = "$wrkdir/arps_nompi.input";
mkinput( \%newpar, $arpsinput, $nompiinput );

#=====================================================================
#
# set the input parameters for arpscvt
# DO NOT need any more
#
#=====================================================================

=ARPSCVT from grd to bin format, no need any more

%arpscvtin = (

  hinfmt => 10,
  nhisfile => 1,

  hdmpfmt => 1,
  thisdmp => 1200.0,
  trstout => 0.0,

  "dirname" => qq('./'),
  exbcdmp => 0,
  varout  => 1,
  mstout  => 1,
  iceout  => 1,
  tkeout  => 1,
  trbout  => 1,
  rainout => 1,
  sfcout  => 1,
  landout => 1,
  prcout  => 1,
  radout  => 1,
  flxout  => 1,
); 

$arpscvtin{"grdbasfn"} = qq('${adasname}.grbgrdbas'),
$arpscvtin{"hisfile(1)"} = qq('${adasname}.grb000000'),

$cvtinput = "$wrkdir/arpscvt.input";
mkinput( \%arpscvtin, $arpscvttempl, $cvtinput );

=cut end of ARPSCVT


#------------------------------------------------------------
#
# Copy data to the work directory
#
#------------------------------------------------------------

if ( $real ) {
  ( -r "$wrkdir/$sfcdata" )
    or cp ( "$tdatdir/$sfcdata", "$wrkdir/$sfcdata" ) == 1
    or die "cp() failed";

  ( -r "$wrkdir/$trndata" )
    or cp ( "$tdatdir/$trndata", "$wrkdir/$trndata" ) == 1
    or die "cp() failed";

#  $inifile = "$exbcname.grb000000" if ( ! -r "$tdatdir/$inifile" );
#  $inigbf  = "$exbcname.grbgrdbas" if ( ! -r "$tdatdir/$inigbf" );

  ( -r "$wrkdir/$exbcname.19980525.000000" )
    or cp ( "$tdatdir/$exbcname.19980525.000000",
	    "$wrkdir/$exbcname.19980525.000000" ) == 1
    or die "cp() failed";

  ( -r "$wrkdir/$exbcname.19980525.030000" )
    or cp ( "$tdatdir/$exbcname.19980525.030000",
	    "$wrkdir/$exbcname.19980525.030000" ) == 1
    or die "cp() failed";

  ( -r "$wrkdir/$exbcname.19980525.060000" )
    or cp ( "$tdatdir/$exbcname.19980525.060000",
	    "$wrkdir/$exbcname.19980525.060000" ) == 1
    or die "cp() failed";

  ( -r "$wrkdir/$soilvar" )
    or cp ( "$tdatdir/$soilvar", "$wrkdir/$soilvar" )  == 1
    or die "cp() failed";
  
  if ($io eq "hdf") {
    ( -r "$wrkdir/$adasname.hdf000000" )
      or cp ( "$tdatdir/$adasname.hdf000000","$wrkdir/$adasname.hdf000000" ) == 1
      or die "$0: cp() failed";

    ( -r "$wrkdir/$adasname.hdfgrdbas" )
      or cp ( "$tdatdir/$adasname.hdfgrdbas", "$wrkdir/$adasname.hdfgrdbas" ) == 1
      or die "$0: cp() failed";
  
  } elsif ($io eq "net") {
    ( -r "$wrkdir/$adasname.net000000" )
      or cp ( "$tdatdir/$adasname.net000000","$wrkdir/$adasname.net000000" ) == 1
      or die "$0: cp() failed";

    ( -r "$wrkdir/$adasname.netgrdbas" )
      or cp ( "$tdatdir/$adasname.netgrdbas", "$wrkdir/$adasname.netgrdbas" ) == 1
      or die "$0: cp() failed";
  
  } else {
    ( -r "$wrkdir/$adasname.bin000000" )
      or cp ( "$tdatdir/$adasname.bin000000","$wrkdir/$adasname.bin000000" ) == 1
      or die "$0: cp() failed";

    ( -r "$wrkdir/$adasname.bingrdbas" )
      or cp ( "$tdatdir/$adasname.bingrdbas", "$wrkdir/$adasname.bingrdbas" ) == 1
      or die "$0: cp() failed";
  }
}

#-----------------------------------------------------------------
#
# make the arps_mpi executables
#
#-----------------------------------------------------------------

print "Making arps_mpi...\n" if ( $verbose );
chdir $srcdir;

if ($split) {
  $makecmd = "./makearps -m linux $run -opt $opt splitfiles";
  system( $makecmd ) ==0 or die "$0: $makecmd failed";

  $makecmd = "./makearps -m linux $run -opt $opt joinfiles";
  system( $makecmd ) ==0 or die "$0: $makecmd failed";
}

# make the executable arps

print "Making arps...\n" if ( $verbose );

$makecmd = "./makearps -m linux $run -opt $opt -io $io $execmd2";
system($makecmd) == 0 or die "$0: $makecmd failed";

$makecmd = "./makearps -m linux $run -opt $opt -io $io $execmd";
system($makecmd) == 0 or die "$0: $makecmd failed";

#=========================================================
#
# make arpscvt
# DO NOT need any more
#
#=========================================================

=ARPSCVT

$makecmd = "./makearps -m linux $run -opt 3 arpscvt";
system($makecmd) == 0 or die "$0: $makecmd failed";

chdir $wrkdir;

# run arpscvt

unless (-s $inifile && -s $inigbf 
     && (-M "$adasname.grb000000" > -M $inifile)
     && (-M "$adasname.grbgrdbas" > -M $inigbf)) {

     print `rm -f $inifile $inigbf`;
     $syscmd = "$srcdir/bin/arpscvt < $cvtinput > arpscvt.output";
     system( $syscmd ) ==0 or die "$0: $syscmd failed";
}
=cut ARPSCVT

#####################################################################
#
#  1. MPI RUN
#
####################################################################

chdir $wrkdir;

if ($split) {
#------------------------------------------------------------------
#  1.1  Splitting files
#------------------------------------------------------------------

  print "Splitting files...\n";
  $syscmd = "$srcdir/bin/splitfiles < $mpiinput > $runname.split.output";
  system( $syscmd ) ==0 or die "$0: $syscmd failed";
}

#------------------------------------------------------------------
#  1.2  Running MPI 
#------------------------------------------------------------------

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $runname started at $time";

$runcmd = 
   "time mpirun -np 4 $srcdir/bin/$execmd < arps.input > $runname.output";
unless ($norun) {
   print 'teste',$runcmd;
  $mpistat = system( $runcmd );
  if ($mpistat != 0) {
    print "$0: WARNING $runcmd returned status $mpistat ($?)\n";
  }
}

#------------------------------------------------------------------
#  1.3 Joining files 
#------------------------------------------------------------------

if ( $real && $split ) {
  print "Joining files...\n";
  $syscmd = "$srcdir/bin/joinfiles < $mpiinput > $runname.join.output";
  $joinstat = system( $syscmd );
  if ($joinstat != 0) {
    print "$0: WARNING $syscmd failed status $joinstat ($?)\n";
  }
}

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $runname finished at $time\n";

#####################################################################
#
#  2. Regular run to compare against (should be bit-for-bit the same)
#
####################################################################

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run no$runname started at $time";

$runcmd = "$srcdir/bin/$execmd2 < $nompiinput > no$runname.output";
unless ($norun) {
  $nompistat = system( $runcmd );
  if ($nompistat != 0) {
    print "$0: WARNING $runcmd failed status $nompistat ($?)\n";
  }
}

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run no$runname finished at $time\n";

exit 0;

####################################################################
#
# Sub Usage
#
####################################################################

sub Usage {
  print "\n";
  print "Usage: $0 [options]\n";
  print "\n";
  print "Options:\n";
  print "  -v           verbose\n";
  print "  -n           do not actually build and execute the program\n";
  print "  -w dir	working directory\n";
  print "  -s dir	ARPS source directory\n";
  print "  -d dir	Data directory\n";
  print "  -io bin/hdf  history dump format\n";
  print "  -opt 0-4     compiler optimization option\n";
  print "  -real        option to run with real data\n";
  print "  -split       option to run splitfiles and joinfiles\n";
  print "\n";
  exit 1;
}
