#!/usr/bin/perl -w
##############################################################################
#
#  Test ARPSSFC
#
#-----------------------------------------------------------------------------
#
#  USAGE: test_realcase.pl [options]
#
#  Options:
#          -h           show this usage
#          -v           verbose
#          -n           do not actually build and execute the program
#          -w dir	working directory
#          -s dir	ARPS source directory
#          -d dir       data directory
#          -io hdf/bin  option for history dump format\n";
#          -opt 0-4     Compiler optimization level\n";
#
#############################################################################

use Cwd;
use File::Basename;
use File::Copy "cp";
use Getopt::Long;

use lib "./scripts";
use mkinput;

$nx_glob = 67;
$ny_glob = 67;
$nz_glob = 43;

#$nxcpu = 2;
#$nycpu = 2;

$nx = $nx_glob;
$ny = $ny_glob;
$nz = $nz_glob;

$runname = "arpssfc";
$execmd = "arpssfc";
$subdir = "arpssfc";

$curdir = cwd();
$topdir = "$ENV{'ARPSDIR'}";
$arpsversion = "$ENV{'ARPSVER'}";


&GetOptions( "h", => \$help,
             "v", => \$verbose,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir,
             "d=s", => \$datadir,
	     "io=s",=> \$io,
	     "opt=s",=> \$opt);

if ( $help ) {  &Usage();  }

( $io  ) or $io  = "bin";
( $opt ) or $opt = "3";

$verb = ""; $run  = "";
( ! $verbose ) or $verb = "-v";
( ! $norun   ) or $run  = "-n";

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

$sfcdir  = "$datadir/arpssfc.data";

$tdatdir = "$datadir/data.test";

$arpsinput = "$srcdir/input/arps.input";
$sfcinput = "$wrkdir/$runname.input";

if ( $verbose ) {
  print "\n";
  print "Current directory: $curdir\n";
  print "Source directory: $srcdir\n";
  print "Working directory: $wrkdir\n";
  print "Data directory: $datadir\n";
  print "\n";
}

$sfcdata = "sfc32km${nx}x${ny}.sfcdata";

#-------------------------------------------------------------
#
# set the input parameters for ext2arps
#
#-------------------------------------------------------------

%newpar = (
   nx => $nx,
   ny => $ny,
   nz => $nz,
   
   runmod  => 1,
   runname => qq('$runname'),
   initime => q('1998-05-25.00:00:00'),

   tstop   => 21600.0,
   dtbig   => 40.0,
   dtsml   => 40.0,

   dx      => 32000.0,
   dy      => 32000.0,
   dz      =>   500.0,
   dzmin   =>    20.0,
   strhopt => 2,
   ctrlat  =>    38.0,
   ctrlon  =>   -98.0,

   mapproj => 2,
   
  trulat1 => 30.0,
  trulat2 => 60.0,
  trulon  => -100.0,

  sfcdat  => 2,
  sfcdtfl => qq('$sfcdata'),

  "dirname" => qq('./'),
  exbcdmp => 1,
  varout  => 1,
  mstout  => 1,
  iceout  => 0,
  tkeout  => 0,
  trbout  => 0,
  rainout => 0,
  sfcout  => 1,
  landout => 0,
  prcout  => 0,
  radout  => 0,
  flxout  => 0,
  qcexout => 0,
  qrexout => 0,
  qiexout => 0,
  qsexout => 0,
  qhexout => 0,

  schmopt => 3,
  sdatopt => 1,
  vdatopt => 1,
  ndatopt => 1,
  vfrcopt => 1,
  nstyp   => 4,
  stypout => 1,
  vtypout => 1,
  laiout  => 1,
  rfnsout => 1,
  vegout  => 1,
  ndviout => 1,
  drawval => 1,
  fstypfl => qq('$sfcdir/soil_1km.data'),
  fvtypfl => qq('$sfcdir/naoge1_01l_1km.img'),
  fndvifl => qq('$sfcdir/namay92ndl_1km.img'),
  bstypfl => qq('$sfcdir/whsoil_1deg.data'),
  bvtypfl => qq('$sfcdir/owe14d_10min.data'),
  bndvifl => qq('$sfcdir/ndvi9005_10min.data'),
  vfrcdr  => qq('$sfcdir'),
);

if ($io eq "hdf") {
  $newpar{sfcdmp} = 3;
} elsif ($io eq "net") {
  $newpar{sfcdmp} = 7;
} else {
  $newpar{sfcdmp} = 1;
}

print "Making the input namelist file...\n" if ( $verbose );
mkinput( \%newpar, $arpsinput, $sfcinput );

#---------------------------------------------------------------
#
# Make executable
#
#---------------------------------------------------------------

chdir $srcdir;

# clean all objects if necessary
#system("makearps clean") ==0 or die "$0: makearps clean failed";

$makecmd = "./makearps -m linux $run -opt $opt -io $io $execmd";
system($makecmd) == 0 or die "$0: $makecmd failed";

#################################################################
#
# Run ARPSSFC
#
#################################################################

chdir $wrkdir;

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $runname started  at $time";

$runcmd = "$srcdir/bin/$execmd < $runname.input > $runname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

#----------------------------------------------------------------
#
# Copy data
#
#----------------------------------------------------------------

( ! -r "arpssfc.sfcdata" )
  or cp( "arpssfc.sfcdata", "$sfcdata" )
  or die "$0: cp() failed";

( ! -f "$wrkdir/$sfcdata" )
  or cp( "$wrkdir/$sfcdata","$tdatdir/$sfcdata" ) == 1
  or die "$0: cp() failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $runname finished at $time";

exit 0;

sub Usage {
  print "\n";
  print "Usage: $0 [options]\n";
  print "\n";
  print "Options:\n";
  print "  -h           show this usage\n";
  print "  -v           verbose\n";
  print "  -n           do not actually build and execute the program\n";
  print "  -s dir       ARPS source directory\n";
  print "  -w dir       working directory\n";
  print "  -d dir       data directory\n";
  print "  -io bin/hdf  history dump format\n";
  print "  -opt 0-4     compiler optimization option\n";
  print "\n";
  exit 1;
}
