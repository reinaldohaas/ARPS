#!/usr/bin/perl -w
##############################################################################
#
#  Test ARPS2WRF using a real case data on May 25, 1998
#
#-----------------------------------------------------------------------------
#
#  USAGE: test_arps2wrf.pl [options]
#
#  Options:
#          -h           show this usage
#          -v           verbose
#          -n           do not actually build and execute the program
#          -w dir	working directory
#          -s dir	ARPS source directory
#          -arpssfc	  option to run arpssfc for surface data
#          -arpstrn   option to run arpstrn for terrain data
#          -arpstern	option to run arpstern for terrain
#          -ext2arps	option to run ext2arps to convert ETA data into ARPS grid
#          -adas	    option to run adas for the initial data set
#          -data	    option to run additional tests for data processing
#          -io hdf/bin  option for history dump format\n";
#          -opt 0-4     Compiler optimization level\n";
#
#############################################################################

use Cwd;
use File::Copy "cp";
use File::Basename;
use Getopt::Long;

use lib "./scripts";
use mkinput;

$nx     = 67;
$ny     = 67;

$nx_wrf = 65;
$ny_wrf = 65;

$runname = "a2w25may1998";
$adasname = "adas25may1998";
$exbcname = "eta25may1998";
$execmd = "arps2wrf";

$curdir = cwd();
$subdir = "arps2wrf";

$topdir = "$ENV{'ARPSDIR'}";
$arpsversion = "$ENV{'ARPSVER'}";

&GetOptions( "h", => \$help,
             "v", => \$verbose,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir,
             "d=s", => \$datadir,
             "arpssfc", => \$arpssfc,
             "arpstern", => \$arpstern,
             "ext2arps", => \$ext2arps,
             "adas", => \$adas,
             "data", => \$data,
             "arpstrn",  => \$arpstrn,
             "io=s",     => \$io,
             "opt=s",    => \$opt);

if ( $help ) {  &Usage();  }

if ( $data ) {
  $arpssfc  = 1;
  $arpstrn  = 1;
  $ext2arps = 1;
  $adas     = 1;
}

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

( -d "$srcdir/data/data.25may1998" ) or 
   die "ARPS data directory not found in $srcdir/data\n";

chdir $curdir;
chdir $srcdir;
$srcdir = cwd();

chdir $curdir;
chdir $wrkdir;
$wrkdir = cwd();

$scrptdir = "$srcdir/scripts";

$tdatdir  = "$datadir/data.test";

$extdir   = "$srcdir/ext2arps";

$arpsinput = "$srcdir/input/arps2wrf.input";

$trndata = "trn32km${nx}x${ny}.trndata";
$sfcdata = "sfc32km${nx}x${ny}.sfcdata";

if ($io eq "hdf") {
  $inifile  = "$adasname.hdf000000";
  $inigbf   = "$adasname.hdfgrdbas";
} else {
  $inifile  = "$adasname.bin000000";
  $inigbf   = "$adasname.bingrdbas";
}

if ( $verbose ) {
  print "\n";
  print "Current directory:   $curdir\n";
  print "Source directory:    $srcdir\n";
  print "Working directory:   $wrkdir\n";
  print "Data directory:      $datadir\n";
  print "Test data directory: $tdatdir\n";
  print "Terrain data file:   $trndata\n";
  print "Surface data file:   $sfcdata\n";
  print "Base state file:     $inifile\n";
  print "History data file:   $inigbf\n";
  print "\n";
}

#
# ARPSSFC
#

if ( $arpssfc ) {
  print "\nCalling test_arpssfc.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpssfc.pl $verb $run -s $srcdir -w $wrkdir -d $datadir -io $io -opt $opt";
  system( $syscmd ) ==0 or die "$0: test_arpssfc.pl failed";
} else {
  ( -r "$wrkdir/$sfcdata" )
    or cp ( "$tdatdir/$sfcdata", "$wrkdir/$sfcdata" ) == 1
    or die "cp() failed";
}

#
# ARPSTERN/ARPSTRN
#

if ( $arpstern ) {
  print "\nCalling test_arpstern.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpstern.pl $verb $run -s $srcdir -w $wrkdir -d $datadir";
  system( $syscmd ) ==0 or die "$0: test_arpstern.pl failed";
} elsif ( $arpstrn ) {
  print "\nCalling test_arpstrn.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpstrn.pl $verb $run -s $srcdir -w $wrkdir -d $datadir -io $io -opt $opt";
  system( $syscmd ) ==0 or die "$0: test_arpstrn.pl failed";
} 

#
# EXT2ARPS
#

if ( $ext2arps ) {
  print "\nCalling test_ext2arps.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_ext2arps.pl $verb $run -s $srcdir -w $wrkdir -d $datadir -io $io -opt $opt";
  system( $syscmd ) ==0 or die "$0: test_ext2arps.pl failed";
} else {
  ( -r "$wrkdir/$exbcname.bin000000" )
    or cp ( "$extdir/$exbcname.bin000000",
            "$wrkdir/$exbcname.bin000000" ) == 1
    or die "cp() failed";

  ( -r "$wrkdir/$exbcname.bingrdbas" )
    or cp ( "$extdir/$exbcname.bingrdbas",
            "$wrkdir/$exbcname.bingrdbas" ) == 1
    or die "cp() failed";

  ( -r "$wrkdir/$exbcname.bin010800" )
    or cp ( "$extdir/$exbcname.bin010800",
            "$wrkdir/$exbcname.bin010800" ) == 1
    or die "cp() failed";

  ( -r "$wrkdir/$exbcname.bingrdbas.01" )
    or cp ( "$extdir/$exbcname.bingrdbas.01",
            "$wrkdir/$exbcname.bingrdbas.01" ) == 1
    or die "cp() failed";

  ( -r "$wrkdir/$exbcname.bin021600" )
    or cp ( "$extdir/$exbcname.bin021600",
            "$wrkdir/$exbcname.bin021600" ) == 1
    or die "cp() failed";

  ( -r "$wrkdir/$exbcname.bingrdbas.02" )
    or cp ( "$extdir/$exbcname.bingrdbas.02",
            "$wrkdir/$exbcname.bingrdbas.02" ) == 1
    or die "cp() failed";

}

#
#   ADAS
#

if ( $adas ) {
  print "\nCalling test_adas.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_adas.pl $verb $run -s $srcdir -w $wrkdir -d $datadir -io $io -opt $opt";
  system( $syscmd ) ==0 or die "$0: test_adas.pl failed";
} else {
  $inifile = "$exbcname.bin000000" if ( ! -r "$tdatdir/$inifile" );
  $inigbf  = "$exbcname.bingrdbas" if ( ! -r "$tdatdir/$inigbf" );
  ( -r "$wrkdir/$inifile" )
    or cp ( "$tdatdir/$inifile","$wrkdir/$inifile" ) == 1
    or die "$0: cp() failed";
  ( -r "$wrkdir/$inigbf" )
    or cp ( "$tdatdir/$inigbf", "$wrkdir/$inigbf" ) == 1
    or die "$0: cp() failed";
}

#------------------------------------------------------------------
#
# Check the necessary files
#
#------------------------------------------------------------------

( -r "$wrkdir/$exbcname.bin000000" )
   or die "$0: $wrkdir/$exbcname.bin000000 not found";
( -r "$wrkdir/$exbcname.bingrdbas" )
   or die "$0: $wrkdir/$exbcname.bingrdbas not found";
( -r "$wrkdir/$exbcname.bin010800" )
   or die "$0: $wrkdir/$exbcname.bin010800 not found";
( -r "$wrkdir/$exbcname.bingrdbas.01" )
   or die "$0: $wrkdir/$exbcname.bingrdbas.01 not found";
( -r "$wrkdir/$exbcname.bin021600" )
   or die "$0: $wrkdir/$exbcname.bin021600 not found";
( -r "$wrkdir/$exbcname.bingrdbas.02" )
   or die "$0: $wrkdir/$exbcname.bingrdbas.02 not found";
( -r "$wrkdir/$sfcdata" ) or die "$0: $wrkdir/$sfcdata not found";
( -r "$wrkdir/$inifile" ) or die "$0: $wrkdir/$inifile not found";
( -r "$wrkdir/$inigbf"  ) or die "$0: $wrkdir/$inigbf not found";

#-----------------------------------------------------------------
#
# make the executable arps2wrf
#
#-----------------------------------------------------------------

print "Making arps2wrf ...\n" if ( $verbose );
chdir $srcdir;

# clean all objects if necessary
#system("makearps clean") ==0 or die "$0: makearps clean failed";

$makecmd = "./makearps -m linux $run -opt $opt $execmd";
system($makecmd) == 0 or die "$0: $makecmd failed";

#----------------------------------------------------------------
#
# set the input parameters for arps
#
#----------------------------------------------------------------

%newpar = (

  hinfmt           => 1,
  grdbasfn         => qq('$inigbf'),
  "hisfile(1)"     => qq('$inifile'),

  sfcinitopt       => qq('ARPS'),
  sfcfmt           => 1,
  sfcdtfl          => qq('$sfcdata'),

  create_bdy       => 1,
  bdyfheader       => qq('./$exbcname'),
  tbgn_bdyin       => 10800,
  tintv_bdyin      => 10800,
  tend_bdyin       => 21600,

  use_arps_grid    => 0,
  nx_wrf           => $nx_wrf,
  ny_wrf           => $ny_wrf,
  mapproj_wrf      => 2,
  sclfct_wrf       => 1.0,
  trulat1_wrf      => 30.0,
  trulat2_wrf      => 60.0,
  trulon_wrf       => -98.0,
  ctrlat_wrf       => 38.0,
  ctrlon_wrf       => -98.0,
  dx_wrf           => 30000.,
  dy_wrf           => 30000.,
 
  ptop             => 5000,


  dyn_opt          => 2,
  diff_opt         => 0,
  km_opt           => 1,
  khdif            => 0.0,
  kvdif            => 0.0,
  mp_physics       => 4,
  ra_lw_physics    => 1,
  ra_sw_physics    => 1,
  bl_sfclay_physics  => 1,
  bl_surface_physics => 2,
  bl_pbl_physics     => 1,
  cu_physics         => 1,

  dt                 => 40,
  spec_bdy_width     => 5,

  iorder             => 2,
  korder             => 2,

  dirname            => qq('./'),
  readyfl            => 0,
  create_namelist    => 1,
  wrfnamelist        => qq('namelist.input')
);

########################################################################
#
#  1. ARPS2WRF run
#
########################################################################

$newinput = "$wrkdir/$runname.input";
mkinput( \%newpar, $arpsinput, $newinput );

chdir $wrkdir;

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $runname started at $time";

$runcmd = "time $srcdir/bin/$execmd < $runname.input > $runname.output";
print $runcmd
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $runname finished at $time\n";

exit 0;

#############################################################################
#
#  SUB Usage
#
#############################################################################

sub Usage {
  print "\n";
  print "Usage: $0 [options]\n";
  print "\n";
  print "Options:\n";
  print "  -h           show this usage\n";
  print "  -v           verbose\n";
  print "  -n           do not actually build and execute the program\n";
  print "  -w dir       working directory\n";
  print "  -s dir       ARPS source directory\n";
  print "  -arpssfc     option to run arpssfc for surface data\n";
  print "  -arpstrn     option to run arpstrn for terrain data\n";
  print "  -arpstern    option to run arpstern for terrain\n";
  print "  -ext2arps    option to run ext2arps to convert ETA data into ARPS grid\n";
  print "  -adas        option to run adas for the initial data set\n";
  print "  -data        option to run additional tests for data processing\n";
  print "  -io hdf/bin  option for history dump format\n";
  print "  -opt 0-4     Compiler optimization level\n";
  print "\n";
  exit 1;
}
