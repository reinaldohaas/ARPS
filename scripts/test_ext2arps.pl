#!/usr/bin/perl -w
##############################################################################
#
#  Test EXT2ARPS using a real case data on May 25, 1998
#
#-----------------------------------------------------------------------------
#
#  USAGE: test_adas.pl [options]
#
#  Options:
#          -h           show this usage
#          -v           verbose
#          -n           do not actually build and execute the program
#          -w dir	working directory
#          -s dir	ARPS source directory
#          -arpstern	option to run arpstern for terrain
#          -arpstrn     option to run arpstrn for terrain data\n";
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

$runname = "eta25may1998";
$execmd = "ext2arps";

$curdir = cwd();
$topdir = "$ENV{'ARPSDIR'}";
$arpsversion = "$ENV{'ARPSVER'}";
$subdir = "ext2arps";


&GetOptions( "h", => \$help,
             "v", => \$verbose,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir,
             "d=s", => \$datadir,
             "arpstern", => \$arpstern,
	     "arpstrn",  => \$arpstrn,
	     "io=s",     => \$io,
	     "opt=s",    => \$opt);

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

$scrptdir = "$srcdir/scripts";

if ($ENV{'HOST'} eq 'paige' || $ENV{'HOST'} eq 'cirrus') {
  $etadir  = "/work/official/data.25may1998/eta40grb";        # paige
} else {
  $etadir  = "$datadir/data.25may1998/eta40grb";             # when $datadir is not too long
}

$tdatdir = "$datadir/data.test";

if( $verbose ) {
  print "\n";
  print "Current directory: $curdir\n";
  print "Source directory: $srcdir\n";
  print "Working directory: $wrkdir\n";
  print "Data directory: $datadir\n";
  print "\n";
}

$arpsinput = "$srcdir/input/arps.input";
$e2ainput = "$wrkdir/$runname.input";

$trndata = "trn32km${nx}x${ny}.trndata";

#
# ARPSTERN/ARPSTRN
#

if ( $arpstern ) {
  print "\nCalling test_arpstern.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpstern.pl $verb $run -s $srcdir -w $wrkdir -d $datadir";
  system( $syscmd ) == 0 or die "$0: test_arpstern.pl failed";
} elsif ( $arpstrn) {
  print "\nCalling test_arpstrn.pl ...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpstrn.pl $verb $run -s $srcdir -w $wrkdir -d $datadir -io $io -opt $opt";
  system( $syscmd ) == 0 or die "$0: test_arpstrn.pl failed";
} elsif ( ! -r $trndata ) {
  ( ! -r "$tdatdir/$trndata" )
    or cp ( "$tdatdir/$trndata", "$wrkdir/$trndata" ) == 1
    or die "cp() failed";
}

( -r $trndata ) or die "$0: $trndata not found";


#----------------------------------------------------------------
#
# set the input parameters for ext2arps
#
#----------------------------------------------------------------

%newpar = (
  nx => $nx,
  ny => $ny,
  nz => $nz,
  
  runmod  => 1,
  runname => qq('$runname'),
  initime => qq('1998-05-25.00:00:00'),

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

  ternopt => 2,
  terndta => qq('$trndata'),

  sfcdmp  => 1,
  soildmp => 1,    # Add sfcdmp and soildmp
  sfcphy  => 0,

  hdmpfmt => 1,
  thisdmp => 3600.0,
  trstout => 3600.0,

  "dirname" => qq('./'),
  exbcdmp  => 1,
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

  dir_extd    => qq('$etadir'),
  extdopt     => 2,
  extdname    => qq('/eta40grb'),
  nextdfil    => 3,
  "extdtime(1)" => qq('1998-05-25.00:00:00+000:00:00'),
  "extdtime(2)" => qq('1998-05-25.00:00:00+003:00:00'),
  "extdtime(3)" => qq('1998-05-25.00:00:00+006:00:00'),
);

if ($io eq "hdf") {
  $newpar{hdmpfmt} = 3;
  $newpar{exbcdmp} = 3;
  $newpar{soildmp} = 3;
  $newpar{ternfmt} = 3;
} elsif ($io eq "net") {
  $newpar{hdmpfmt} = 7;
  $newpar{exbcdmp} = 7;
  $newpar{soildmp} = 7;
  $newpar{ternfmt} = 7;
} else {
  $newpar{hdmpfmt} = 1;
  $newpar{exbcdmp} = 1;
  $newpar{soildmp} = 1;
  $newpar{ternfmt} = 1;
}

$newpar{tinitebd} = $newpar{initime};
print "Making the input namelist file...\n" if ( $verbose );
mkinput( \%newpar, $arpsinput, $e2ainput );

#-------------------------------------------------------------
#
#  Making EXT2ARPS
#
#-------------------------------------------------------------

chdir $srcdir;

# clean all objects if necessary
#system("makearps clean") ==0 or die "$0: makearps clean failed";

$makecmd = "./makearps -m linux $run -opt $opt -io $io $execmd";
system($makecmd) == 0 or die "$0: $makecmd failed";

##############################################################
#
# Run EXT2ARPS
#
##############################################################

chdir $wrkdir;

$runcmd = "$srcdir/bin/$execmd < $runname.input > $runname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

print "Ext2arps was successful\n";

#--------------------------------------------------------------
#
# Copy data to data directory
#
#--------------------------------------------------------------

print "Copying data file to data directory...\n";

#( ! -f "$wrkdir/$runname.gridcntl" ) 
#  or cp( "$wrkdir/$runname.gridcntl", "$tdatdir/$runname.gridcntl" ) == 1
#  or die "cp() failed";

if ($io eq "hdf") {

  ( ! -f "$wrkdir/$runname.hdfgrdbas" )
    or cp( "$wrkdir/$runname.hdfgrdbas", "$tdatdir/$runname.hdfgrdbas" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.hdf000000" )
    or cp( "$wrkdir/$runname.hdf000000", "$tdatdir/$runname.hdf000000" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.hdf010800" )
    or cp( "$wrkdir/$runname.hdf010800", "$tdatdir/$runname.hdf010800" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.hdf021600" )
    or cp( "$wrkdir/$runname.hdf021600", "$tdatdir/$runname.hdf021600" ) == 1
    or die "cp() failed";
  
} elsif ($io eq "net") {

  ( ! -f "$wrkdir/$runname.netgrdbas" )
    or cp( "$wrkdir/$runname.netgrdbas", "$tdatdir/$runname.netgrdbas" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.net000000" )
    or cp( "$wrkdir/$runname.net000000", "$tdatdir/$runname.net000000" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.net010800" )
    or cp( "$wrkdir/$runname.net010800", "$tdatdir/$runname.net010800" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.net021600" )
    or cp( "$wrkdir/$runname.net021600", "$tdatdir/$runname.net021600" ) == 1
    or die "cp() failed";
  
} else {
  
  ( ! -f "$wrkdir/$runname.bingrdbas" )
    or cp( "$wrkdir/$runname.bingrdbas", "$tdatdir/$runname.bingrdbas" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.bin000000" )
    or cp( "$wrkdir/$runname.bin000000", "$tdatdir/$runname.bin000000" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.bin010800" )
    or cp( "$wrkdir/$runname.bin010800", "$tdatdir/$runname.bin010800" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.bin021600" )
    or cp( "$wrkdir/$runname.bin021600", "$tdatdir/$runname.bin021600" ) == 1
    or die "cp() failed";

}

( ! -f "$wrkdir/$runname.soilvar.000000" )
    or cp( "$wrkdir/$runname.soilvar.000000", "$tdatdir/$runname.soilvar.000000" ) == 1
    or die "cp() failed";
  
( ! -f "$wrkdir/$runname.soilvar.010800" )
    or cp( "$wrkdir/$runname.soilvar.010800", "$tdatdir/$runname.soilvar.010800" ) == 1
    or die "cp() failed";
  
( ! -f "$wrkdir/$runname.soilvar.021600" )
    or cp( "$wrkdir/$runname.soilvar.021600", "$tdatdir/$runname.soilvar.021600" ) == 1
    or die "cp() failed";
  
( ! -f "$wrkdir/$runname.19980525.000000" )
    or cp( "$wrkdir/$runname.19980525.000000", "$tdatdir/$runname.19980525.000000" ) == 1
    or die "cp() failed";
  
( ! -f "$wrkdir/$runname.19980525.030000" )
    or cp( "$wrkdir/$runname.19980525.030000", "$tdatdir/$runname.19980525.030000" ) == 1
    or die "cp() failed";
  
( ! -f "$wrkdir/$runname.19980525.060000" )
    or cp( "$wrkdir/$runname.19980525.060000", "$tdatdir/$runname.19980525.060000" ) == 1
    or die "cp() failed";

exit 0;

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
  print "  -arpstern    option to run arpstern for terrain data\n";
  print "  -arpstrn     option to run arpstrn for terrain data\n";
  print "  -io hdf/bin  option for history dump format\n";
  print "  -opt 0-4     Compiler optimization level\n";
  print "\n";
  exit 1;
}
