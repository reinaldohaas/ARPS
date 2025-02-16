#!/usr/bin/perl -w
##############################################################################
#
#  Test ARPSTRN
#
#-----------------------------------------------------------------------------
#
#  USAGE: test_arpstrn.pl [options]
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
use File::Copy 'cp';
use File::Basename;
use Getopt::Long;

use lib "./scripts";
use mkinput;

$nx = 67;
$ny = 67;

$subdir = "arpstrn";
$execmd = "arpstrn";
$rname = "arpstrn";

$curdir = cwd();
$topdir = "$ENV{'ARPSDIR'}";
$arpsversion = "$ENV{'ARPSVER'}";

&GetOptions( "h", => \$help,
             "v", => \$verbose,
             "n", => \$norun,
             "w=s",   => \$wrkdir,
             "s=s",   => \$srcdir,
             "d=s",   => \$datadir,
             "io=s",  => \$io,
             "opt=s", => \$opt);

if ( $help ) {  &Usage();  }

( $io  ) or $io  = "bin";
( $opt ) or $opt = "3";

$verb= ''; $run = '';
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

$trndir  = "$datadir/arpstern.data";
$tdatdir = "$datadir/data.test";

if ( $verbose ) {
  print "\n";
  print "Current directory: $curdir\n";
  print "Source directory: $srcdir\n";
  print "Working directory: $wrkdir\n";
  print "Data directory: $datadir\n";
  print "\n";
}

$arpstrninput = "$srcdir/input/arpstrn.input";
$trninput = "$wrkdir/$rname.input";

$trndata = "trn32km${nx}x${ny}.trndata";

#--------------------------------------------------------
#
# set the input parameters for arpstrn
#
#--------------------------------------------------------

%newpar = (
  nx => $nx,
  ny => $ny,

  runname => "'trn32km${nx}x${ny}'",

  dx      => 32000.0,
  dy      => 32000.0,
  ctrlat  =>    38.0,
  ctrlon  =>   -98.0,

  mapproj => 2,
  trulat1 => 30.0,
  trulat2 => 60.0,
  trulon  => -100.0,
  
  trnanxopt => 1,
  trndataopt => 2,
  dir_trndata => "'$trndir'",
  fn_trndata => "'tbase_global_5min.data'",

  qw(mapfile(1)) => "'$trndir/us_state.mapdata'",
  qw(mapfile(2)) => "'$trndir/world_us_country.mapdata'"
);

if ($io eq "net") {
  $newpar{terndmp} = 7
} elsif ($io eq "hdf") {
  $newpar{terndmp} = 3
} else {
  $newpar{terndmp} = 1
}

print "Making the input namelist file...\n";
mkinput( \%newpar, $arpstrninput, $trninput );

#------------------------------------------------------------
#
# Make the executable
#
#------------------------------------------------------------

chdir $srcdir;

# Clean object files before makearps
#system("makearps clean") ==0 or die "$0: makearps clean failed";

$makecmd = "./makearps -m linux $run -opt $opt -io $io $execmd";
system($makecmd) == 0 or die "$0: $makecmd failed";

################################################################
#
# Run ARPSTRN
#
###############################################################

chdir $wrkdir;

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $rname started  at $time";

$runcmd = "$srcdir/bin/$execmd < $rname.input > $rname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $rname finished at $time";

#------------------------------------------------------------
#
# Copy data
#
#-------------------------------------------------------------

(! -f "$wrkdir/$trndata" )
  or cp( "$wrkdir/$trndata","$tdatdir/$trndata" ) == 1
  or die "$0: ARPSTRN test fail." ;

exit 0;

###################################################################
#
# subroutine Usage
#
###################################################################

sub Usage {
  print "\n";
  print "Usage: $0 [options]\n";
  print "\n";
  print "Options:\n";
  print "  -h           show this usage\n";
  print "  -v           verbose\n";
  print "  -s dir       ARPS source directory\n";
  print "  -w dir       working directory\n";
  print "  -d dir       data directory\n";
  print "  -io hdf/bin  option for history dump format\n";
  print "  -opt 0-4     Compiler optimization level\n";
  print "\n";
  exit 1;
}
