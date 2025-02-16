#!/usr/bin/perl -w
##############################################################################
#
#  Test ARPS using a real case data on May 25, 1998
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
#          -agri	option to run arpsagri instead of arps
#          -arpssfc	option to run arpssfc for surface data
#          -arpstern	option to run arpstern for terrain
#          -ext2arps	option to run ext2arps to convert ETA data into ARPS grid
#          -adas	option to run adas for the initial data set
#          -data	option to run additional tests for data processing
#          -arpstrn     option to run arpstrn for terrain data\n";
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

$nx_glob = 67;
$ny_glob = 67;
$nz_glob = 43;

$nx = $nx_glob;
$ny = $ny_glob;
$nz = $nz_glob;

$runname = "arps25may1998";
$adasname = "adas25may1998";
$exbcname = "eta25may1998";
$execmd = "arps";

$curdir = cwd();
$subdir = "realcase";

$topdir = "$ENV{'ARPSDIR'}";
$arpsversion = "$ENV{'ARPSVER'}";

&GetOptions( "h", => \$help,
             "v", => \$verbose,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir,
             "d=s", => \$datadir,
             "agri", => \$agri,
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
  $arpstrn = 1;
  $ext2arps = 1;
  $adas     = 1;
}

if ( $agri ) {
  $subdir = "agri";
  $execmd = "arpsagri";
  $runname = "agri25may1998";
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

$tdatdir = "$datadir/data.test";

$arpsinput = "$srcdir/input/arps.input";

$trndata = "trn32km${nx}x${ny}.trndata";
$sfcdata = "sfc32km${nx}x${ny}.sfcdata";

if ($io eq "hdf") {
  $inifile  = "$adasname.hdf000000";
  $inigbf   = "$adasname.hdfgrdbas";
} elsif ($io eq "net") {
  $inifile  = "$adasname.net000000";
  $inigbf   = "$adasname.netgrdbas";
} else {
  $inifile  = "$adasname.bin000000";
  $inigbf   = "$adasname.bingrdbas";
}
$soilvar = "$exbcname.soilvar.000000";

if ( $verbose ) {
  print "\n";
  print "Current directory:   $curdir\n";
  print "Source directory:    $srcdir\n";
  print "Working directory:   $wrkdir\n";
  print "Data directory:      $datadir\n";
  print "Test data directory: $tdatdir\n";
  print "Terrain data file:   $trndata\n";
  print "Surface data file:   $sfcdata\n";
  print "Soil variable file:  $soilvar\n";
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
} else {
  ( -r "$wrkdir/$trndata" )
    or cp ( "$tdatdir/$trndata", "$wrkdir/$trndata" ) == 1
    or die "cp() failed";
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
  if ($io eq "hdf"){
    $inifile = "$exbcname.hdf000000" if ( ! -r "$tdatdir/$inifile" );
    $inigbf  = "$exbcname.hdfgrdbas" if ( ! -r "$tdatdir/$inigbf" );
  } elsif ($io eq "net"){
    $inifile = "$exbcname.net000000" if ( ! -r "$tdatdir/$inifile" );
    $inigbf  = "$exbcname.netgrdbas" if ( ! -r "$tdatdir/$inigbf" );
  } else {
    $inifile = "$exbcname.bin000000" if ( ! -r "$tdatdir/$inifile" );
    $inigbf  = "$exbcname.bingrdbas" if ( ! -r "$tdatdir/$inigbf" );
  }
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

( -r "$wrkdir/$sfcdata" ) or die "$0: $wrkdir/$sfcdata not found";
( -r "$wrkdir/$trndata" ) or die "$0: $wrkdir/$trndata not found";
( -r "$wrkdir/$exbcname.19980525.000000" )
   or die "$0: $wrkdir/$exbcname.19980525.000000 not found";
( -r "$wrkdir/$exbcname.19980525.030000" )
   or die "$0: $wrkdir/$exbcname.19980525.030000 not found";
( -r "$wrkdir/$exbcname.19980525.060000" )
   or die "$0: $wrkdir/$exbcname.19980525.060000 not found";
( -r "$wrkdir/$exbcname.soilvar.000000" )
   or die "$0: $wrkdir/$exbcname.soilvar.000000 not found";
( -r "$wrkdir/$inifile" ) or die "$0: $wrkdir/$inifile not found";
( -r "$wrkdir/$inigbf"  ) or die "$0: $wrkdir/$inigbf not found";

#-----------------------------------------------------------------
#
# make the executable arps
#
#-----------------------------------------------------------------

print "Making arps...\n" if ( $verbose );
chdir $srcdir;

# clean all objects if necessary
#system("makearps clean") ==0 or die "$0: makearps clean failed";

$makecmd = "./makearps  -m linux $run -opt $opt -io $io $execmd";
system($makecmd) == 0 or die "$0: $makecmd failed";

#----------------------------------------------------------------
#
# set the input parameters for arps
#
#----------------------------------------------------------------

$rstime = 10800.0;
$rstime_str = sprintf "%6.6d", $rstime;

%newpar = (
  nx => $nx,
  ny => $ny,
  nz => $nz,
 
  runold        => qq('$runname'),
  rstime        => $rstime,
  nxc           => $nx,
  nyc           => $ny,
  nzc           => $nz,
  grdsrt        => qq(.true.),
  nfinelv       => 1,
  "ngrdnew(1)"  => 1,
  "ixc(1,1)"    => 34.0,
  "jyc(1,1)"    => 32.5,
  "ixln(1,1)"   => 22.0,
  "jyln(1,1)"   => 21.0,
  "gangle(1,1)" => 0.0,

  runmod  => 1,

  initime => qq('1998-05-25.00:00:00'),
  initopt => 3,
  inifmt  => 1,
  inifile => qq('$inifile'),
  inigbf  => qq('$inigbf'),
  rstinf  => qq('$runname.rst$rstime_str'),

  tstop   => 21600.0,
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
  trbvimp  => 1,

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
  thisdmp => 3600.0,
  trstout => $rstime,

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
} elsif ($io eq "net") {
  $newpar{hdmpfmt} = 7;
  $newpar{inifmt}  = 7;
  $newpar{exbcfmt} = 7;
  $newpar{soilfmt} = 7;
  $newpar{sfcfmt}  = 7;
  $newpar{ternfmt} = 7;
} else {
  $newpar{hdmpfmt} = 1;
  $newpar{inifmt}  = 1;
  $newpar{exbcfmt} = 1;
  $newpar{soilfmt} = 1;
  $newpar{sfcfmt}  = 1;
  $newpar{ternfmt} = 1;
}

########################################################################
#
#  1. ARPS run
#
########################################################################

$newpar{tinitebd} = $newpar{initime};
$newpar{runname} = qq('$runname');

$newinput = "$wrkdir/$runname.input";
mkinput( \%newpar, $arpsinput, $newinput );

chdir $wrkdir;

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $runname started at $time";

$runcmd = "time $srcdir/bin/$execmd < $runname.input > $runname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $runname finished at $time\n";

########################################################################
#
#  2. A restart run
#
########################################################################

if ( $agri ) {
  $newpar{rstart} = qq(.true.);
}
else {
  $newpar{initopt} = 2;
}

$runname  = "${runname}_rst";
$newpar{runname} = qq("$runname");

$newinput = "$wrkdir/$runname.input";
mkinput( \%newpar, $arpsinput, $newinput );

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Restart run $runname started at $time";

$runcmd = "time $srcdir/bin/$execmd < $runname.input > $runname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Restart run $runname finished at $time\n";

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
  print "  -agri        option to run arpsagri instead of arps\n";
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
