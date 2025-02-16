#!/usr/bin/perl -w
##############################################################################
#
#  Test ADAS using a real case data on May 25, 1998
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
#          -arpstern	option to run arpstern for terrain
#          -arpstrn     option to run arpstrn for terrain data\n";
#          -ext2arps	option to run ext2arps to convert ETA data into ARPS grid
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

#$nxcpu = 2;
#$nycpu = 2;

$nx = $nx_glob;
$ny = $ny_glob;
$nz = $nz_glob;

$runname = "adas25may1998";
$etaname = "eta25may1998";
$execmd = "adas";

$curdir = cwd();
$subdir = "adas";

$topdir = "$ENV{'ARPSDIR'}";
$arpsversion = "$ENV{'ARPSVER'}";

&GetOptions( "h", => \$help,
             "v", => \$verbose,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir,
             "d=s", => \$datadir,
             "arpstern", => \$arpstern,
             "ext2arps", => \$ext2arps,
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

$tdatdir = "$datadir/data.test";

$scrptdir = "$srcdir/scripts";

$lapsdir = "$datadir/data.25may1998/lapsprd";
$nidsdir = "$tdatdir/nids/nids2arps";
$satdir  = "$datadir/data.25may1998/sat";


#$raddata  = "radsite.1998052500";
$mesodata1 = "okmeso/okmeso199805250000.lso";
$mesodata2 = "okmeso/okmeso199805250100.lso";
$saodata1 = "sao/sao199805250000.lso";
$saodata2 = "sao/sao199805250100.lso";
$snddata  = "raob/raob199805250000.snd";
#$prodata  = "pro/pro199805250000.pro";

$adasinput = "$srcdir/input/arps.input";
$runinput  = "$wrkdir/$runname.input";

$trndata  = "trn32km${nx}x${ny}.trndata";

if ( $verbose ) {
  print "\n";
  print "Current directory:   $curdir\n";
  print "Source directory:    $srcdir\n";
  print "Working directory:   $wrkdir\n";
  print "Data directory:      $datadir\n";
  print "Test data directory: $tdatdir\n";
  print "Terrain data file:   $trndata\n";
#  print "Surface data file:   $sfcdata\n";
#  print "Soil variable file:  $soilvar\n";
#  print "Base state file:     $inifile\n";
#  print "History data file:   $inigbf\n";
  print "\n";
}

#
# ARPSTERN/ARPSTRN
#

if ( $arpstern ) {
  print "\nCalling test_arpstern.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpstern.pl $verb $run -s $srcdir -w $wrkdir -d $datadir";
  system( $syscmd ) == 0 or die "$0: test_arpstern.pl failed";
}elsif ( $arpstrn ) {
  print "\nCalling test_arpstrn.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpstrn.pl $verb $run -s $srcdir -w $wrkdir -d $datadir -io $io -opt $opt";
  system( $syscmd ) == 0 or die "$0: test_arpstrn.pl failed";
} else {
  ( ! -f "$tdatdir/$trndata" )
    or cp ( "$tdatdir/$trndata", "$wrkdir/$trndata" ) == 1
    or die "cp() failed";
}

( -r $trndata ) or die "$0: $trndata not found";


#
# EXT2ARPS
#

if ( $ext2arps ) {
  print "\nCalling test_ext2arps.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_ext2arps.pl $verb $run -s $srcdir -w $wrkdir -d $datadir -io $io -opt $opt";
  system( $syscmd ) == 0 or die "$0: test_ext2arps.pl failed";
} else {

  if ($io eq "hdf") {
    ( -r "$wrkdir/$etaname.hdf000000" )
      or cp ( "$tdatdir/$etaname.hdf000000","$wrkdir/$etaname.hdf000000") == 1
      or die "cp() $tdatdir/$etaname.hdf000000 $wrkdir/$etaname.hdf000000 failed";
  
    ( -r "$wrkdir/$etaname.hdfgrdbas" )
      or cp ( "$tdatdir/$etaname.hdfgrdbas","$wrkdir/$etaname.hdfgrdbas") == 1
      or die "cp() failed";

  } elsif ($io eq "net") {

    ( -r "$wrkdir/$etaname.net000000" )
      or cp ( "$tdatdir/$etaname.net000000","$wrkdir/$etaname.net000000") == 1
      or die "cp() $tdatdir/$etaname.net000000 $wrkdir/$etaname.net000000 failed";
  
    ( -r "$wrkdir/$etaname.netgrdbas" )
      or cp ( "$tdatdir/$etaname.netgrdbas","$wrkdir/$etaname.netgrdbas") == 1
      or die "cp() failed";

  } else {
    ( -r "$wrkdir/$etaname.bin000000" )
      or cp ( "$tdatdir/$etaname.bin000000","$wrkdir/$etaname.bin000000") == 1
      or die "cp() $tdatdir/$etaname.bin000000 $wrkdir/$etaname.bin000000 failed";
  
    ( -r "$wrkdir/$etaname.bingrdbas" )
      or cp ( "$tdatdir/$etaname.bingrdbas","$wrkdir/$etaname.bingrdbas") == 1
      or die "cp() failed";
  }
}

if ($io eq "hdf") {
  ( -r "$wrkdir/$etaname.hdf000000" )
     or die "$0: $wrkdir/$etaname.hdf000000 not found";
  ( -r "$wrkdir/$etaname.hdfgrdbas" )
     or die "$0: $wrkdir/$etaname.hdfgrdbas not found";
} elsif ($io eq "net") {
  ( -r "$wrkdir/$etaname.net000000" )
     or die "$0: $wrkdir/$etaname.net000000 not found";
  ( -r "$wrkdir/$etaname.netgrdbas" )
     or die "$0: $wrkdir/$etaname.netgrdbas not found";
} else {
  ( -r "$wrkdir/$etaname.bin000000" )
     or die "$0: $wrkdir/$etaname.bin000000 not found";
  ( -r "$wrkdir/$etaname.bingrdbas" )
     or die "$0: $wrkdir/$etaname.bingrdbas not found";
}

chdir $srcdir;

#----------------------------------------------------------------
#
# Making the executable ADAS
#
#----------------------------------------------------------------

# clean all objects if necessary
#system("makearps clean") ==0 or die "$0: makearps clean failed";

$makecmd = "./makearps -m linux  $run -opt $opt -io $io $execmd";
system($makecmd) == 0 or die "$0: $makecmd failed";

#----------------------------------------------------------------
#
# set the input parameters for ADAS
#
#----------------------------------------------------------------

%newpar = (
   nx => $nx,
   ny => $ny,
   nz => $nz,
  
  runname => qq('$runname'),
  initime => q('1998-05-25.00:00:00'),
  initopt => 3,
  inifmt  => 1,
  inifile => qq('$etaname.bin000000'),
  inigbf  => qq('$etaname.bingrdbas'),

  dx      => 32000.0,
  dy      => 32000.0,
  dz      =>   500.0,
  dzmin   =>    20.0,
  ctrlat  =>    38.0,
  ctrlon  =>   -98.0,

  mapproj => 2,
  trulat1 => 30.0,
  trulat2 => 60.0,
  trulon  => -100.0,

  ternopt => 2,
  strhopt => 2,
  terndta => qq('$trndata'),

  hdmpfmt => 1,

  "dirname" => qq('./'),
  exbcdmp => 0,
  varout  => 1,
  mstout  => 1,
  sfcout  => 1,

  obropt  => 11,
  obrzero => 15000.0,

  "ianxtyp(3)" => 21,

  "xyrange(2)" => 100.E03,
  "xyrange(3)" => 80.E03,
  "xyrange(4)" => 80.E03,

  "zrange(1)" => 400.,
  "zrange(2)" => 400.,
  "zrange(3)" => 300.,
  "zrange(4)" => 200.,

  raduvobs  => 0,
  radrhobs  => 0,
  radcldopt => 0,
  radqvopt  => 0,
  radqcopt  => 0,
  radqropt  => 0,
  radptopt  => 0,

  range_cld => 100.0e03,

  rh_thr2  => 0.95,
  smth_opt => 1,

  ir_fname => qq('$satdir/mci2arps/sat25may1998.9805280015.goes08.cttemp'),
  vis_fname => qq('$satdir/mci2arps/sat25may1998.9805280015.goes08.albedo'),
  calib_fname => qq('$srcdir/data/adas/ircalib.adastab'),
  backerrfil => qq('$srcdir/data/adas/ruc3herr.adastab'),

  nsngfil => 2,
  "sngfname(1)" => qq('$lapsdir/$saodata1'),
  "sngtmchk(1)" => qq('$lapsdir/$saodata2'),
  "sngfname(2)" => qq('$lapsdir/$mesodata1'),
  "sngtmchk(2)" => qq('$lapsdir/$mesodata2'),

  blackfil => qq('$srcdir/data/adas/blacklist.sfc'),
  "srcsng(1)" => qq('SA'),
  "srcsng(2)" => qq('MESO'),
  "srcsng(3)" => qq('MOBLM'),
  "srcsng(4)" => qq('ARMMN'),
  "srcsng(5)" => qq('SYNOP'),
  "srcsng(6)" => qq('SHIP'),
  "srcsng(7)" => qq('BUOY'),
  "sngerrfil(1)" => qq('$srcdir/data/adas/saoerr.adastab'),
  "sngerrfil(2)" => qq('$srcdir/data/adas/mesoerr.adastab'),
  "sngerrfil(3)" => qq('$srcdir/data/adas/moblmerr.adastab'),
  "sngerrfil(4)" => qq('$srcdir/data/adas/armmnerr.adastab'),
  "sngerrfil(5)" => qq('$srcdir/data/adas/synoperr.adastab'),
  "sngerrfil(6)" => qq('$srcdir/data/adas/shiperr.adastab'),
  "sngerrfil(7)" => qq('$srcdir/data/adas/buoyerr.adastab'),

  nuafil => 1,
  "uafname(1)" => qq('$lapsdir/$snddata'),
  "srcua(1)" => qq('NWS RAOB'),
  "srcua(2)" => qq('FCLASS'),
  "srcua(3)" => qq('MCLASS'),
  "srcua(4)" => qq('WPDN PRO'),
  "uaerrfil(1)" => qq('$srcdir/data/adas/raoberr.adastab'),
  "uaerrfil(2)" => qq('$srcdir/data/adas/fclasserr.adastab'),
  "uaerrfil(3)" => qq('$srcdir/data/adas/mclasserr.adastab'),
  "uaerrfil(4)" => qq('$srcdir/data/adas/profilerr.adastab'),

# skip nids2arps nradfil => 0, otherwise nradfil => 2
  nradfil => 0,
  "srcrad(1)" => qq('88D-NIDS'),
  "srcrad(2)" => qq('88D-NIDS'),
  "radfname(1)"  => qq('$nidsdir/KAMA.980525.0001'),
  "radfname(2)"  => qq('$nidsdir/KARX.980525.0000'),
  "raderrfil(1)" => qq('$srcdir/data/adas/radnidserr.adastab'),
  "raderrfil(2)" => qq('$srcdir/data/adas/radnidserr.adastab'),
  "raderrfil(3)" => qq('$srcdir/data/adas/radcasaerr.adastab'),
  "iuserad(1,4)" => 1,
  "iuserad(2,4)" => 1,

  "srcret(1)"    => qq('88D-RET'),
  "reterrfil(1)" => qq('$srcdir/data/adas/ret88Derr.adastab'),
);

if ($io eq "hdf") {
  $newpar{hdmpfmt} = 3;
  $newpar{inifmt}  = 3;
  $newpar{ternfmt} = 3;
  $newpar{inifile} = qq('$etaname.hdf000000'),
  $newpar{inigbf}  = qq('$etaname.hdfgrdbas'),
} elsif ($io eq "net") {
  $newpar{hdmpfmt} = 7;
  $newpar{inifmt}  = 7;
  $newpar{ternfmt} = 7;
  $newpar{inifile} = qq('$etaname.net000000'),
  $newpar{inigbf}  = qq('$etaname.netgrdbas'),
} else {
  $newpar{hdmpfmt} = 1;
  $newpar{inifmt}  = 1;
  $newpar{ternfmt} = 1;
  $newpar{inifile} = qq('$etaname.bin000000'),
  $newpar{inigbf}  = qq('$etaname.bingrdbas'),
}

print "Making the input namelist file...\n";
mkinput( \%newpar, $adasinput, $runinput );


##############################################################
#
# Run ADAS
#
##############################################################

if ( $newpar{ternopt} == 2 ) {
  print "Copying ARPS terrain data file...\n";
  ( ! -f "$tdatdir/$trndata" ) 
    or cp( "$tdatdir/$trndata", "$wrkdir/$trndata" ) == 1
    or die "cp() failed";
}

chdir $wrkdir;

$runcmd = "time $srcdir/bin/$execmd < $runname.input > $runname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

print "ADAS was successful.\n";

#----------------------------------------------------------------
#
# Copy file to data directory
#
#----------------------------------------------------------------

print "Copying data file to data directory...\n";

( ! -f "$wrkdir/$runname.gribcntl" )
    or cp( "$wrkdir/$runname.gribcntl", "$tdatdir/$runname.gribcntl" ) == 1
    or die "cp() failed";
  
if ($io eq "hdf") {
  ( ! -f "$wrkdir/$runname.hdfgrdbas" )
    or cp( "$wrkdir/$runname.hdfgrdbas", "$tdatdir/$runname.hdfgrdbas" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.hdf000000" )
    or cp( "$wrkdir/$runname.hdf000000", "$tdatdir/$runname.hdf000000" ) == 1
    or die "cp() failed";
} elsif ($io eq "net") {
  ( ! -f "$wrkdir/$runname.netgrdbas" )
    or cp( "$wrkdir/$runname.netgrdbas", "$tdatdir/$runname.netgrdbas" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.net000000" )
    or cp( "$wrkdir/$runname.net000000", "$tdatdir/$runname.net000000" ) == 1
    or die "cp() failed";
} else {
  ( ! -f "$wrkdir/$runname.bingrdbas" )
    or cp( "$wrkdir/$runname.bingrdbas", "$tdatdir/$runname.bingrdbas" ) == 1
    or die "cp() failed";
  
  ( ! -f "$wrkdir/$runname.bin000000" )
    or cp( "$wrkdir/$runname.bin000000", "$tdatdir/$runname.bin000000" ) == 1
    or die "cp() failed";
}

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
  print "  -d dir       data directory\n";
  print "  -arpstern    option to run arpstern for terrain\n";
  print "  -ext2arps    option to run ext2arps to convert ETA data into ARPS grid\n";
  print "  -arpstrn     option to run arpstrn for terrain data\n";
  print "  -io hdf/bin  option for history dump format\n";
  print "  -opt 0-4     Compiler optimization level\n";
  print "\n";
  exit 1;
}
