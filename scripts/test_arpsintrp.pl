#!/usr/bin/perl -w

use Cwd;
use File::Copy "cp";
use File::Basename;
use Getopt::Long;

use lib "./scripts";
use mkinput;

#$nx_glob = 67;
#$ny_glob = 67;
#$nz_glob = 43;

$nx = 67;
$ny = 67;
$nz = 43;

$runname = "intrp25may1998";
$execmd = "arpsintrp";

$curdir = cwd();
$subdir = "arpsintrp";

$topdir = "$ENV{'ARPSDIR'}";
$arpsversion = "$ENV{'ARPSVER'}";

goto end_of_usage;

usage: {
  print "\n";
  print "Usage: $0 [options]\n";
  print "\n";
  print "Options:\n";
  print "  -v           verbose\n";
  print "  -n           do not actually build and execute the program\n";
  print "  -w dir       working directory\n";
  print "  -s dir       ARPS source directory\n";
  print "  -d dir       data directory\n";
  print "\n";
  exit 1;
}

end_of_usage:

&GetOptions( "v", => \$verbose,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir,
             "d=s", => \$datadir);

$verb = '';
$run  = '';
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

$intrpinput = "$srcdir/input/arpsintrp.input";
$runinput  = "$wrkdir/$runname.input";

$trndata  = "trn32km${nx}x${ny}.trndata";
$sfcdata  = "sfc32km${nx}x${ny}.sfcdata";

$inifile = "arps25may1998.hdf000000";
$inigbf  = "arps25may1998.hdfgrdbas";

if ( $verbose ) {
  print "\n";
  print "Current directory:   $curdir\n";
  print "Source directory:    $srcdir\n";
  print "Working directory:   $wrkdir\n";
  print "Data directory:      $datadir\n";
  print "Terrain data file:   $trndata\n";
  print "Surface data file:   $sfcdata\n";
  print "Base state file:     $inigbf\n";
  print "History data file:   $inifile\n";
  print "\n";
}

( ! -f "$tdatdir/$trndata" )
    or cp ( "$tdatdir/$trndata", "$wrkdir/$trndata" ) == 1
    or die "cp() failed";

( -r $trndata ) or die "$0: $trndata not found";

( ! -f "$tdatdir/$sfcdata" )
    or cp ( "$tdatdir/$sfcdata", "$wrkdir/$sfcdata" ) == 1
    or die "cp() failed";

( -r $sfcdata ) or die "$0: $sfcdata not found";


( -r "$wrkdir/$inifile" )
   or die "$0: $wrkdir/$inifile not found";
( -r "$wrkdir/$inigbf" )
   or die "$0: $wrkdir/$inigbf not found";

chdir $srcdir;

# clean all objects if necessary
#system("makearps clean") ==0 or die "$0: makearps clean failed";

$makecmd = "makearps $run -opt 3 $execmd";
system($makecmd) == 0 or die "$0: $makecmd failed";

# set the input parameters for arpsintrp

%newpar = (
   hinfmt      => 3,
   hdmpinopt   => 2,
   grdbasfn    => $inigbf,
   nhisfilei   => 1,
   "hisfile(1)"=> $inifile,

# set dimensional parameters for ARPS
   nx1 => $nx,
   ny1 => $ny,
   nz1 => $nz,
  
  runname => qq('$runname'),

  z_or_p   => 1,
  xy_or_ll => 2,
  ctrlat1  => 36.0,
  ctrlon1  => -100,

  snap_to_grid => 0,
  same_res     => 0,

  dx1  =>  30000,
  dy1  =>  30000,
  dz1  =>    400,

  zrefsfc1  =>  0,
  strhopt1  =>  1,
  dzmin1    =>  20,
  dlayer11  =>  0,
  dlayer21  =>  100000,
  strhtune1 =>  0.2,
  zflat1    =>  100000,
  nstyp1    =>  9,
  noutgrds  =>  0,

  ternopt1  => 2,
  ntmerge   => 12,
  mntopt1   => 0,
  terndta1  => $trndata,
  ternfmt1  => 0,

  intrphopt => 2,
  intrpvopt => 1,

  bglopt    => 4,
  misvalue  => -9999.0,

  ntagopt => 1,
  aghght  => 1000.0,

  sfcdat  => 0,

  "dirname"  => './',
  hdmpfmt    =>  3,
  grbpkbit   => 16,
  hdfcompr   =>  5,

  filcmprs => 0,
  readyfl  => 1,

  basout   => 0,
  grdout   => 0,
  varout   => 1,
  mstout   => 1,
  iceout   => 1,
  tkeout   => 0,
  trbout   => 0,
  rainout  => 0,
  prcout   => 0,

  sfcout=>1,
  landout  => 1,
  radout   => 0,
  flxout   => 0,

  exbcdmp=>0,
     qcexout  => 1,
     qrexout  => 1,
     qiexout  => 1,
     qsexout  => 1,
     qhexout  => 1,

  ngbrz=>12,
     zbrdmp => 12000.0,

  sfcdmp  => 0,
  soildmp=>2,
  terndmp => 0
);

print "Making the input namelist file...\n";
mkinput( \%newpar, $intrpinput, $runinput );

if ( $newpar{ternopt} == 2 ) {
  print "Copying ARPS terrain data file...\n";
  ( ! -f "$tdatdir/$trndata" ) 
    or cp( "$tdatdir/$trndata", "$wrkdir/$trndata" ) == 1
    or die "cp() failed";
}

chdir $wrkdir;

$runcmd = "time $srcdir/bin/$execmd < $runname.input > $runname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

print "ARPSINTRP was successful.\n";

exit 0
