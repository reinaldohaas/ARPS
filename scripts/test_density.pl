#!/usr/bin/perl -w

use Cwd;
use File::Basename;
use File::Copy "cp";
use Getopt::Long;

use lib "./scripts";
use mkinput;

$nx_glob = 256;
$ny_glob = 04;
$nz_glob = 64;

#$nxcpu = 2;
#$nycpu = 2;

$nx = $nx_glob;
$ny = $ny_glob;
$nz = $nz_glob;

$runname = "density";
$subdir = "density";
$execmd = "arps";

$curdir = cwd();
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
  print "\n";
  exit 1;
}

end_of_usage:

&GetOptions( "v", => \$verbose,
             "n", => \$norun,
             "s=s", => \$srcdir,
             "w=s", => \$wrkdir );

$verb = "";
$run  = "";
( ! $verbose ) or $verb = "-v";
( ! $norun   ) or $run  = "-n";

( $srcdir ) or $srcdir = "$topdir/arps${arpsversion}";
( $wrkdir ) or $wrkdir = "$srcdir/$subdir";

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
  print "\n";
}

cp( "$srcdir/scripts/density_grads.r", "$wrkdir/density_grads.r" );

$arpsinput = "$srcdir/input/arps.input";

# make the executable arps if necessary

unless(-x "$srcdir/bin/$execmd"){

  print "Making arps...\n" if ( $verbose );
  chdir $srcdir;

  $makecmd = "makearps $run -opt 3 $execmd";
  system($makecmd) == 0 or die "$0: $makecmd failed";
}

# set the input parameters for ext2arps

%newpar = (
# set dimensional parameters for ARPS
   nx => $nx,
   ny => $ny,
   nz => $nz,

  runmod  => 2,
  inibasopt => 2,
  ubar0   => 0.0,
  vbar0   => 0.0,
  pt0opt  => 6,
  ptpert0 => -15.0,
  pt0radx => 4000.0,
  pt0rady => 0.0,
  pt0radz => 2000.0,
  pt0ctrx => 0.0,
  pt0ctry => 0.0,
  pt0ctrz => 3000.0,

  tstop   => 900.0,
  dtbig   => 1.0,
  dtsml   => 0.25,

  dx      => 100.0,
  dy      => 100.0,
  dz      => 100.0,
  dzmin   => 100.0,

  wbc     => 1,
  ebc     => 1,
  nbc     => 2,
  sbc     => 2,
  tbc     => 1,
  bbc     => 1,

  madvopt => 1,
  sadvopt => 1,

  tmixopt => 0,
  cmix2nd => 1,
  cfcm2h  => 0.0075,
  cfcm2v  => 0.0075,
  cmix4th => 0,
  raydmp  => 0,
  flteps  => 0.05,

  moist   => 0,
  mphyopt => 0,

  "dirname" => qq('./'),
  thisdmp => 150.0,
  hdmpfmt => 9,
  mstout  => 0,
  tkeout  => 0,
  trbout  => 0,
  rainout => 0,
);

chdir $wrkdir;

@advopt = ( 1, 1, 2, 2, 2, 4, 2, 5 );

for ( $i = 0; $i < @advopt; $i+=2 ) {
  $newpar{madvopt} = $advopt[$i];
  $newpar{sadvopt} = $advopt[$i+1];

  $rname = "$runname$newpar{madvopt}$newpar{sadvopt}a";
  $newpar{runname} = qq('$rname');

  $runnum = ($i+2)/2;
  print "    $runnum) Run $rname\n";

  $newinput = "$wrkdir/$rname.input";
  mkinput( \%newpar, $arpsinput, $newinput );

  $time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
  print "       started at $time";

  $runcmd = "$srcdir/bin/$execmd < $rname.input > $rname.output";
  ( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";
  ( $norun ) or system( qq(grads -bl -c "exec density_grads.r $rname"));

  $time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
  print "       finished at $time\n";

}

exit 0
