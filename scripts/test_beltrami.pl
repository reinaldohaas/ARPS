#!/usr/bin/perl -w

use Cwd;
use File::Basename;
use File::Copy "cp";
use Getopt::Long;

use lib "./scripts";
use mkinput;

$nx_glob = 92;
$ny_glob = 62;
$nz_glob = 47;

#$nxcpu = 2;
#$nycpu = 2;

$nx = $nx_glob;
$ny = $ny_glob;
$nz = $nz_glob;

$runname = "beltrami0a";
$subdir = "beltrami";
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
$run = "";
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

cp( "$srcdir/scripts/beltrami_grads.r", "$wrkdir/beltrami_grads.r" );

$arpsinput = "$srcdir/input/arps.input";

# make the executable arps if necessary

unless(-x "$srcdir/bin/$execmd"){
  ( ! $verbose) or print "Making arps...\n";
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
  
# 3-d run mode
   runmod  => 1,

# initialization
   inibasopt => 6,
   initopt => 1,
   pt0opt  => 7,

# grid set up
   dx      => 3.0,
   dy      => 2.0,
   dz      => 1.0,

# time set up
   tstop   => 82.0,
   dtbig   => 0.2,
   dtsml   => 0.002,

# boundary conditions
   wbc     => 2,
   ebc     => 2,
   nbc     => 2,
   sbc     => 2,
   tbc     => 2,
   bbc     => 2,
   vimplct => 0,

# dynamics and physics configurations
   tmixopt => 1,
   tmixcst => 1.0,
   trbvimp => 1,
   cmix2nd => 0,
   cmix4th => 0,
   raydmp  => 0,
   moist   => 0,
   mphyopt => 0,

# output parameters
   "dirname" => qq('./'),
   thisdmp => 41.0,
   hdmpfmt => 9,
   tmaxmin => 1.0,
   mstout  => 0,
   tkeout  => 0,
   trbout  => 0,
   rainout => 0,
);

$newinput = "$wrkdir/$runname.input";

$newpar{runname} = qq('$runname');
mkinput( \%newpar, $arpsinput, $newinput );

chdir $wrkdir;

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       started  at $time";

$runcmd = "$srcdir/bin/$execmd < $runname.input > $runname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";
( $norun ) or system( qq(grads -lb -c "exec beltrami_grads.r $runname 1") );
( $norun ) or system( qq(grads -lb -c "exec beltrami_grads.r $runname 2") );

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       finished at $time";

exit 0
