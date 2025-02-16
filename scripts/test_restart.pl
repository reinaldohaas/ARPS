#!/usr/bin/perl -w
# ==================================================================
# Two variables need to check or change in this script
#
#  $topdir, $arpsversion
# ==================================================================

use Cwd;
use File::Basename;
use File::Copy "cp";
use Getopt::Long;

use lib "./scripts";
use mkinput;

$nx_glob = 67;
$ny_glob = 67;
$nz_glob = 43;

$nxcpu = 2;
$nycpu = 2;

$nx = $nx_glob;
$ny = $ny_glob;
$nz = $nz_glob;

$runname = "restart";
$subdir = "restart";
$execmd = "arps";
$execmd2 = "arps_mpi";

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
$verb = "";

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

$arpsinput = "$srcdir/input/arps.input";

# make the executable arps if necessary

unless(-x "$srcdir/bin/$execmd"){
  print "Making arps...\n" if ( $verbose );
  chdir $srcdir;

  $makecmd = "makearps $run -opt 3 $execmd";
  system($makecmd) == 0 or die "$0: $makecmd failed";
}

# make the executable arps_mpi if necessary

unless(-x "$srcdir/bin/$execmd2"){
  print "Making arps_mpi...\n" if ( $verbose );
  chdir $srcdir;

  $makecmd = "makearps $run -opt 3 $execmd2";
  system($makecmd) == 0 or die "$0: $makecmd failed";
}

chdir $srcdir;

$makecmd = "makearps $run -opt 3 splitfiles";
system( $makecmd ) ==0 or die "$0: $makecmd failed";

$makecmd = "makearps $run -opt 3 joinfiles";
system( $makecmd ) ==0 or die "$0: $makecmd failed";

# set the input parameters for ext2arps

%newpar = (
# set dimensional parameters for ARPS
  nx => $nx,
  ny => $ny,
  nz => $nz,

  nproc_x => $nxcpu,
  nproc_y => $nycpu,
  max_fopen => $nxcpu*$nycpu,   

# 3-d run mode
  runmod  => 1,

# initialization
  initopt   => 1,
  pt0opt    => 1,
  ptpert0   => 4.0,
  sndfile   => qq('$srcdir/sounding/may20.snd'),

# grid set up
  dx      => 1000.0,
  dy      => 1000.0,
  dz      => 500.0,

# time set up
  tstop   => 7200.0,
  dtbig   => 6.0,
  dtsml   => 1.0,

# boundary conditions
  wbc     => 4,
  ebc     => 4,
  nbc     => 4,
  sbc     => 4,
  tbc     => 1,
  bbc     => 1,
  vimplct => 1,

# dynamics and physics configurations
  tmixopt => 4,
  tmixcst => 1.0,

  cmix2nd => 1,
  cfcm2h  => 0.0,
  cfcm2v  => 4.0e-4,
  cmix4th => 1,
  cfcm4h  => 1.0e-3,
  cfcm4v  => 0.0,

  raydmp  => 1,

  moist   => 1,
  mphyopt => 1,

# output parameters
  "dirname" => qq('./'),
  thisdmp => 600.0,
  hdmpfmt => 9,
  tmaxmin => 60.0,
  varout  => 1,
  mstout  => 1,
  iceout  => 0,
  tkeout  => 1,
  trbout  => 1,
  rainout => 1,
  trstout => 3600.0,
);

####################################################################
#
#  1. Run normal model to 7200 seconds
#
####################################################################

chdir $wrkdir;

$rname = "${runname}";
$newpar{runname} = qq('$rname');

$newinput = "$wrkdir/$rname.input";
mkinput( \%newpar, $arpsinput, $newinput );

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $rname started at $time";

$runcmd = "time $srcdir/bin/$execmd < $rname.input > $rname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $rname finished at $time\n";


####################################################################
#
#  2. Run normal restart from 3600 seconds to 7200 seconds
#
####################################################################

chdir $wrkdir;

$newpar{initopt} = 2;
$newpar{rstinf}  = qq('$rname.rst003600');   # use $rname from step 1

$rname = "${runname}_rst";
$newinput = "$wrkdir/$rname.input";
$newpar{runname} = qq('$rname');
mkinput( \%newpar, $arpsinput, $newinput );

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Rerun $rname started at $time";

$runcmd = "time $srcdir/bin/$execmd < $rname.input > $rname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Rerun $rname finished at $time\n";

####################################################################
#
#  3. Run MPI restart from 3600 seconds to 7200 seconds
#
####################################################################

chdir $wrkdir;

$newpar{initopt} = 2;         
$rname = "${runname}_mpirst";        # step 1 and 2 must run properly

$newinput = "$wrkdir/$rname.input";
$newpar{runname} = qq('$rname');
mkinput( \%newpar, $arpsinput, $newinput );

print "Splitting files...\n";
$syscmd = "$srcdir/bin/splitfiles < $rname.input > $rname.split.output";
system( $syscmd ) ==0 or die "$0: $syscmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Rerun $rname started at $time";

$runcmd = "time mpirun -np 4 $srcdir/bin/$execmd2 < $rname.input > $runname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Rerun $rname finished at $time\n";

print "Joining files...\n";
$syscmd = "$srcdir/bin/joinfiles < $rname.input > $rname.join.output";
$joinstat = system( $syscmd );
if ($joinstat != 0) {
    print "$0: WARNING $syscmd failed status $joinstat ($?)\n";
}
  
exit 0
