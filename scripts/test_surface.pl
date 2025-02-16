#!/usr/bin/perl -w

use Cwd;
use File::Basename;
use File::Copy "cp";
use Getopt::Long;

use lib "./scripts";
use mkinput;

$nx_glob = 04;
$ny_glob = 04;
$nz_glob = 35;

#$nxcpu = 2;
#$nycpu = 2;

$nx = $nx_glob;
$ny = $ny_glob;
$nz = $nz_glob;

#$runname = "surface0a";
$subdir = "surface";
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
  print "  -version ver arps version\n";
  print "  -n           do not actually build and execute the program\n";
  print "  -s dir       ARPS source directory\n";
  print "  -w dir       working directory\n";
  print "\n";
  exit 1;
}

end_of_usage:

&GetOptions( "v", => \$verbose,
             "version=s", => \$arpsversion,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir );

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

$arpsinput = "$srcdir/input/arps.input";

cp( "$srcdir/scripts/surface_grads.r", "$wrkdir/surface_grads.r" );

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

   runmod  => 4,
   initopt => 1,
   pt0opt  => 0,

   dx => 32000.0,
   dy => 32000.0,
   dz => 500.0,

   tstop => 172800.0,
   dtbig => 30.0,
   dtsml => 30.0,

   wbc => 2,
   ebc => 2,
   nbc => 2,
   sbc => 2,
   tbc => 1,
   bbc => 1,

   vimplct => 1,

   tmixopt => 4,
   tkeopt  => 3,
   tmixcst => 0.0,
   cmix2nd => 1,
   cfcm2h  => 0.0,
   cfcm2v  => 0.0004,
   cmix4th => 0,
   raydmp  => 1,

   coriopt  => 1,
   coriotrm => 2,

   moist   => 1,
   mphyopt => 1,

   radopt  => 2,
   sfcphy  => 4,
   sfcdiag => 1,
   pbldopt => 2,

   "dirname" => qq('./'),
   thisdmp => 600.0,
   hdmpfmt => 9,
   varout  => 1,
   mstout  => 1,
   tkeout  => 1,
   trbout  => 1,
   rainout => 1,
   sfcout  => 1,
   landout => 1,
   prcout  => 1,
   radout  => 1,
   flxout  => 1,
);

print "\n";
print "    1) FIFE case\n";

$rname = "fife0a";
$newpar{runname} = qq('$rname');
$newpar{sndfile} = qq('$srcdir/sounding/fife.snd'); 
$newpar{initime} = qq('1987-07-11.12:00:00');
$newpar{ctrlat}  = 39.11;
$newpar{ctrlon}  = -96.56;
$newpar{styp}    = 7;
$newpar{vtyp}    = 3;
$newpar{lai0}    = 2.8;
$newpar{roufns0} = 0.065;
$newpar{veg0}    = 0.99;
$newpar{ptslnd0} = 295.15;
$newpar{tsoil0}  = 297.15;
$newpar{wetsfc0} = 0.27;
$newpar{wetdp0}  = 0.2556;

$newinput = "$wrkdir/$rname.input";
mkinput( \%newpar, $arpsinput, $newinput );

chdir $wrkdir;

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       started at $time";

$runcmd = "time $srcdir/bin/$execmd < $rname.input > $rname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";
( $norun ) or 
   system ( qq(grads -bl -c "exec surface_grads.r $rname $arpsversion" ));

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       finished at $time\n";

print "    2) Wangara case\n";

$rname  = "wangara0a";
$newpar{runname} = qq('$rname');
$newpar{sndfile} = qq('$srcdir/sounding/wangara.snd');   #'
$newpar{initime} = qq('1967-08-16.23:00:00');
$newpar{tmaxmin} = $newpar{tstop} * 2,
$newpar{ctrlat}  = -34.5;
$newpar{ctrlon}  = 144.93;
$newpar{styp}    = 5;
$newpar{vtyp}    = 1;
$newpar{lai0}    = 0.1;
$newpar{roufns0} = 0.24;
$newpar{veg0}    = 0.05;
$newpar{ptslnd0} = 273.683;
$newpar{tsoil0}  = 276.0;
$newpar{wetsfc0} = 0.1533;
$newpar{wetdp0}  = 0.1555;

$newinput = "$wrkdir/$rname.input";
mkinput( \%newpar, $arpsinput, $newinput );

chdir $wrkdir;

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       started at $time";

$runcmd = "time $srcdir/bin/$execmd < $rname.input > $rname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";
( $norun ) or 
   system ( qq(grads -bl -c "exec surface_grads.r $rname $arpsversion" ));

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       finished at $time\n";

exit 0
