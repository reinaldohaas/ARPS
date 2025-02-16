#!/usr/bin/perl -w

use Cwd;
use File::Basename;
use File::Copy "cp";
use Getopt::Long;

use lib "./scripts";
use mkinput;

$nx = 67;
$ny = 67;
$nz = 43;

$runname = "sat25may1998";
$execmd = "mci2arps";

$curdir = cwd();
$topdir = "$ENV{'ARPSDIR'}";
$arpsversion = "$ENV{'ARPSVER'}";
$subdir = "mci2arps";

goto end_of_usage;

usage: {
  print "\n";
  print "Usage: $0 [options]\n";
  print "\n";
  print "Options:\n";
  print "  -v           verbose\n";
  print "  -n           do not actually build and execute the program\n";
  print "  -w dir       working directory\n";
  print "  -d dir       data directory\n";
  print "  -s dir       ARPS source directory\n";
#  print "  -arpstern    option to run arpstern for terrain data\n";
  print "\n";
  exit 1;
}

end_of_usage:

&GetOptions( "v", => \$verbose,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir,
             "d=s", => \$datadir
#             "arpstern", => \$arpstern 
              );

$verb = "";
$run  = "";

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

#$scrptdir = "$srcdir/scripts";

$satdir  = "$datadir/data.25may1998/sat";

#$tdatdir = "$datadir/data.test";

if( $verbose ) {
  print "\n";
  print "Current directory: $curdir\n";
  print "Source directory: $srcdir\n";
  print "Working directory: $wrkdir\n";
  print "Data directory: $datadir\n";
  print "\n";
}

$arpsinput = "$srcdir/input/arps.input";

$satinput = "$wrkdir/$runname.input";

# set dimensional parameters for ARPS

#$trndata = "trn32km${nx}x${ny}.trndata";
#$sfcdata = "sfc32km${nx}x${ny}.sfcdata";

#$soilvar = "$runname.soilvar.0000";

unless (-s "palgrey.hdf") { print `cp $satdir/mci2arps/palgrey.hdf .`;}


# set dimensional parameters for ARPS
#%dims = ( nx => $nx,
#          ny => $ny,
#          nz => $nz,
#          nxebc => 01,
#          nyebc => 01,
#          nzebc => 01,
#          nx_radiat => 01,
#          ny_radiat => 01,
#          nz_radiat => 01,
#);

# set the input parameters for mci2arps
%newpar = (
  nx => $nx,
  ny => $ny,
  nz => $nz,

  runmod  => 1,
  runname => qq('$runname'),
  initime => qq('1998-05-25.00:00:00'),
  initopt => 3,
  inifmt  => 10,

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

  hdmpfmt => 10,

  dirname => qq('./')

);

$newpar{tinitebd} = $newpar{initime};
print "Making the input namelist file...\n" if ( $verbose );
mkinput( \%newpar, $arpsinput, $satinput );

chdir $srcdir;

$makecmd = "makearps $run -opt 3 $execmd";
system($makecmd) == 0 or die "$0: $makecmd failed";

chdir $wrkdir;

my @satname = qw(goes08ch4_19980528_0015 goes08vis_19980528_0015);

foreach $sat ( @satname ) {

  $runcmd = "$srcdir/bin/$execmd $satdir/$sat < $runname.input > $runname.output.$sat";
  print "running command: $runcmd\n";

  ( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";


#  add copy to data dir

}

print `cp $runname.*.cttemp $runname.*.albedo $satdir/mci2arps/.`;

foreach $hdf ( glob("$runname.*.hdf")) {

  my $gif = $hdf;
  $gif =~ s/hdf/gif/;
  $runcmd = "imscale -infile $hdf -outfile $gif";
  print "running command: $runcmd\n";
  ($norun) or system ($runcmd) == 0 or warn "$0: $runcmd failed";

}

exit 0;
