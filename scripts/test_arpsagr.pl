#!/usr/bin/perl -w

use Cwd;
use File::Copy "cp";
use File::Basename;
use Getopt::Long;

use lib "./scripts";
use mkinput;
use updinc;
use getincpar;

$nx_glob = 67;
$ny_glob = 67;
$nz_glob = 43;

$nxcpu = 2;
$nycpu = 2;

$nx = $nx_glob;
$ny = $ny_glob;
$nz = $nz_glob;

$runname = "arpsagr25may1998";
$adasname = "adas25may1998";
$exbcname = "eta25may1998";
$execmd = "arpsagr";

$topdir = "/home/haas";
$curdir = cwd();
$subdir = "arpsagr";

$arpsversion = "5.3.4";

goto end_of_usage;

usage: {
  print "\n";
  print "Usage: $0 [options]\n";
  print "\n";
  print "Options:\n";
  print "  -v           verbose\n";
  print "  -n           do not actually build and execute the program\n";
  print "  -w dir	working directory\n";
  print "  -s dir	ARPS source directory\n";
  print "  -arpssfc	option to run arpssfc for surface data\n";
  print "  -arpstern	option to run arpstern for terrain\n";
  print "  -ext2arps	option to run ext2arps to convert ETA data into ARPS grid\n";
  print "  -adas	option to run adas for the initial data set\n";
  print "  -data	option to run additional tests for data processing\n";
  print "\n";
  exit 1;
}

end_of_usage:

&GetOptions( "v", => \$verbose,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir,
             "d=s", => \$datadir,
             "arpssfc", => \$arpssfc,
             "arpstern", => \$arpstern,
             "ext2arps", => \$ext2arps,
             "adas", => \$adas,
             "data", => \$data );

if ( $data ) {
  $arpssfc  = 1;
  $arpstern = 1;
  $ext2arps = 1;
  $adas     = 1;
}

$verb='';
$run ='';
( ! $verbose ) or $verb = "-v";
( ! $norun   ) or $run  = "-n";

( $srcdir  ) or $srcdir  = "$topdir/arps${arpsversion}";
( $wrkdir  ) or $wrkdir  = "$srcdir/$subdir";
( $datadir ) or $datadir = "/home/haas/data";

( -d "$srcdir/src/arps" ) or die "ARPS root directory not found in $srcdir\n";
( -d $wrkdir ) or mkdir $wrkdir, 0755;

chdir $curdir;
chdir $srcdir;
$srcdir = cwd();

chdir $curdir;
chdir $wrkdir;
$wrkdir = cwd();

$scrptdir = "$srcdir/scripts";

$sfcdir  = "$datadir/data.arps/arpssfc";
$trndir  = "$datadir/data.arps/arpstern";
$etadir  = "$datadir/data.25may1998/eta40grb";
$tdatdir = "$datadir/data.test";

$arpsinput = "$srcdir/input/arps.input";
$dimsinc = "$srcdir/include/dims.inc";

# set dimensional parameters for ARPS
%dims = ( nx => $nx,
          ny => $ny,
          nz => $nz,
          nxebc => "nx",
          nyebc => "ny",
          nzebc => "nz",
          nx_radiat => "nx",
          ny_radiat => "ny",
          nz_radiat => "nz",
);

print "Updating dims.inc...\n" if ( $verbose );
updinc( \%dims, $dimsinc );

$trndata = "trn32km${nx}x${ny}.trndata";
$sfcdata = "sfc32km${nx}x${ny}.sfcdata";
$inifile  = "$adasname.grb000000";
$inigbf   = "$adasname.grbgrdbas";
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

    print"$tdatdir/$sfcdata", "$wrkdir/$sfcdata\n";
if ( $arpssfc ) {
  print "\nCalling test_arpssfc.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpssfc.pl $verb $run -s $srcdir -w $wrkdir -d $datadir";
  system( $syscmd ) ==0 or die "$0: test_arpssfc.pl failed";
}
else {
  ( -r "$wrkdir/$sfcdata" )
    or cp ( "$tdatdir/$sfcdata", "$wrkdir/$sfcdata" ) == 1
    or die "cp() failed";
}

if ( $arpstern ) {
  print "\nCalling test_arpstern.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpstern.pl $verb $run -s $srcdir -w $wrkdir -d $datadir";
  system( $syscmd ) ==0 or die "$0: test_arpstern.pl failed";
}
else {
  ( -r "$wrkdir/$trndata" )
    or cp ( "$tdatdir/$trndata", "$wrkdir/$trndata" ) == 1
    or die "cp() failed";
}

if ( $ext2arps ) {
  print "\nCalling test_ext2arps.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_ext2arps.pl $verb $run -s $srcdir -w $wrkdir -d $datadir";
  system( $syscmd ) ==0 or die "$0: test_ext2arps.pl failed";
}
else {
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

if ( $adas ) {
  print "\nCalling test_adas.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_adas.pl $verb $run -s $srcdir -w $wrkdir -d $datadir";
  system( $syscmd ) ==0 or die "$0: test_adas.pl failed";
}
else {
  $inifile = "$exbcname.grb000000" if ( ! -r "$tdatdir/$inifile" );
  $inigbf  = "$exbcname.grbgrdbas" if ( ! -r "$tdatdir/$inigbf" );
  ( -r "inifile" )
    or cp ( "$tdatdir/$inifile","$wrkdir/$inifile" ) == 1
    or die "$0: cp() failed";
  ( -r "inigbf" )
    or cp ( "$tdatdir/$inigbf", "$wrkdir/$inigbf" ) == 1
    or die "$0: cp() failed";
}

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

# make the executable arps

print "Making arpsagr...\n" if ( $verbose );
chdir $srcdir;

$makecmd = "./makearps -m linux $run -opt 3 $execmd";
system($makecmd) == 0 or die "$0: $makecmd failed";

# set the input parameters for ext2arps

%newpar = (
  nxc           => $nx,
  nyc           => $ny,
  nzc           => $nz,
  grdsrt        => qq(.true.),
  intrat        => 3,
  intratt       => 3,
  nfinelv       => 1,
  "ngrdnew(1)"  => 1,
  "ixc(1,1)"    => 34.0,
  "jyc(1,1)"    => 32.5,
  "ixln(1,1)"   => 22.0,
  "jyln(1,1)"   => 21.0,
  "gangle(1,1)" => 0.0,

  runmod  => 1,

  initime => qq('1998-05-25.00:00'),
  initopt => 3,
  inifmt  => 10,
  inifile => qq('$inifile'),
  inigbf  => qq('$inigbf'),

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

  hdmpfmt => 10,
  thisdmp => 3600.0,
  trstout => 0.0,

  dirname => qq('./'),
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

$newpar{runname} = qq("$runname");
$newpar{tinitebd} = $newpar{initime};

print "Creating input file...\n" if ( $verbose );
$newinput = "$wrkdir/$runname.input";
mkinput( \%newpar, $arpsinput, $newinput );

# estimate the total space for all grids

print "Getting the numbers of arrays from agrigrid.inc...\n" if ( $verbose );
$gridh  = "$srcdir/include/agrigrid.inc";
@vkeys = ( "nsint", "nsreal", "nslogc", "ns1d",  "nx1d",   "ny1d",
           "nz1d",  "nxy2d",  "nxz2d",  "nyz2d", "nxyz3d", "nexbc3d" );

%gridhpar = getincpar( \@vkeys, "$gridh" );

print "Getting the coarse grid information...\n" if ( $verbose );
$nxc = $newpar{nxc};
$nyc = $newpar{nyc};
$nzc = $newpar{nzc};

$nxyc  = $nxc * $nyc;
$nxzc  = $nxc * $nzc;
$nyzc  = $nyc * $nzc;
$nxyzc = $nxyc * $nzc;

print "Getting the fine grid information...\n" if ( $verbose );
$nx1 = $newpar{q(ixln(1,1))} * $newpar{intrat} + 1;
$ny1 = $newpar{q(jyln(1,1))} * $newpar{intrat} + 1;
$nz1 = $nzc;
$nxy1  = $nx1 * $ny1;
$nxz1  = $nx1 * $nzc;
$nyz1  = $ny1 * $nzc;
$nxyz1 = $nxy1 * $nzc;

$nxmax = $nxc > $nx1 ? $nxc : $nx1;
$nymax = $nyc > $ny1 ? $nyc : $ny1;
$nzmax = $nzc > $nz1 ? $nzc : $nz1;
$nxymax = $nxmax * $nymax;
$nxyzmax = $nxymax * $nzmax;
$nxyz0max = ($nxmax+1) * ($nymax+1) * ($nzmax+1);

if ( $verbose ) {
  print "nxc = $nxc\n";
  print "nyc = $nyc\n";
  print "nzc = $nzc\n";
  print "nxyc = $nxyc\n";
  print "nxzc = $nxzc\n";
  print "nyzc = $nyzc\n";
  print "nxyzc = $nxyzc\n";

  print "nx1 = $nx1\n";
  print "ny1 = $ny1\n";
  print "nz1 = $nz1\n";
  print "nxy1 = $nxy1\n";
  print "nxz1 = $nxz1\n";
  print "nyz1 = $nyz1\n";
  print "nxyz1 = $nxyz1\n";

  print "nsint = $gridhpar{nsint}\n";
  print "nsreal = $gridhpar{nsreal}\n";
  print "nslogc = $gridhpar{nslogc}\n";
  print "ns1d = $gridhpar{ns1d}\n";
  print "nx1d = $gridhpar{nx1d}\n";
  print "ny1d = $gridhpar{ny1d}\n";
  print "nz1d = $gridhpar{nz1d}\n";
  print "nxy2d = $gridhpar{nxy2d}\n";
  print "nxz2d = $gridhpar{nxz2d}\n";
  print "nyz2d = $gridhpar{nyz2d}\n";
  print "nxyz3d = $gridhpar{nxyz3d}\n";
}

print "Estimating the size of storage set in agrialloc.inc...\n" if ( $verbose );

$allocsize = 0;

print "Adding 4 3-d arrays for AGRI temporary space\n" if ( $verbose );
$allocsize = $allocsize + $nxyzmax * 4;

print "alloc_permanent_array = $allocsize\n" if ( $verbose );

print "Adding permanent arrays for all ARPS grids\n" if ( $verbose );
$allocsize = $allocsize
           + $gridhpar{nsint} * 2
           + $gridhpar{nsreal} * 2
           + $gridhpar{nslogc} * 2
           + $gridhpar{ns1d} * 2
           + $gridhpar{nx1d} * ($nxc+$nx1)
           + $gridhpar{ny1d} * ($nyc+$ny1)
           + $gridhpar{nz1d} * ($nzc+$nz1)
           + $gridhpar{nxy2d} * ($nxyc+$nxy1)
           + $gridhpar{nxz2d} * ($nxzc+$nxz1)
           + $gridhpar{nyz2d} * ($nyzc+$nyz1)
           + $gridhpar{nxyz3d} * ($nxyzc+$nxyz1) ;

print "alloc_permanent_array = $allocsize\n" if ( $verbose );

print "Adding permanent space for base grid EXBC\n" if ( $verbose );
if ( $newpar{lbcopt} == 2 ) { $allocsize = $allocsize + $nxyzc * 22 };

print "alloc_exbc = $allocsize\n" if ( $verbose );

print "Adding temporary space for radiation\n" if ( $verbose );
if ( $newpar{radopt} == 2 ) {
  $radcstinc  = "$srcdir/include/radcst.inc";
  @vkeys = ( "n2d_radiat", "n3d_radiat" );
  %radcstpar = getincpar( \@vkeys, "$radcstinc" );

  $allocsize = $allocsize
             + $radcstpar{n2d_radiat} * $nxymax
             + $radcstpar{n3d_radiat} * $nxyzmax ;
}

print "alloc_radiation = $allocsize\n" if ( $verbose );

print "Adding temporary space for ARPS grids\n" if ( $verbose );
$allocsize = $allocsize + $nxyzmax * 26;
$allocsize = $allocsize + $nxyz0max * 3;
$allocsize = $allocsize + $nxymax * 10;

print "alloc_temporary = $allocsize\n" if ( $verbose );

$allocsize = ( int $allocsize/1000000 + 1 ) * 1000000;

print "alloc_final = $allocsize\n" if ( $verbose );

print "Updating agrialloc.inc with (lstore=$allocsize)...\n" if ( $verbose );
$allocinc = "$srcdir/include/agrialloc.inc";
$allocpar{lstore} = $allocsize;
updinc( \%allocpar, "$allocinc" );

chdir $wrkdir;

chomp( $time = `date +%H:%M' '%a' '%h' '%d' '%Y` );
print "\n";
print "       Run $runname started at $time\n";

$runcmd = "time $srcdir/bin/$execmd < $runname.input > $runname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

chomp( $time = `date +%H:%M' '%a' '%h' '%d' '%Y` );
print "       Run $runname finished at $time\n";

exit 0
