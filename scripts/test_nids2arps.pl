#!/usr/bin/perl -w
#####################################################################################
#
#    Usage: test_nids2arps.pl [options]
#     
#       Options:
#         -v           verbose
#         -n           do not actually build and execute the program
#         -w dir       working directory
#         -d dir       data directory
#         -s dir       ARPS source directory
#         -io bin/hdf  history dump format
#         -opt 0-4     compiler optimization option
#         -arpstern    option to run arpstern for terrain data
#         -arpstrn     option to run arpstrn for terrain data
#       
#
#####################################################################################

use Cwd;
use File::Basename;
use File::Copy "cp";
use Getopt::Long;

use lib "./scripts";
use mkinput;

$nx = 67;
$ny = 67;
$nz = 43;

$runname = "nids25may1998";
$execmd = "nids2arps";

$curdir = cwd();
$topdir = "$ENV{'ARPSDIR'}";
$arpsversion = "$ENV{'ARPSVER'}";
$subdir = "nids2arps";

@radname = ( AMA, ARX, DDC, DFX, DYX, EWX, FDR, FDX, FWS,
             GRK, HGX, ICT, INX, LBB, LOT, LSX, LZK,
             MAF, PUX, SGF, SHV, SJT, SRX, TLX, VNX );

&GetOptions( "v", => \$verbose,
             "h", => \$help,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir,
             "d=s", => \$datadir,
             "io=s",  => \$io,
             "opt=s", => \$opt,
             "arpstern", => \$arpstern,
             "arpstrn",  => \$arpstrn );
&Usage if ($help);

( $io ) or $io = "bin";
( $opt) or $opt= 3;

$verb = ""; $run  ="";
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

$scrptdir = "$srcdir/scripts";

$nidsdir  = "$datadir/data.25may1998/nids";
$etadir  = "$datadir/data.25may1998/eta40grb";

$tdatdir = "$datadir/data.test";

if( $verbose ) {
  print "\n";
  print "Current directory: $curdir\n";
  print "Source directory: $srcdir\n";
  print "Working directory: $wrkdir\n";
  print "Data directory: $datadir\n";
  print "\n";
}

$arpsinput = "$srcdir/input/nids2arps.input";
$nidsinput = "$wrkdir/$runname.input";

system("cat $srcdir/input/arps.input $srcdir/input/nids2arps.input_incomplete".
       " > $srcdir/input/nids2arps.input") unless (-r "$srcdir/input/nids2arps.input");

$trndata = "trn32km${nx}x${ny}.trndata";

if ( $arpstrn ) {
  print "\nCalling test_arpstrn.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpstrn.pl $verb $run -io $io -opt $opt -s $srcdir -w $wrkdir -d $datadir";
  system( $syscmd ) == 0 or die "$0: test_arpstrn.pl failed";
}elsif ( $arpstern) {
  print "\nCalling test_arpstern.pl...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpstern.pl $verb $run -s $srcdir -w $wrkdir -d $datadir";
  system( $syscmd ) == 0 or die "$0: test_arpstern.pl failed";
} elsif ( ! -r $trndata ) {
  (-r "$wrkdir/$trndata" )
    or cp ( "$tdatdir/$trndata", "$wrkdir/$trndata" ) == 1
    or die "cp() failed";
}

chdir $wrkdir;
( -r $trndata ) or die "$0: $trndata not found";

#----------------------------------------------------------------------------
#
# set the input parameters for nids2arps
#
#----------------------------------------------------------------------------

%newpar = (
  nx => $nx,
  ny => $ny,
  nz => $nz,
  
  runmod  => 1,
  runname => qq('$runname'),
  initime => qq('1998-05-25.00:00:00'),
  initopt => 3,
  inifmt  => 10,
  inifile => qq('$runname.grb000000'),
  inigbf  => qq('$runname.grbgrdbas'),

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

  ternopt => 2,
  terndta => qq('$trndata'),

  hdmpfmt => 10,
  thisdmp => 3600.0,
  trstout => 3600.0,

  "dirname" => qq('./'),
  exbcdmp  => 1,
  varout  => 1,
  mstout  => 1,
  iceout  => 0,
  tkeout  => 0,
  trbout  => 0,
  rainout => 0,
  sfcout  => 1,
  landout => 0,
  prcout  => 0,
  radout  => 0,
  flxout  => 0,
  qcexout => 0,
  qrexout => 0,
  qiexout => 0,
  qsexout => 0,
  qhexout => 0,

  dir_extd    => qq('$etadir'),
  extdopt     => 2,
  extdname    => qq('/eta40grb'),
  nextdfil    => 3,
  "extdtime(1)" => qq('1998-05-25.00:00:00+000:00:00'),
  "extdtime(2)" => qq('1998-05-25.00:00:00+003:00:00'),
  "extdtime(3)" => qq('1998-05-25.00:00:00+006:00:00'),

# Added nids_data Namelist  on July 26 2001
  radar_name => qq('KTLX'),
  reflistfn  => qq('reflist'),
  vellistfn  => qq('vellist'),
  map_flg    => 0,
  radar_map_file => qq('xxxxx'),
  dir_name   => qq('./'),

  etfn => qq('/work/ekemp/decodenids/TLX199905040009.ETP'),
  et_remapopt => 2,
  et_radius => 4000.,
  vilfn => qq('/work/ekemp/decodenids/TLX199905040009.VIL'),
  vil_remapopt => 2,
  vil_radius => 4000.,
  dpafn => qq('/work/ekemp/decodenids/TLX199905040012.DPA'),
  dpa_remapopt => 2,
  dpa_radius => 4000.,
);

$newpar{ternfmt} = ($io eq "hdf")?1:0;
$newpar{tinitebd} = $newpar{initime};
print "Making the input namelist file...\n" if ( $verbose );
mkinput( \%newpar, $arpsinput, $nidsinput );

#----------------------------------------------------------------
#
# Make nids2arps
#
#----------------------------------------------------------------

chdir $srcdir;

$makecmd = "makearps $run -opt $opt $execmd";
system($makecmd) == 0 or die "$0: $makecmd failed";

#---------------------------------------------------------------
#
# Run nids2arps for each radar
#
#---------------------------------------------------------------

chdir $wrkdir;

foreach $rad ( @radname ) {

  open RFLIST, "$nidsdir/$rad/$rad\_rf_list" or 
     die "ERROR opening nids list file $nidsdir/$rad/$rad\_rf_list ($!)";
  @rflist = <RFLIST>;
  close RFLIST;
  open VLLIST, "$nidsdir/$rad/$rad\_vl_list" or 
     die "ERROR opening nids list file $nidsdir/$rad/$rad\_vl_list ($!)";
  @vllist = <VLLIST>;
  close VLLIST;

  $rflist = "$rad\_rf_list";
  $vllist = "$rad\_vl_list";
  open RFLIST, ">$rflist" or die "ERROR creating $rflist ($!)";
  open VLLIST, ">$vllist" or die "ERROR creating $vllist ($!)";
  foreach (@rflist) {
    $_ =~ s/.*25may1998\/data(.*)/$datadir\/data.25may1998$1/;
    print RFLIST "$_";
  }
  foreach (@vllist) {
    $_ =~ s/.*25may1998\/data(.*)/$datadir\/data.25may1998$1/;
    print VLLIST;
  }
  close RFLIST;
  close VLLIST;

  $newpar{"radar_name"} = qq('K$rad');
#  $newpar{"reflistfn"}  = qq('$rflist');
#  $newpar{"vellistfn"}  = qq('$vllist');
#
#  @etpdir=`ls $nidsdir/$rad/ETP`;
#  chop $etpdir[0];
#  if( ! $etpdir[0]) { $etpdir[0]='NONE';  }
#  $newpar{"etfn"}  = qq('$nidsdir/$rad/ETP/$etpdir[0]');
#
#  @vildir=`ls $nidsdir/$rad/VIL`;
#  chop $vildir[0];
#  if( ! $vildir[0]) { $vildir[0]='NONE';  }
#  $newpar{"vilfn"}  = qq('$nidsdir/$rad/VIL/$vildir[0]');
#
#  @dpadir=`ls $nidsdir/$rad/DPA`;
#  chop $dpadir[0];
#  if( ! $dpadir[0]) { $dpadir[0]='NONE';  }
#  $newpar{"dpafn"}  = qq('$nidsdir/$rad/DPA/$dpadir[0]');

  mkinput( \%newpar, $arpsinput, $nidsinput );

  if ($io eq "hdf") {
    $runcmd = "$srcdir/bin/$execmd K$rad $rflist $vllist -hdf 5  < $runname.input > $runname.output.K$rad";
  } else {
    $runcmd = "$srcdir/bin/$execmd K$rad $rflist $vllist -binary < $runname.input > $runname.output.K$rad";
  }
  print "Running for $rad: \n$runcmd\n";

  ( $norun ) or ( system($runcmd) == 0 ) or die "$0: $runcmd failed";

  print "$rad data processed successfully\n";
  
#  system( "mv 'nids25may1998.NIDSdp000000' 'K$rad.NIDSdp0000' ");
#  system( "mv 'nids25may1998.NIDSet000000' 'K$rad.NIDSet0000' ");
#  system( "mv 'nids25may1998.NIDSvl000000' 'K$rad.NIDSvl0000' ");

  @n2afiles = <$wrkdir/K$rad.*>;
  foreach $n2afile ( @n2afiles ) {
    $nfile = basename( "$n2afile", "" );
    cp( "$wrkdir/$nfile", "$tdatdir/nids/nids2arps/$nfile" ) == 1
      or die "cp() failed";
  }
}

exit 0;

################################################################################
#
# Sub Usage
#
################################################################################

sub Usage {
  print "\n";
  print "Usage: $0 [options]\n";
  print "\n";
  print "Options:\n";
  print "  -h           help\n";
  print "  -v           verbose\n";
  print "  -n           do not actually build and execute the program\n";
  print "  -w dir       working directory\n";
  print "  -d dir       data directory\n";
  print "  -s dir       ARPS source directory\n";
  print "  -io bin/hdf  history dump format\n";
  print "  -opt 0-4     compiler optimization level\n";
  print "  -arpstern    option to run arpstern for terrain data\n";
  print "  -arpstrn     option to run arpstrn for terrain data\n";
  print "\n";
  exit 1;
}
