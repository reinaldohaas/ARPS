#!/usr/bin/perl -w
#######################################################################
#
# USAGE:
#  
#   testarps.pl [options]
#
#   Options:
#     -h           help
#     -v           verbose
#     -n           do not actually build and execute the program
#     -a           archive in $arcdir
#     -s dir       ARPS source directory
#     -w dir       working directory
#     -d dir       data directory
#     -io bin/hdf  history dump options
#     -opt 0-4     compiler optimization option
#
#     -symmetry    run symmetry test
#     -surface     run surface tests
#     -density     run density current tests
#     -beltrami    run beltrami test
#     -warmbubble  run warm bubble test
#   
#     -arpstern    run arpstern for terrain
#     -arpstrn     run arpstrn for terrain
#     -arpssfc     run arpssfc for surface data
#     -ext2arps    run ext2arps to convert ETA data into ARPS grid
#     -adas        run adas for the initial data set
#     -realcase    run arps using 25 May 1998 SAMEX case
#     -arpsmpi     same as realcase but run with arpsmpi (4 PEs)
#
#     -all         run all tests including those to process data
#     -arps        run all tests for ARPS model
#     -mci2arps    run mci2arps tests         # not done
#     -nids2arps   run nids2arps tests
#     -agri        same as realcase but run with arpsagri  # not done
#     -data        run all tests for data processing
#
####################################################################

use Cwd;
use File::Copy "cp";
use File::Basename;
use Getopt::Long;

use lib "./scripts";
use mkinput;

$curdir = cwd();

###################################################################
#
# Set Envrionment variables first
#
#      setenv ARPSDIR /work/ywang
#      setenv ARPSVER 5.0.0Beta7
#
##################################################################

$arcdir = "$ENV{'ARPSDIR'}/std_test";
$arpsversion = "$ENV{'ARPSVER'}";
$topdir = "$ENV{'ARPSDIR'}";
$def_datadir = "$topdir/arps$arpsversion/data";

##################################################################
#
#  DONOT change below if you do not know PERL language
#
##################################################################

#($sec,$min,$hr,$mdy,$mon,$yr,$wdy,$ydy,$isdst) = localtime(time);
($dum,$dum,$dum,$dum,$mon,$yr,$dum,$dum,$dum) = localtime(time);
$month = (qw(jan feb mar apr may jun jul aug sep oct nov dec))[$mon];
$year = 1900 + $yr;
$year += 100 if ($yr < 80);   # Y2K patch
#$today = "$mdy$month$year";
$today = "$month$year";

$subdir = "test-arps$arpsversion-$today";


if ( $#ARGV < 0 ) { &Usage; }

&GetOptions( "h"      => \$help,
             "v",     => \$verbose,
             "n",     => \$norun,
             "a",     => \$archive,
             "s=s",   => \$srcdir,
             "w=s",   => \$wrkdir,
             "d=s",   => \$datadir,
             "io=s",  => \$io,
             "opt=s", => \$opt,
             "symmetry",   => \$symmetry,
             "surface",    => \$surface,
             "nids2arps",  => \$nids2arps,
             "mci2arps",   => \$mci2arps,
             "density",    => \$density,
             "beltrami",   => \$beltrami,
             "warmbubble", => \$warmbubble,
             "realcase",   => \$realcase,
             "arpssfc",    => \$arpssfc,
             "arpstern",   => \$arpstern,
             "arpstrn",    => \$arpstrn,
             "ext2arps",   => \$ext2arps,
             "adas",       => \$adas,
             "agri",       => \$agri,
             "arpsmpi",    => \$arpsmpi,
             "all",        => \$all,
             "arps!",      => \$arps,
             "data!",      => \$data );

if ( $help ) {  &Usage;     }

if ( $all ) {
  $beltrami   = 1;
  $density    = 1;
  $warmbubble = 1;
  $symmetry   = 1;
  $surface    = 1;
  $restart    = 1;
  $nids2arps  = 1;
#  $mci2arps   = 1;
  $arpssfc    = 1;
  $arpstern   = 1;
  $arpstrn    = 1;
  $ext2arps   = 1;
  $adas       = 1;
#  $agri       = 1;
  $realcase   = 1;
  $arpsmpi    = 1;
}

if ( $arps ) {
  $symmetry   = 1;
  $density    = 1;
  $beltrami   = 1;
  $warmbubble = 1;
  $surface    = 1;
#  $agri       = 1;
  $realcase   = 1;
  $arpsmpi    = 1;
}

if ( $data ) {
  $arpssfc    = 1;
  $arpstrn    = 1;
  $ext2arps   = 1;
  $adas       = 1;
#  $arpstern   = 1;
#  $nids2arps  = 1;
#  $mci2arps   = 1;
}

( $io  ) or $io  = "bin";
( $opt ) or $opt = "3";

$verb=''; $run='';
( ! $verbose ) or $verb = "-v";
( ! $norun   ) or $run  = "-n";

( $srcdir  ) or $srcdir  = "$topdir/arps${arpsversion}";
( $wrkdir  ) or $wrkdir  = "$srcdir/$subdir";
( $datadir ) or $datadir = $def_datadir;

print "wrkdir = $wrkdir\n";
( -d "$srcdir/src/arps" ) or die "ARPS root directory not found in $srcdir\n";
( -d $wrkdir ) or `mkdir -p $wrkdir`;

chdir $curdir;
chdir $srcdir;
$srcdir = cwd();

chdir $curdir;
chdir $wrkdir;
$wrkdir = cwd();
$subdir = basename( "$wrkdir", "" );

$scrptdir = "$srcdir/scripts";

print "\n";
print "Current directory: $curdir\n";
print "Source directory: $srcdir\n";
print "Working directory: $wrkdir\n";
print "Data directory: $datadir\n";
print "\n";

if ($data) {
  ( -d "$srcdir/data/arpstern.data" ) or 
   die "arpstern.data directory not found in $srcdir/data\n";

  ( -d "$srcdir/data/arpssfc.data" ) or 
   die "arpssfc.data directory not found in $srcdir/data\n";
}

####################################################################

if ( $arpstrn ) {
  print "\nTesting ARPSTRN...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpstrn.pl $verb $run -s $srcdir ".
            "-w $wrkdir/arpstrn -d $datadir -io $io -opt $opt";
  system( $syscmd ) == 0 or print "$0: test_arpstrn.pl failed\n";
}

if ( $arpstern ) {
  print "\nTesting ARPSTERN...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpstern.pl $verb $run ".
            "-s $srcdir -w $wrkdir/arpstern -d $datadir";
  system( $syscmd ) == 0 or print "$0: test_arpstern.pl failed\n";
}

if ( $ext2arps ) {
  print "\nTesting EXT2ARPS...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_ext2arps.pl $verb $run -io $io -opt $opt ".
            "-s $srcdir -w $wrkdir/ext2arps -d $datadir";
  system( $syscmd ) == 0 or print "$0: test_ext2arps.pl failed\n";
}

if ( $nids2arps ) {
  print "\nTesting nids2arps...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_nids2arps.pl $verb $run ".
            "-s $srcdir -w $wrkdir/nids2arps -d $datadir";
  system( $syscmd ) == 0 or print "$0: test_nids2arps.pl failed\n";
}

if ( $mci2arps ) {
  print "\nTesting mci2arps...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_mci2arps.pl $verb $run ".
            "-s $srcdir -w $wrkdir/mci2arps -d $datadir";
  system( $syscmd ) == 0 or print "$0: test_mci2arps.pl failed\n";
}

if ( $adas ) {
  print "\nTesting ADAS...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_adas.pl $verb $run -io $io -opt $opt ".
            "-s $srcdir -w $wrkdir/adas -d $datadir";
  system( $syscmd ) == 0 or print "$0: test_adas.pl failed\n";
}

if ( $arpssfc ) {
  print "\nTesting ARPSSFC...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpssfc.pl $verb $run -io $io -opt $opt ".
            "-s $srcdir -w $wrkdir/arpssfc -d $datadir";
  system( $syscmd ) == 0 or print "$0: test_arpssfc.pl failed\n";
}

if ( $arpsagr ) {
  print "\nTesting ARPSAGR...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpsagr.pl $verb $run -io $io -opt $opt ".
            "-s $srcdir -w $wrkdir/arpsqgr -d $datadir";
  system( $syscmd ) == 0 or print "$0: test_arpsagrc.pl failed\n";
}



if ( $symmetry ) {
  print "\nTesting ARPS with built-in symmetric setup...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_symmetry.pl $verb $run ".
            "-s $srcdir -w $wrkdir/symmetry";
  system( $syscmd ) == 0 or print "$0: test_symmetry.pl failed\n";
}

if ( $beltrami ) {
  print "\nTesting ARPS with Beltrami flow...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_beltrami.pl $verb $run ".
            "-s $srcdir -w $wrkdir/beltrami";
  system( $syscmd ) == 0 or print "$0: test_beltrami.pl failed\n";
}

if ( $density ) {
  print "\nTesting ARPS with density currents...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_density.pl $verb $run ".
            "-s $srcdir -w $wrkdir/density";
  system( $syscmd ) == 0 or print "$0: test_density.pl failed\n";
}

if ( $warmbubble ) {
  print "\nTesting ARPS with warm bubble using sounding may20.snd...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_warmbubble.pl $verb $run ".
            "-s $srcdir -w $wrkdir/warmbubble";
  system( $syscmd ) == 0 or print "$0: test_warmbubble.pl failed\n";
}

if ( $surface ) {
  print "\nTesting soil model in 1-D with FIFE and Wangara cases...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_surface.pl $verb $run -s $srcdir ".
            "-w $wrkdir/surface -version $arpsversion";
  system( $syscmd ) == 0 or print "$0: test_surface.pl failed\n";
}

if ( $restart ) {
  print "\nTesting restart using sounding may20.snd ...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_restart.pl $verb $run -s $srcdir ".
            "-w $wrkdir/surface -version $arpsversion";
  system( $syscmd ) == 0 or print "$0: test_restart.pl failed\n";
}

if ( $realcase ) {
  print "\nTesting ARPS using real case on 25 May 1998...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_realcase.pl $verb $run -io $io -opt $opt ".
            "-s $srcdir -w $wrkdir/realcase -d $datadir";
  system( $syscmd ) == 0 or print "$0: test_realcase.pl failed\n";
}

if ( $agri ) {
  print "\nTesting AGRI using real case on 25 May 1998...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_agri.pl $verb $run ".
            "-s $srcdir -w $wrkdir/agri -d $datadir";
  system( $syscmd ) == 0 or print "$0: test_agri.pl failed\n";
}

if ( $arpsmpi ) {
  print "\nTesting ARPSMPI using real case on 25 May 1998...\n";
  chdir $srcdir;
  $syscmd = "$scrptdir/test_arpsmpi.pl $verb $run -io $io -opt $opt ".
            "-s $srcdir -w $wrkdir/arpsmpi -d $datadir -real";
  system( $syscmd ) == 0 or print "$0: test_arpsmpi.pl failed\n";
}

chdir $wrkdir;
chdir "..";

( ! $archive ) or system( "cp -pr $subdir $arcdir" );

exit 0;

####################################################################
#
#  SUB Usage
#
####################################################################

sub Usage {
  print "\n";
  print "Usage: $0 [options]\n";
  print "\n";
  print "Options:\n";
  print "  -h           help\n";
  print "  -v           verbose\n";
  print "  -n           do not actually build and execute the program\n";
  print "  -a           archive in $arcdir\n";
  print "  -s dir       ARPS source directory\n";
  print "  -w dir       working directory\n";
  print "  -d dir       data directory\n";
  print "  -io bin/hdf  history dump option\n";
  print "  -opt 0-4     compiler optmization option\n";
  print "\nARPS simple tests\n";
  print "  -symmetry    run symmetry test\n";
  print "  -surface     run surface tests\n";
  print "  -density     run density current tests\n";
  print "  -beltrami    run beltrami test\n";
  print "  -warmbubble  run warm bubble test\n";
  print "\nA real case tests\n";
  print "  -arpstern    run arpstern for terrain\n";
  print "  -arpstrn     run arpstrn for terrain\n";
  print "  -arpssfc     run arpssfc for surface data\n";
  print "  -ext2arps    run ext2arps to convert ETA data into ARPS grid\n";
  print "  -adas        run adas for the initial data set\n";
  print "  -realcase    run arps using 25 May 1998 SAMEX case\n";
  print "  -arpsmpi     same as realcase but run with arpsmpi (4 PEs)\n";
  print "\nMisc. tests\n";
  print "  -all         run all tests including those to process data\n";
  print "  -arps        run all tests for ARPS model\n";
  print "  -data        run all tests for data processing\n";
  print "  -agri        same as realcase but run with arpsagri\n";
  print "  -nids2arps   run nids2arps tests\n";
  print "  -mci2arps    run mci2arps tests\n";
  print "\n";
  exit 1;
}
