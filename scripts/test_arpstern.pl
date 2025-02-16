#!/usr/bin/perl -w

use Cwd;
use File::Copy 'cp';
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

$subdir = "arpstern";
$execmd = "arpstern";
$rname = "arpstern";
$run="";

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
  print "  -s dir       ARPS source directory\n";
  print "  -w dir       working directory\n";
  print "  -d dir       data directory\n";
  print "\n";
  exit 1;
}

end_of_usage:

&GetOptions( "v", => \$verbose,
             "n", => \$norun,
             "w=s", => \$wrkdir,
             "s=s", => \$srcdir,
             "d=s", => \$datadir );

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

#$scrptdir = "$srcdir/scripts";

#$etadir  = "$datadir/data.25may1998/eta40grb";
#$sfcdir  = "$datadir/arpssfc.data";
$trndir  = "$datadir/arpstern.data";

$tdatdir = "$datadir/data.test";

if ( $verbose ) {
  print "\n";
  print "Current directory: $curdir\n";
  print "Source directory: $srcdir\n";
  print "Working directory: $wrkdir\n";
  print "Data directory: $datadir\n";
  print "\n";
}

$arpsterninput = "$srcdir/input/arpstern.input";
$trninput = "$wrkdir/$rname.input";

$trndata = "trn32km${nx}x${ny}.trndata";
#$sfcdata = "sfc32km${nx}x${ny}.sfcdata";

# set the input parameters for ext2arps

%newpar = (
# set dimensional parameters for ARPS
  nx => $nx,
  ny => $ny,
  nz => $nz,
   
  itertype => 2,
  comtype  => 1,
  tdatadir => qq('$tdatdir'),
  terndir  => qq('$trndir'),

  mapproj => 2,
  trulat1 => 30.0,
  trulat2 => 60.0,
  trulon  => -100.0,
  sclfct  => 1.0,

  dx      => 32000.0,
  dy      => 32000.0,
  ctrlat  =>    38.0,
  ctrlon  =>   -98.0,
);

print "Making the input namelist file...\n";
mkinput( \%newpar, $arpsterninput, $trninput );

#----------- ADD BY YUNHENG WANG -------------------------------
# run dir1deg, dir5min, dir30sec separately if necessary

%dirdat=("dir1deg.dat"   => "dir1deg",
         "dir5min.dat"   => "dir5min",
	 "dir30sec.dat"  => "dir30sec");

foreach (keys %dirdat){
  unless(-r "$tdatdir/$_"){
    chdir $srcdir;
    
    unless(-x "$srcdir/bin/$dirdat{$_}"){
      print "Making $dirdat{$_}\n" if( $verbose );
      $makecmd = "makearps $run -opt 3 $dirdat{$_}";
      system($makecmd) == 0 or die "$0: $makecmd failed";
    }
    
    print "Running $dirdat{$_}\n" if( $verbose );
    chdir $wrkdir;
    $runcmd = "$srcdir/bin/$dirdat{$_} < $rname.input > $dirdat{$_}.output";
    ( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";
    
# copy to $tdatdir -- output directory
    ( ! -f "$wrkdir/$_" )
       or (cp( "$wrkdir/$_","$tdatdir/$_" ) == 1
           and cp("$wrkdir/$dirdat{$_}.hdr","$tdatdir/$dirdat{$_}.hdr")==1 )
       or die "$0: cp() failed" ;

  }  
}

#----------- END OF ADDITION --------------------------------------

chdir $srcdir;

unless(-x "$srcdir/bin/$execmd"){
  $makecmd = "makearps $run -opt 3 $execmd";
  system($makecmd) == 0 or die "$0: $makecmd failed";
}

chdir $wrkdir;

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $rname started  at $time";

$runcmd = "$srcdir/bin/$execmd < $rname.input > $rname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $rname finished at $time";

( ! -f "arpstern.dat" )
  or cp( "arpstern.dat", "$trndata" ) == 1
  or die "$0: cp() failed";

# using arpstrn terrain file now:
#( ! -f "$wrkdir/$trndata" )
#  or cp( "$wrkdir/$trndata","$tdatdir/$trndata" ) == 1
#  or die "$0: cp() failed" ;

exit 0
