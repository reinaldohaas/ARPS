#!/usr/bin/perl -w

use Cwd;
use File::Basename;
use File::Copy "cp";
use Getopt::Long;

use lib "./scripts";
use mkinput;

$nx_glob = 67;
$ny_glob = 67;
$nz_glob = 43;

$nx = $nx_glob;
$ny = $ny_glob;
$nz = $nz_glob;

$runname = "arps25may1998";
$subdir = "intrp";
$execmd = "arpsintrp";
$arpsplt = "arpspltpost";

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

$intrpinput = "$srcdir/input/arpsintrp.input";
$pltinput = "$srcdir/input/arpsplt.input";

#
# make the executable arpsplt if necessary
#

unless(-x "$srcdir/bin/$arpsplt"){
  print "Making arpsplt ...\n" if ( $verbose );
  chdir $srcdir;

  $makecmd = "makearps $run -opt 3 $arpsplt";
  system($makecmd) == 0 or die "$0: $makecmd failed";
}

#
# make the executable arpsintrp if necessary
#

unless(-x "$srcdir/bin/$execmd"){
  print "Making arpsintrp...\n" if ( $verbose );
  chdir $srcdir;

  $makecmd = "makearps $run -opt 3 $execmd";
  system($makecmd) == 0 or die "$0: $makecmd failed";
}

#
# set the input parameters for arpsintrp
#

%newpar = (

# &history_data
   hinfmt       => 1,
   hdmpinopt    => 1,
   hdmpfheader  => "'arps25may1998'",
   tintv_dmpin  => 3600.0,
   tbgn_dmpin   => 0.0,
   tend_dmpin   => 21600.0,
   grdbasfn     => "'./arps25may1998.bingrdbas'",

# &output_dims
    nx1 => $nx,
    ny1 => $ny,
    nz1 => $nz,

# &newgrid_soil
   soilscheme => 0,
   nzsoil1    => 2,
   dzsoil1    => 1.0,

# &newgrid
   z_or_p   => 1,
   xy_or_ll => 1,
   xctr1    =>  712000.0,
   yctr1    => 1336000.0,

   snap_to_grid => 0,
   same_res     => 0,

   dx1       =>  20000.0,
   dy1       =>  20000.0,
   dz1       =>  500.0,

   zrefsfc1  => 0.0,
   strhopt1  => 0,

   noutgrds      => 0,
   "xctr_grd(1)" => 0.,
   "yctr_grd(1)" => 0.,
   "clat_grd(1)" => 0.,
   "clon_grd(1)" => 0.,
   "name_grd(1)" => "'NULL'",

# &newterrain
   ternopt1 => 3,
   terndta1 => "'arpstern.data'",
   ternfmt1 => 0,

# &bgloption
   intrphopt => 2,
   intrpvopt => 1,

   bglopt    => 2,
   misvalue  => -9999.0,

# &sfc_data
   sfcdat   => 0,
   sfcdtfl  => "'../jun08/jun08B.sfcdata'",
   sfcfmt   => 0,

# &output
   hdmpfmt  => 1,
   readyfl  => 0,

   sfcout   => 1,
   landout  => 1,
   radout   => 1,
   flxout   => 1,

   qcexout  => 1,
   qrexout  => 1,
   qiexout  => 1,
   qsexout  => 1,
   qhexout  => 1,

   exbcdmp  => 0,
   sfcdmp   => 0,
   soildmp  => 0,
   terndmp  => 0
);

chdir $wrkdir;

#
# set the input parameters for arpsplt
#

%newplt = (

# &history_data
   hinfmt      => 1,
   hdmpinopt   => 1,
   hdmpfheader => "'./arps25may1998'",
   tintv_dmpin => 3600.0,
   tbgn_dmpin  => 0.0,
   tend_dmpin  => 21600.0,
   grdbasfn    => "'arps25may1998.bingrdbas'"
);

#
# Copy data
#

print "Coping data from $srcdir/realcase ..." if ($verbose);

$real = "$srcdir/realcase";
( -d "$real") or die "Can not find real case.";

( -r "$real/arps25may1998.bingrdbas") or die "Can not find base file!";
cp ( "$real/arps25may1998.bingrdbas", "$wrkdir/" ) == 1
                          or die "cp() failed";

$bgn = $newpar{tbgn_dmpin};
$end = $newpar{tend_dmpin};
$inc = $newpar{tintv_dmpin};

while ($bgn <= $end) {
  $fn = sprintf "%s/arps25may1998.bin%06d", $real,$bgn;  
  ( -r "$fn") or die "Can not find file $fn!";
  cp("$fn", "$wrkdir/") == 1 or die "cp() failed";
  $bgn += $inc;
}

#
# Plot real case
#

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run arpspltpost for realcase started at $time";

$newplt{hdmpfheader} = qq('${runname}');
$newplt{grdbasfn}    = qq('${runname}.bingrdbas');

$newinput = "$wrkdir/${runname}_plt.input";
mkinput( \%newplt, $pltinput, $newinput );

$runcmd = "$srcdir/bin/$arpsplt < ${runname}_plt.input > ${runname}_plt.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run arpsplt finished at $time\n";

#
# 1. NW case
#

$rname = "${runname}_nw";
$newpar{runname} = qq('$rname');
$newpar{xctr1}  =  712000.0,
$newpar{yctr1}  = 1336000.0,

$newinput = "$wrkdir/$rname.input";
mkinput( \%newpar, $intrpinput, $newinput );

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $rname started at $time";

$runcmd = "time $srcdir/bin/$execmd < $rname.input > $rname.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

# plot

print "       Run arpsplt\n";
$newplt{hdmpfheader} = qq('${rname}');
$newplt{grdbasfn}    = qq('${rname}.bingrdbas');

$newinput = "$wrkdir/${rname}_plt.input";
mkinput( \%newplt, $pltinput, $newinput );

$runcmd = "$srcdir/bin/$arpsplt < ${rname}_plt.input > ${rname}_plt.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $rname finished at $time\n";

#
# 2. SE case
#

$rname = "${runname}_se";
$newpar{runname} = qq('$rname');
$newpar{xctr1}  = 1336000.0,
$newpar{yctr1}  =  712000.0,

chdir $wrkdir;

$newinput = "$wrkdir/$rname.input";
mkinput( \%newpar, $intrpinput, $newinput );

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $rname started at $time";

$runcmd = "time $srcdir/bin/$execmd < ${rname}.input > ${rname}.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

# Plot

print "       Run arpsplt\n";
$newplt{hdmpfheader} = qq('${rname}');
$newplt{grdbasfn}    = qq('${rname}.bingrdbas');

$newinput = "$wrkdir/${rname}_plt.input";
mkinput( \%newplt, $pltinput, $newinput );

$runcmd = "$srcdir/bin/$arpsplt < ${rname}_plt.input > ${rname}_plt.output";
( $norun ) or system( $runcmd ) == 0 or die "$0: $runcmd failed";

$time = `date +%H:%M' '%a' '%h' '%d' '%Y`;
print "       Run $rname finished at $time\n";

exit 0
