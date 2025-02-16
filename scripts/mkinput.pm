#!/usr/bin/perl -w

#
# subroutine mkinput is to read in an ARPS namelist file, substitute
# the values of parameters specified in a hash list, newpar, and then
# save to a new file. The input and output filenames and the
# parameters/values pair hash should be passed through argument list.
#
# Usage: mkinput( \%newpar, $infn, $outfn );
#
# where
#        %newpar is a hash containing all parameters and their 
#

package mkinput;
require Exporter;
@ISA = qw( Exporter );
@EXPORT = qw( mkinput );

sub mkinput {

  my $newpar = shift @_;
  my $infn  = shift @_;
  my $outfn = shift @_;  
  my $chgpar = {};

  open INFN,  "<$infn"  or die "$0: cannot find $infn: $!\n";
  open OUTFN, ">$outfn" or die "$0: cannot find $outfn: $!\n";

  my ( $line, $tmp, $var, $vkey, $val );

  while ( defined (INFN) && ($line = <INFN>) ) {
    chomp ( $tmp = $line );

#   print "tmp = $tmp\n";
    if ( $tmp =~ /^\s+?                           # start with whitespaces
                  (\w+\(??\s*?\d*?\s*?,??\s*?\d*?\s*?\)??)     # keyword(*)
                  \s*?=\s*?(.*?)$/x ) {                        # = whatever
      print OUTFN "   ";

#     print "tmp = $tmp\n";
      while ( $tmp ) {

        ($var, $tmp) = ($1, $2) if (
          $tmp =~ /^\s*?
                   (\w+\(??\s*?\d*?\s*?,??\s*?\d*?\s*?\)??)   # \1
                   \s*?=\s*?
                   (.*?)$/x );                                # \2

# capture pattern /'text'/
        ($val, $tmp) = ($1, $2) if ( $tmp =~ /^\s*?('.*?'),(.*?)$/ );

# capture pattern /.text./
        ($val, $tmp) = ($1, $2) if ( $tmp =~ /^\s*?(\.\w+\.),(.*?)$/ );

# capture a number, e.g., 123, or 1.0, or 1, or 1.e-1, or 0.1-e02.
        ($val, $tmp) = ($1, $2)
          if ( $tmp =~ /^\s*?([-+]?[0-9]*?\.?[0-9]*?[eE]?[-+]?[0-9]*?),(.*?)$/ );

# capture pattern /val,val, ...,/
        ($val, $tmp) = ("$val, $1", $2) while ( $tmp =~ /^\s*([^=,()]+)\s*,(.*?)$/ );

# capture pattern /text\s+$/
        ($val, $tmp) = ($1, $2) if ( $tmp =~ /^\s*(\W*?\w+?\W*?\w*?)(\s*)$/ );

#        print "var = $var,	value = $val\n";
        foreach $vkey (keys %{$newpar} ) {
          if ( $vkey eq $var ) {
            $val = ${$newpar}{$vkey};
#           print "key = $vkey,	value = $val\n";
          }
        }
        $var = $1 if ( $var =~ /^\s*?(.*?)\s*?$/ );
        $val = $1 if ( $val =~ /^\s*?(.*?)\s*?$/ );
#        print "var = $var,	value = $val\n";
        print OUTFN "$var = $val,";

        $var = "";
        $val = "";
        $tmp = $1 if ( $tmp =~ /\s*?(.*?)\s*?$/ );
        $tmp = "" if ( $tmp =~ /^\s*?$/ );
#       print "tmp = $tmp\n";
      }
      print OUTFN "\n";
    }
    elsif ( $tmp =~ /^\s+\&end/i ) {
      print OUTFN "$line\n";
    }
    elsif ( $tmp =~ /^\s+\&\w+/ ) {
      print OUTFN "$line";
    }
    else {
      print OUTFN "$line";
    }
  }

  close OUTFN;

  return $chgpar;

}
