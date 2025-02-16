#!/usr/bin/perl

#
# subroutine updinc is to read in an ARPS dimension file, substitute
# the values of parameters specified in a hash list, %incpar, and then
# save to the same file. The filename and the parameters/values hash
# should be passed through argument list.
#
# Usage: updinc( \%incpar, $incfn );
#
# where
#   %incpar is a hash containing all parameters and their values
#

package updinc;
require Exporter;
@ISA = qw( Exporter );
@EXPORT = qw( updinc );

use File::Copy 'cp';

sub updinc {

  my $incpar = shift @_;
  my $incfn  = shift @_;

  open INCFN, "<$incfn"  or die "$0: cannot open $incfn: $!\n";
  open TMPFN, ">$incfn.tmp"  or die "$0: cannot open $incfn.tmp: $!\n";

  my ( $line, $tmp, $var, $vkey, $val );

  my $update = 0;

  while ( $line = <INCFN> ) {
    chomp ( $tmpline = $line );

    if ( $tmpline =~ /^      \s*?[p|P][a|A][r|R][a|A][m|M][e|E][t|T][e|E][r|R]\s*?\((.*?)$/ ) {

      $tmp = $1;
      print TMPFN "      PARAMETER( ";

      while ( $tmp ) {

        if ( $tmp =~ /^\s*?(\w+?)\s*=\s*(.*?)\s*?,+?(.*?)$/ ) {
          ($var, $val, $tmp) = ($1, $2, $3);
          foreach $dkey (keys %{$incpar} ) {
            if ( $dkey eq $var ) {
              if ( $val ne ${$incpar}{$dkey} ) {
                $val = ${$incpar}{$dkey};
                $update ++;
              }
            }
          }
          print TMPFN "$var = $val, ";
        }     # key = val,
        elsif ( $tmp =~ /^\s*?(\w+?)\s*=\s*(.*?)\s*?\)+?(.*?)$/ ) {
          ($var, $val, $tmp) = ($1, $2, $3);
          foreach $dkey (keys %{$incpar} ) {
            if ( $dkey eq $var ) {
              if ( $val ne ${$incpar}{$dkey} ) {
                $val = ${$incpar}{$dkey};
                $update ++;
              }
            }
          }
          print TMPFN "$var = $val )";
          print TMPFN "$tmp\n";
          $tmp = "";
        }     # key = val )

        $var = "";
        $val = "";
      }
    }
    else {
      print TMPFN "$line";
    }
  }

  close INCFN;
  close TMPFN;

  if ( $update > 0 ) {
    cp( "$incfn.tmp", "$incfn" ) == 1 or die "$0: cp() failed";
    unlink "$incfn.tmp";
  }

  return $update;

}
