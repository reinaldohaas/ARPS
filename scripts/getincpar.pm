#!/usr/bin/perl -w

#
# subroutine getincpar is to get parameters from an ARPS include file.
#
# Usage: getincpar( \%pars, $incfn );
#
# where
#   %pars  - a hash containing all parameter keys and their values
#   %incfn - filename of ARPS include file
#

package getincpar;
require Exporter;
@ISA = qw( Exporter );
@EXPORT = qw( getincpar );

use File::Copy 'cp';

sub getincpar {

  local ( *vkeys ) = shift @_;
  my $incfn = shift @_;

  my %pars = ();

#print "$incfn\n";
  open INCFN, "<$incfn"  or die "$0: cannot open $incfn: $!\n";

  my ( $line, $tmp, $var, $val );

  while ( $line = <INCFN> ) {
    chomp ( $tmp = $line );

    if ( $tmp =~ /^     \s+?[p|P][a|A][r|R][a|A][m|M][e|E][t|T][e|E][r|R]\s*?\((.*?)$/ or
         $tmp =~ /^     \S\s*?(.*?)$/ ) {

      $tmp = $1;
#print "tmp0 = $tmp\n";

      while ( $tmp ) {
        if ( $tmp =~ /^\s*?(\w+?)\s*?=\s*?(.*?)\s*?\)(.*?)$/ ) {
          ($var, $val, $tmp) = ($1, $2, $3);
#print "tmp1 = $tmp, var = $var, val = $val\n";
          foreach $k ( @vkeys ) {
            if ( $k eq $var ) {
              $pars{$k} = $val;
#print "$k = $val\n";
            }
          }
          $tmp = "";
        }     # key = val )
        elsif ( $tmp =~ /^\s*?(\w+?)\s*?=\s*?(.*?)\s*?,(.*?)$/ ) {
          ($var, $val, $tmp) = ($1, $2, $3);
#print "tmp2 = $tmp, var = $var, val = $val\n";
          foreach $k ( @vkeys ) {
            if ( $k eq $var ) {
              $pars{$k} = $val;
#print "$k = $val\n";
            }
          }
        }
        else {
          $tmp = "";
        }

        $var = "";
        $val = "";
      }
    }
  }

  close INCFN;

  return %pars;
}
