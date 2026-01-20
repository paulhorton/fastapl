#!/usr/bin/perl -w
#  Author: Paul Horton
#  Copyright (C) 2011, 2017. Paul Horton, All rights reserved.
#  Created: 20110214
#  Updated: 20170805
#
#  Purpose: Provide basic information regarding Amino acids for
#           convenient use with fastapl.
#
#  Description:  See pod below.
#
package fastaplAmino;
use Exporter 'import';

sub AA();          #  Returns list of standard amino acid characters.
sub AAU();         #  Like AA, but includes U for selenocysteine.
sub prop(*@);      #  prop($name) returns a hash reference representing the given amino acid propensity, such as molecular weight.

@EXPORT  =  qw(AA AAU prop propID propName);
use strict;
use feature qw(say switch);


#  List of standard amino acid characters
sub  AA(){   qw(A C    N D Q E  F G H I L K M P R S T V W Y)   }
sub AAU(){   qw(A C U  N D Q E  F G H I L K M P R S T V W Y)   }


my %KYTJ820101 =
    ( name => 'Kyte Doolittle',
      id   => 'KYTJ820101',
      A =>  1.8, R => -4.5, N => -3.5, D => -3.5, B => -3.5, C =>  2.5,
      Q => -3.5, E => -3.5, Z => -3.5, G => -0.4, H => -3.2, I =>  4.5, L =>  3.8, K => -3.9,
      M =>  1.9, F =>  2.8, P => -1.6, S => -0.8, T => -0.7, W => -0.9,
      Y => -1.3, V =>  4.2, X => -0.19 );


#  Computed by P.H. based on WoLF PSORT dataset of proteins from animal, plant, fungi with subcellular localization site annotation.
my %eukaryoticAbundance =
    ( name => 'Typical Eukaryotic Abundance',
      id   => 'euAbundance',
      A => 0.072,
      C => 0.020,
      U => 0.0001,   #  selenocysteine, absent in plants and fungi.
      N => 0.042,
      D => 0.049,
      Q => 0.041,
      E => 0.0629,
      F => 0.043,
      G => 0.066,
      H => 0.023,
      I => 0.053,
      K => 0.059,
      L => 0.096,
      M => 0.025,
      P => 0.055,
      R => 0.052,
      S => 0.079,
      T => 0.055,
      V => 0.064,
      W => 0.012,
      Y => 0.031,
);


my %molWeight =
    ( _name => 'Molecular weight',
      _id   => 'molWeight',
      A =>  89.1,
      R => 174.2,
      N => 132.1,
      D => 133.1,
      B => 132.6,
      C => 121.2,
      U => 168.1,
      Q => 146.2,
      E => 147.1,
      Z => 146.6,
      G =>  75.1,
      H => 155.2,
      I => 131.2,
      L => 131.2,
      K => 146.2,
      M => 149.2,
      F => 165.2,
      P => 115.1,
      S => 105.1,
      T => 119.1,
      W => 204.2,
      Y => 181.2,
      V => 117.1,
      X => 129.5,
    );



my %toProp  =
    (
      'euabundance'    => \%eukaryoticAbundance,
      'kyte doolittle' => \%KYTJ820101,
      'kytj820101'     => \%KYTJ820101,
      'kd'             => \%KYTJ820101,
      'molweight'      => \%molWeight,
    );


sub propID{
    my $key  =  lc $_[0];
    return ''  unless exists $toProp{$key};
    return ''  unless exists ${ $toProp{$key} }{_id};
    return  ${ $toProp{$key} }{_id};
}


sub propName{
    my $key  =  lc $_[0];
    return ''  unless exists $toProp{$key};
    return ''  unless exists ${ $toProp{$key} }{_name};
    return  ${ $toProp{$key} }{_name};
}


#  Return hash mapping amino acid one-letter code to property $propName
#  $propName can be abbreviated as a prefix of a prop name if it is long enough to match only one prop.
#  $propName may optionally be followed by a list of keys.
#    e.g.  prop 'molWeight', qw(R K H)  returns a hash mapping R,K,H to the molecular weights of Arginine, Lysine and Histidine.
#  Prototype allows bareword argement as in:  prop molWeight, qw(R K H);
sub prop(*@){
    my ($propNameOrPrefix, @aaSubset)  =  @_;

    my $propName   =  _propNameCheck( $propNameOrPrefix );
    my $propRf  =  $toProp{$propName};

    @aaSubset   or   return %$propRf;


    my %propSubset;
    for( @aaSubset ){
        my $key  =  uc $_;
        exists $propRf->{$key}   or   die  "Propensity '$propName', has no (amino acid) key '$key'";
        $propSubset{$key}  =  $propRf->{$key};
    }

    return %propSubset;
}



#  Return hash mapping amino acid one-letter code to property $propName.
#  _propNameCheck( propName )
#  If propName matches (wholly or by a unique prefix) a %toProp key, return that key.
#  otherwise die with error.
#
#  Technical detail: if multiple names match, but all point to an identical value in %toProp,
#                    the longest of those names is returned.
sub _propNameCheck(*){
    my $propName   =  $_[0];

    #  Canonicalize propensity name to allow matching insensive to case and white space.
    $propName  =  lc $propName;  $propName =~ s/\s//g;


    return $propName    if exists $toProp{$propName};   #<---  EARLY EXIT POINT.


    #  $propName is not an exact match, so check for prefix mathing.
    my @matchingKey;  #To hold matching propNames;


    my %toProp_valSeen;
    for  my $toProp_key  (sort {length $a <=> length $b} keys %toProp){
        if(   $propName  eq  lc substr( $toProp_key, 0, length $propName )   ){
            $toProp_valSeen{ $toProp{$toProp_key} }++   or   push @matchingKey, $toProp_key;
        }
    }

    @matchingKey   or   die  "No fastaplAmino prop matching '$propName'";

    @matchingKey == 1   or
        die  "Ambiguous prefix match to fastaplAmino prop; '$propName' matches { @matchingKey }";

    return  $matchingKey[0];
}



1;



=pod


=head1  NAME

fastaplAmino


=head1  SYNOPSIS

  #  print hydrophobicity of Glycine.
  my %h  =  $properties{'Kyte Doolittle'};
  print  $h{'G'}, "\n";


=head1  DESCRIPTION

Module to provide some functions on amino acids sequences useful
with B<fastapl>.  Currently very minimal.  Will probably be
expanded in the future.


=head1  LICENSE

You may use this program under the GNU public license.


=head1  AUTHOR

Paul Horton <paulh@iscb.org>


=head1  COPYRIGHT

Copyright (C) 2011.

=cut
