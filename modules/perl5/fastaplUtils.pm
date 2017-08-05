#!/usr/bin/perl -w
#  Author: Paul Horton
#  Copyright (C) 2014, 2017.  Paul Horton, All rights reserved.
#  Creation Date: 2014.3.8
#  Last Modification: 20170802
#
#  Description: utils used by fastapl which might be useful in other contexts as well.
#
package fastaplUtils;
use Exporter 'import';

use feature 'say';

@EXPORT  =  qw(astr astrf hstr hstrf maxidx maxSeg maxSegs say2 printinHtml printinXterm);

=pod

=head1 NAME

fastaplUtils - utils used by fastapl which might be useful in other contexts as well.

=head1 SYNOPSIS

  sub astr;          #Convert list or array to a tab separated string.
  sub astrf;         #Like astr(), but with customizable separator given by the first arg.
  sub hstr;          #Convert hash reference or even-length list to string.
  sub hstrf;         #Like hstf(), but with customizable format stipulated by first arg.
  sub maxidx;        #Returns argmax of list; or in list context (argmax, max).
  sub maxSeg;        #Returns a maximal sum interval of @_, or (0,0,0).
  sub maxSegs;       #Returns the maximal sum intervals of @_ or ().
  sub say2;          #say to STDERR, and also STDOUT unless both are connected to a terminal.
  
  #  Print with markup; e.g. printin( 'blue,ul', @arg ) prints "@arg" in blue, underlined text.
  sub printinHtml;   #printin using HTML tags
  sub printinXterm;  #printin xterm escape cods

=cut

sub astr;          #Convert list or array to a tab separated string.
sub astrf;         #Like astr(), but with customizable separator given by the first arg.
sub hstr;          #Convert hash reference or even-length list to string.
sub hstrf;         #Like hstf(), but with customizable format stipulated by first arg.
sub maxidx;        #Returns argmax of list; or in list context (argmax, max).
sub maxSeg;        #Returns a maximal sum interval of @_, or (0,0,0).
sub maxSegs;       #Returns the maximal sum intervals of @_ or ().
sub say2;          #say to STDERR, and also STDOUT unless both are connected to a terminal.

#  Print with markup; e.g. printin( 'blue,ul', @arg ) prints "@arg" in blue, underlined text.
sub printinHtml;   #printin using HTML tags
sub printinXterm;  #printin xterm escape cods


#  Convert LIST or ARRAY to a tab separated string.
sub astr{   astrf "\t\n", @_;   }


#  Convert LIST or ARRAY possibly including references to other arrays (which are not allowed to contain references themselves)
#  to a string according to FORMAT.
#  Ordinary elements are separated by $sep1
#  Array references are expanded to $sep1 separated strings, which are separated from other elements by $sep2.
#  $sep1 is the first char of FORMAT and $sep2 is the remainder.
sub astrf{
    defined wantarray   or   die  'useless call of astrf in void context';
    @_   or   die  'astrf() expected at least one arg but got none';
    my ($format, @elem)  =  @_;

    my ($sep1, $sep2)  =  $format =~ /^(.)(.*)$/s
        or   die  "astrf() expected first arg to be a string of 2 or more characters, but got '$format'";
    @elem   or   return '';  #Empty list, so return empty string.

    #  If a @elem is a single hash reference, create string  key1 . $sep . val1 . $sep key2...
    if(  ref $elem[0]  eq  'HASH'  ){
        my @a  =  %{$elem[0]};
        return  _astrf( $sep1, \@a );
    }

    my $retVal  =  '';

    my $count  =  0;
    my $prevWasArray;
    for my $elem (@elem){
        my $elemType  =  ref $elem;
        if(  $elemType  ){
            $elemType eq 'ARRAY'  or   die  "astrf() received reference to '$elemType' which it cannot handle";
            $retVal  .=  $sep2    if $count;
            $retVal  .=  _astrf( $sep1, $elem );
            $prevWasArray = 1;
        }
        else{
            my $sep  =  $prevWasArray? $sep2 : $sep1;
            $retVal  .=  $sep    if $count;
            $elem  =  'undef'   if !defined $elem;
            $retVal  .=  $elem;
            $prevWasArray = 0;
        }
        ++$count;
    }

    return $retVal;
}#END: astrf()


#  Auxillary function for astrf, should receive exactly two args ($sep, $aRef).
sub _astrf{
    my ($sep, $aRef)  =  @_;

    my $retVal  =  '';
    for my $elem (@$aRef){
        $retVal .=  $sep     unless  $retVal eq '';
        $retVal .=  $elem // 'undef';
    }
    return $retVal;

}# END: _astrf()



#  Convert HASH to a string as key1\tval1\nkey2\val2...
sub hstr{    hstrf "\t\n", @_;    }


#  hstrf( FORMAT, HASH ) converts HASH to string as stipulated by FORMAT.
#  String is formatted as <key1 sep1 val1 sep2 key2 sep1 val2...>
#  The first character of FORMAT is $sep1, the rest is $sep2.
#  HASH may take the form of a (even length) list, or a reference to a hash or (even length) array.
sub hstrf{
    @_   or   die  'hstrf expected at least one arg but got none';
    my $format  =  shift;
    (my $sep1, my $sep2)  =  $format =~ /^(.)(..?)$/s
        or   die  "hstrf() expected first arg to be a string of 2 or more characters, but got '$format'";

    @_   or   return;  # empty hash, so just return.
    my $argType  =  ref $_[0];

    my $hRef;  #  to hold reference to hash.

    #  Set $hRef from @_, details depend on $argType.
    #
    #  ─────  Passed list  ─────
    if(  !$argType  ){
        die  'In hstrf(), expected hash, but got list with odd number of elements'   if @_ % 2;
        my %h  =  @_;
        $hRef  =  \%h;
    }
    #  ─────  Passed hash reference  ─────
    elsif(  $argType  eq  'HASH'  ){
        @_ == 1
            or   die  'hstrf() only expected one argument when its first is a hash reference, but got ', 0+@_, ' args.';
        $hRef  =  shift;
    }
    #  ─────  Passed array reference  ─────
    elsif(  $argType  eq  'ARRAY'  ){
        @_ == 1
            or   die  'hstrf() only expected one argument when its first is an array reference, but got ', 0+@_, ' args.';
        my %h  =  @$_;
        $hRef  =  \%h;
    }    
    else{
        die  "hstrf() did not expect its first argument to be a reference of type '$argType'.";
    }

    #  Sort like: AA, Aa, aA, aa, AC, Ac, ...
    my @key   =   sort { lc $a cmp lc $b  ||  $a cmp $b  }   keys %$hRef;

    my $retVal  =  '';
    for my $key (@key){
        $retVal .=  $sep2    unless $retVal eq '';
        $retVal .=  $key . $sep1 . $$hRef{$key};
    }
    return $retVal;

}#END: hstrf()


#  Returns argmax of input list; or in list context (argmax, max).
#  The smallest index is returned in the case of ties.
sub maxidx{
    @_   or   die 'maxidx called with empty list';
    my $maxIdx =  0;
    my $max    =  $_[0];

    for ( 1..$#_ ) {
        ($maxIdx, $max)  =  ($_, $_[$_])   if  $_[$_] > $max;
    }

    wantarray?   ($maxIdx, $max)  :  $maxIdx;
}#END: maxidx()



#  Returns an interval (begin, end, maxSum),
#  such that [begin, end] give the closed interval of a maximum positive sum (= maxSum) segment of @_.  Ties go to the left-most interval.
#
#  Example illustrating ties:
#
#      maxSeg  0 +2 -6 +2 -2 +2;   #  returns (1, 1, 2).
#
#  Returns (0, 0, 0) when no positive elements are found in @_.
sub maxSeg{

    wantarray   or   die  'maxSegment() called in non-list context';
    @_          or   die  return 0, 0, 0;   #Empty list.

    my $maxIdx;
    my ($maxSum, $sum)  =  (0, 0);
    for my $i (0..$#_) {
        $sum  +=  $_[$i];
        if(  $maxSum < $sum  ){   ($maxIdx, $maxSum)  =  ($i, $sum)   }
        elsif(     0 > $sum  ){                $sum   =   0           }
    }

    $maxSum > 0   or   return 0, 0, 0;


    #  At this point,
    #  $maxIdx should hold the rightmost index of $maxSum, the maximum value in @_

    $sum  =  0;
    for my $i (reverse 0..$maxIdx) {
        $sum += $_[$i];
        return  $i, $maxIdx, $maxSum    if  $sum >= $maxSum;
    }

    die  'code should never reach here!';
}#END: maxSeg()


#  Returns a list of minimal length maximal positive scoring segments -- i.e. returns a list of all segments such that
#    1)  The sum is positive
#    2)  No extension can increase the sum
#    3)  No retraction can attain an equal sum, i.e. zero-sum ends are not included.
#
#  Each element of the list is [$begin, $end, $maxSum] where [begin, end] are indices of a maximum sum segment of @_ and $maxSum its sum.
#  When no positive elements exist in @_ the empty list is returned.
sub maxSegs{

    wantarray   or   die  'maxSegment() called in non-list context';
    @_          or   return ();  #Empty list.

    my @max;
    my @retVal;

    my ($maxIdx, $maxSum, $sum)  =  (undef, 0, 0);

    #  Push all positive maximal segments, except the one at the end of @_
    for  my $i  ( 0..$#_ )  {
        $sum  +=  $_[$i];
        if(  $maxSum <  $sum  ){
            ($maxSum, $maxIdx)  =  ($sum, $i);
        }
        elsif(   0  >=  $sum  ){
            #Sum went negative, so save info for most recent maxSeg and reset $sum, $maxSum.
            push  @max,  [$maxIdx, $maxSum]    if  $maxSum > 0;
            $sum  =  $maxSum  =  0;
        }
    }

    push  @max,  [$maxIdx, $maxSum]    if  $maxSum > 0;

    #  At this point, we know the sum and right hand end of each maximal region.
    #  So we work back from the end to find the start.
    for (@max) {
        ($maxIdx, $maxSum)  =  @$_;

        my $i  =  $maxIdx;
        for(  $sum  =  $_[$i];
              $sum != $maxSum;
              $sum  +=  $_[--$i]
            ){
         }

        push @retVal, [$i, $maxIdx, $maxSum];
    }

    return @retVal;

    die  'Execution should never reach here.';
}#END: maxSegs()



#  Print string in @_ to STDERR and also STDOUT unless both are connected to a terminal.
sub say2{
    for( @_ ){
        say STDERR;
        say STDOUT  unless  -t *STDERR && -t *STDOUT;
    }
}



#  ──────────  printin functions  ──────────

my %colorName_RGB =
    qw(black 000000  red CD0000  green 00CD00  yellow CDCD00  blue 0000CD  magenta CD00CD  cyan 00CDCD  gray E5E5E5  red_bold FF0000  green_bold 00FF00  yellow_bold FFFF00  blue_bold 5C5CFF  magenta_bold FF00FF  cyan_bold 00FFFF  white FFFFFF);

my %nameToXColorCode  =  qw(default 0  black 30  red 31  green 32  yellow 33  blue 34  magenta 35  cyan 36  gray 37);


sub printinHtml{
    my $markupSpec  =  lc shift;
    my @printArg    =  @_;

    my ($opn, $clz)  =   ('','');
    my @directive  =  split ',' => $markupSpec;
    for( @directive ){
        if( /blink/         ){$opn .= '<BLINK>';  $clz = '</BLINK>' . $clz;  next}
        if( /bold/          ){$opn .= '<B>'    ;  $clz = '</B>'     . $clz;  next}   
        if( /it(?:alics?)?/ ){$opn .= '<I>'    ;  $clz = '</I>'     . $clz;  next}
        if( /iv|inverse/    ){                                               next}
        if( /ul|underline/  ){$opn .= '<U>'    ;  $clz = '</U>'     . $clz;  next}   
        if( /plain/         ){                                               next}

        my $colorHex;
        if(  exists $colorName_RGB{$_}  ){
            $colorHex  =  $colorName_RGB{$_};
        }else{
            ($colorHex)  =  /^\#?([0-9A-F]{6})$/i
                or   die  "trouble parsing markup spec '$markupSpec', term '$_' in sub printinHtml";
        }

        $opn .= qq(<FONT color="#$colorHex">);
        $clz  =  '</FONT>' . $clz;
    }

    print  "$opn@printArg$clz";
}#END  printinHtml


sub printinXterm{
    my $markupSpec  =  lc shift;
    my @printArg    =  @_;

    my ($opn, $clz)  =   ('','');
    my @directive  =  split ',' => $markupSpec;
    for( @directive ){
        if( /blink/         ){$opn .= "\e[5m";  $clz = "\e[25m" . $clz;  next}
        if( /bold/          ){$opn .= "\e[1m";  $clz = "\e[0m"  . $clz;  next}   
        if( /it(?:alics?)?/ ){                                           next}
        if( /iv|inverse/    ){$opn .= "\e[7m";  $clz = "\e[27m" . $clz;  next}
        if( /ul|underline/  ){$opn .= "\e[4m";  $clz = "\e[24m" . $clz;  next}   
        if( /plain/         ){                                           next}
        if(  exists $nameToXColorCode{$_}  ){
            $opn  .=  "\e[" . $nameToXColorCode{$_} . 'm';
            $clz   =  "\e[39m" . $clz;
        }else{
            die  "sub printinXterm does not support color in hexadecimal"    if /^\#?([0-9A-F]{6})$/i;
            die  "trouble parsing markup spec '$markupSpec', term '$_' in sub printinXterm";
        }
    }

    print  "$opn@printArg$clz";
}#END  printinXterm

1;


=pod

=head1 AUTHOR

Paul Horton
