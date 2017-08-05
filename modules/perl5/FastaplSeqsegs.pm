#!/usr/bin/perl -w
#  Author: Paul Horton
#  Copyright (C) 2014,2016,2017  Paul Horton, All rights reserved.
#  Creation Date: 2014.5.9
#  Last Modification: 20170802
#
#  Description:  Class to hold an array of sequence segments.
#
#  Context: fastapl
#
package FastaplSeqsegs;
use Carp  'confess';
use feature  'say';
use strict;

sub SIZ(){  0  }  #Number of coords held; or equivalently, twice the number of segments held.
sub NAM(){  1  }  #Name of this collection of matches.                                       
sub ITR(){  2  }  #Number for built-in iterator ab();                                        
sub VEC(){  3  }  #Vec string holding data.                                                  


$::seq    if 0;   #Suppress used only once warnings

use overload
    bool   =>  sub{  ${$_[0]}[SIZ]  },
    '""'   =>  \&asString
    ;


sub new{
    my( $class, $name, @coord )  =  @_;

    $name //= '';

    @coord % 2
        and   die  "When attempting to construct fastaplSeqSegs object named '$name', expected an even sized list of (up, dw) pairs, but got list of size ", 0+@coord;

    !@coord   or
        $coord[0] =~ /^-?\d+$/   or   die  "When attempting to construct fastaplSeqSegs object named '$name', expected integer argument but got '$coord[0]'";


    my $size  =  @coord / 2;
    my $itr = 0;
    my $vec  =  pack 'l*', @coord;

    my $self  =  [$size, $name, $itr, $vec];
    return  bless $self;
}


sub name{  ${$_[0]} [NAM]  }
sub size{  ${_[0]} [SIZ]      }
sub  fin{  ${_[0]} [SIZ] - 1  }


#  Forward Iterator.
sub itr{
    my( $self, $i  )  =  ( $_[0], 0 );

    return sub{
        defined $i   or   return ();  #return empy list, false in list context.

        $i < $$self[SIZ]   or   do{  undef $i;  return};

        return    unpack 'll',  substr $$self[VEC], 8*$i++, 8;
    }

}#  itr()


sub rangeItr{  $::seqIsCirc? goto &rangeItr_circular : goto &rangeItr_linear  }

sub rangeItr_linear{
    my( $self, $i, $L  )  =  ( $_[0], 0, ::len() );


    return sub{
        defined $i   or   confess  'Attempted to use exhausted rangeItr iterator';
        $i < $$self[SIZ]   or   do{  undef $i;  return};

        my( $u, $d )  =   unpack 'll',  substr $$self[VEC], 8*$i++, 8;

        if(  $d < $u  ){
            if(  0 <= $d  ){                      #   0 < $d < $u
                +$L >= $u  or  goto DIE_RANGE;    return 1,     $d,    $u-1
            }
            if(  0 >= $u  ){                      #  $d < $u ≦ 0
                -$L <= $d  or  goto DIE_RANGE;    return 1,  $L+$d, $L+$u-1
            }
            goto DIE_WRAPD;    #  $d < 0 < $u
        }{#else $u ≦ $d
            if(  0 <= $u  ){                      #   0 ≦ $u ≦ $d
                +$L >= $d  or  goto DIE_RANGE;    return 0,     $u,    $d-1
            }
            if(  0 >= $d  ){                      #  $u ≦ $d ≦ 0
                -$L <= $u  or  goto DIE_RANGE;    return 0,  $L+$u, $L+$d-1
            }
            goto DIE_WRAPD;    #  $u < 0 < $d
        }

      DIE_RANGE:
        die  'In ', (caller(1))[3], " handling seqsegs '", $self->name(), "' the ${i}th seqseg ($u,$d) is incompatible with \$seq length $L\n";
      DIE_WRAPD:
        die  'In ', (caller(1))[3], " handling seqsegs '", $self->name(), "' the ${i}th seqseg ($u,$d) crosses zero, but \$seq is not circular\n";
    }#END sub to return
}#END rangeItr_linear

sub rangeItr_circular{
    my( $self, $i, $L  )  =  ( $_[0], 0, ::len() );


    return sub{
        defined $i   or   confess  'Attempted to use exhausted rangeItr iterator';
        $i < 2*$$self[SIZ]   or   do{  undef $i;  return};

        my( $u, $d );
        if(  $i % 2  ){
            ( $u, $d )  =   unpack 'll',  substr $$self[VEC], 4*($i-1), 8;
            $i++;
            return 1, $L+$d, $L-1    if $d < $u;
            return 0, $L+$u, $L-1
        }

        ( $u, $d )  =   unpack 'll',  substr $$self[VEC], 4*$i, 8;
        if(  $d < $u  ){
            if(  0 <= $d  ){                            #   0 < $d < $u
                +$L >= $u  or  goto DIE_RANGE;  $i+=2;  return 1,     $d,    $u-1
            }
            if(  0 >= $u  ){                            #  $d < $u ≦ 0
                -$L <= $d  or  goto DIE_RANGE;  $i+=2;  return 1,  $L+$d, $L+$u-1
            }
            $L >= $u-$d    or  goto DIE_RANGE;  $i++;   return 1,     0,     $u-1
        }{#else $u ≦ $d
            if(  0 <= $u  ){                            #   0 ≦ $u ≦ $d
                +$L >= $d  or  goto DIE_RANGE;  $i+=2;  return 0,     $u,    $d-1
            }
            if(  0 >= $d  ){                            #  $u ≦ $d ≦ 0
                -$L <= $u  or  goto DIE_RANGE;  $i+=2;  return 0,  $L+$u, $L+$d-1
            }
            $L >=  $d-$u   or  goto DIE_RANGE;  $i++;   return 0,      0,    $d-1
        }

      DIE_RANGE:
        die  'In ', (caller(1))[3], " handling seqsegs '", $self->name(), "' the ${i}th seqseg ($u,$d) is incompatible with \$seq length $L\n";
    }#END sub to return
}#END: rangeItr_circular


#  Built in iterator.
#  Set $a, $b to next coordinate pair and return true
#  Unless all pairs have been iterated already, in which case reset and return 0.
sub ab{
    my $self  =  $_[0];

    if(  $$self[ITR] == $$self[SIZ]  ){
        $::a = $::b =  undef;
        return  $$self[ITR] = 0;
    }

    ($::a, $::b)  =   unpack 'll',  substr $$self[VEC], 8*$$self[ITR], 8;
    
    return  ++$$self[ITR];
}



#  Return seqSeg $i strand (0 or 1) followed by one or two ranges covering the $seq indices in seqSeg.
#  Circular sequences return two ranges if the seqSeg crosses zero, otherwise one.
#  A range is a pair of non-negative integers.
#
#    for( 0..$seqSegs->fin() ){
#        my( $isRev, @ranges )  =  $seqSegs->ranges( $_ );
#        while(  ($a, $b)  =  splice @ranges, 0, 2  ){
#             if( $isRev ){   vec( $$revVec, $_, 1 )  =  1   for  $a..b   }
#             else        {   vec( $$fwdVec, $_, 1 )  =  1   for  $a..b   }
#        }
#    }
#
#  Does ranges checking ($seq length may have changed).
sub ranges{
    my( $self, $i )  =  @_;

    #  Short var names increase readability below.
    my( $L, $u, $d )  =   (   ::len(),   unpack 'll',  substr $$self[VEC], 8*$i, 8   );

    if(  $d < $u  ){
        if(  0 <= $d  ){                      #   0 < $d < $u
            +$L >= $u  or  goto DIE_RANGE;    return 1,     $d,    $u-1
        }
        if(  0 >= $u  ){                      #  $d < $u ≦ 0
            -$L <= $d  or  goto DIE_RANGE;    return 1,  $L+$d, $L+$u-1
        }
        $::seqIsCirc   or  goto DIE_WRAPD;    #  $d <  0 < $u
        $L >=  $u-$d   or  goto DIE_RANGE;    return 1,      0,    $u-1    ,$L+$d,  $L-1
    }{#else $u ≦ $d
        if(  0 <= $u  ){                      #   0 ≦ $u ≦ $d
            +$L >= $d  or  goto DIE_RANGE;    return 0,     $u,    $d-1
        }
        if(  0 >= $d  ){                      #  $u ≦ $d ≦ 0
            -$L <= $u  or  goto DIE_RANGE;    return 0,  $L+$u, $L+$d-1
        }
        $::seqIsCirc   or  goto DIE_WRAPD;    #  $u <  0 < $d
        $L >=  $d-$u   or  goto DIE_RANGE;    return 0,      0,    $d-1    ,$L+$u,  $L-1
    }

  DIE_RANGE:
    die  'In ', (caller(1))[3], " handling seqsegs '", $self->name(), "' the ${i}th seqseg ($u,$d) is incompatible with \$seq length $L\n";

  DIE_WRAPD:
    die  'In ', (caller(1))[3], " handling seqsegs '", $self->name(), "' the ${i}th seqseg ($u,$d) crosses zero, but \$seq is not circular\n";
}#END: ranges



sub widen{
    my $self  =  $_[0];
    @_ < 4   or   die  "widen expected at most 2 arguments, but was called with '@_'";
    my $upMargin  =  @_ > 1?  $_[1] : 5;
    my $dwMargin  =  @_ > 2?  $_[2] : $upMargin;

    $upMargin >= 0   or   die  "widen does not currently allow negative arguments, but got '($_[1], $_[2])'";
    $dwMargin >= 0   or   die  "widen does not currently allow negative arguments, but got '($_[1], $_[2])'";

    die  "Error: widen('@_'); negative margins not currently allowed"   if  $upMargin < 0 || $dwMargin < 0;

    my $len  =  length $::seq;
    my( $newUp, $newDw );

    if(  $::seqIsCirc  ){
        #  Circular so $newUp, $newDw can be of mixed sign. 
        #  But both should be ∈ [-len,+len] with |$newUp-$newDw| ≦ len
        my $lenHalf  =  $len / 2;
        for my $i (0..$$self[SIZ]-1){
            my( $up, $dw )  =   unpack 'll',  substr $$self[VEC], 8*$i, 8;
            if(  $up > $dw  ){
                $newUp  =  $up + $upMargin;
                $newDw  =  $dw - $dwMargin;
                if(  $newUp - $newDw > $len  ){
                    my $center  =  ($newUp + $newDw) / 2;
                    $center  -=  $len * int $center/$len;
                    my $newUp  =   ceil( $center + $lenHalf );
                    my $newDw  =  floor( $center - $lenHalf );
                }
                else{
                    if( $newUp > +$len ){  $_ -= $len  for $newUp, $newDw  }
                    if( $newDw < -$len ){  $_ += $len  for $newUp, $newDw  }
                }
            }else{#  $up ≦ $dw
                $newUp  =  $up - $upMargin;
                $newDw  =  $dw + $dwMargin;
                if(  $newDw - $newUp > $len  ){
                    my $center  =  ($newUp + $newDw) / 2;
                    $center  -=  $len * int $center/$len;
                    my $newUp  =  floor( $center + $lenHalf );
                    my $newDw  =   ceil( $center - $lenHalf );
                }
                else{
                    if( $newDw > +$len ){  $_ -= $len  for $newUp, $newDw  }
                    if( $newUp < -$len ){  $_ += $len  for $newUp, $newDw  }
                }
            }
            substr  $$self[VEC], 8*$i, 8  =>  pack 'll', $newUp, $newDw;
        }#END for $i
    }
    else{
        #  Linear sequence.  $up, $dw should have the same sign.
        #  Widening is not allowed to change the sign of coordinates or exceed length of sequence.
        for my $i (0..$$self[SIZ]-1){
            my( $up, $dw )  =   unpack 'll',  substr $$self[VEC], 8*$i, 8;

            if(  $up > $dw  ){
                $newUp  =  $up + $upMargin;
                $newDw  =  $dw - $dwMargin;
                my $floor =  $dw <  0?    -$len  :    0;
                my $ceil  =  $floor + $len;
                $newUp = $ceil    if $newUp > $ceil;
                $newDw = $floor   if $newDw < $floor;
            }else{#  $up ≦ $dw
                $newUp  =  $up - $upMargin;
                $newDw  =  $dw + $dwMargin;
                my $floor =  $up <  0?    -$len  :    0;
                my $ceil  =  $floor + $len;
                $newUp = $floor   if $newUp < $floor;
                $newDw = $ceil    if $newDw > $ceil;
            }
            substr  $$self[VEC], 8*$i, 8  =>  pack 'll', $newUp, $newDw;
        }
    }
    return $self;
}#  widen( upMargin, dwMargin )



#  Transfer all reverse strand segments to forward strand, in place.
sub allToFwd{
    my $self  =  $_[0];

    for my $i (0..$$self[SIZ]-1){
        my( $up, $dw)  =   unpack 'll',  substr $$self[VEC], 8*$i, 8;
         substr  $$self[VEC], 8*$i, 8  =>  pack 'll', $dw, $up    if $dw < $up;
    }
}


#  Transfer all forward strand segments to reverse strand, in place.
sub allToRev{
    my $self  =  $_[0];

    for my $i (0..$$self[SIZ]-1){
        my( $up, $dw )  =   unpack 'll',  substr $$self[VEC], 8*$i, 8;
        substr  $$self[VEC], 8*$i, 8   =>  pack 'll', $dw, $up    if $up < $dw;
    }
}



# @_  -->  $self, LIST
# push LIST onto $$self[VEC].  LIST must be even sized.
sub push{
    my $self  =  shift;
    @_ % 2  and  die  "in push(), expected an even sized list of (up, dw) pairs, but got list of size ", 0+@_;

    $$self[VEC]  .=   pack 'l*', @_;
    $$self[SIZ]  +=   @_ / 2;
}


#  Add $addend to every coordinate.
sub add{
    my( $self, $addend )  =  @_;

    for my $i (0..$$self[SIZ]-1){
        my( $up, $dw)  =   unpack 'll',  substr $$self[VEC], 8*$i, 8;
        substr  $$self[VEC], 8*$i, 8   =>  pack 'll', $up+$addend, $dw+$addend;
    }
}


#  Set each coordinate equal to $minuend minus the coordinates original value.
sub subtractFrom{
    my( $self, $minuend, $start  )  =  @_;

    for my $i ($start..$$self[SIZ]-1){
        my( $up, $dw )  =   unpack 'll',  substr $$self[VEC], 8*$i, 8;
        substr  $$self[VEC], 8*$i, 8  =>  pack 'll', $minuend-$up, $minuend-$dw;
    }
}


#  Append $seqSeg to $self
sub append{
    my( $self, $seqSeg )  =  @_;

    $$self[VEC]  .=  $$seqSeg[VEC];
    $$self[SIZ]  +=  $$seqSeg[SIZ];
}


#  Return coordinates of ith seqSeg
sub get{
    my( $self, $idx )  =  @_;
    $idx < $$self[SIZ]   or   die  "$idx > ", $$self->fin();

    return   unpack 'll',  substr $$self[VEC], 8*$idx, 8;
}


sub dump{  say  join ' ',  unpack 'l*', ${_[0]}[VEC]  }

#  Return coordinates as flat list  ( up0, dw0, up1, dw1, ⋯  )
sub list{   unpack 'l*', ${_[0]}[VEC]   }


sub asString{
    my $self  =  $_[0];

    my $retVal  =  "$$self[NAM]:";

    for my $i (0..$$self[SIZ]-1){
        my( $up, $dw )  =   unpack 'll',  substr $$self[VEC], 8*$i, 8;
        $retVal  .=  "\t$up $dw";
    }

    return $retVal;
}

1;
