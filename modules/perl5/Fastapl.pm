#!/usr/bin/perl -w
#  Author: Paul Horton
#  Copyright (C) 2011, 2017. Paul Horton, All rights reserved.
#  Creation Date: 2014.3.8
#  Last Modification: 20170802
#
#  Functions used by fastapl
#
package Fastapl;
use utf8;
use Exporter 'import';
use feature qw(say state);
use Carp;
use Getopt::Long qw(:config posix_default bundling permute);
use List::Util qw(min max);
use List::MoreUtils qw(all any minmax uniq);
use POSIX;
use FastaplSeqsegs;
use fastaplUtils;

@EXPORT  =  qw(
$seqIsCirc $seq @seq seq seg $seq1 $seq2 $head @head %head $id $comment $commentAtStartOfStream
fin fin1 fin2 len len1 len2
trim rc matches expandRegex widen
codonAA ORFs translate translations runs kmers addKmers asHash
pr ps Ps
printin downcase upcase blink bold color colorize iv ul black blue default gray green magenta red yellow
_dieWithUsage _looksLikeFastaFile
%opt_
$numSeqLines_ $maxSeqLineLenSeen_
);

use strict;
use warnings;


# ──────────  Declare vars shared between fastapl.pm and fastapl  ──────────
our %opt_;  #To hold command line options


# ──────────  Initialize Command Line Options  ──────────

#  Scripts passed directly from command line.
$opt_{mainScript}  =  '';  #Record script + (optionally) BEGIN and/or END blocks.
$opt_{begScript}   =  '';  #Begin script, alternative way to give BEGIN block.
$opt_{endScript}   =  '';  #End script, alternative way to give  END  block.

$opt_{modules}  =  undef;  #Modules passed by user for use in script.  Array ref.

$opt_{countP}     =  0;  #If true, output number of records for which user script returns true.
$opt_{grepP}      =  0;  #If true, print record when user record script returns true.
$opt_{grepQuietP} =  0;  #Quiet version of grepMode.  If matching record found returns 0, otherwise 1.

$opt_{generateStandaloneP}  =  0;  #If true, output program with scripts embedded.
$opt_{outputFilename}       = '';  #Output filename.  If omitted output goes to STDOUT.


$opt_{dnaP}             =  0;  #the user has declared input sequences are DNA?
my  $masking_flag_      =  0;  #When true, functions such as ORFs treat lower case chars in $seq as masked.
$opt_{sortP}            =  0;  #Sort records by user given script?
$opt_{printEachRecDS_P} =  0;  #Print each record as double stranded?
$opt_{printEachRecSS_P} =  0;  #Print each record as single strand?

#  Mutually exclusive ways for user to stipulate output sequence line length.
$opt_{printSeqOn1LineP} =  0;  #Output sequences on one line?
$opt_{seqLineLen}       =  0;  #Length to use when printing sequences.  0 value means guess sequence length from input.  Ignored when $opt_{printSeqOn1LineP} given.

$opt_{rulerInlineP}     =  0;  #Display seqs with inline ruler?
$opt_{rulerTrackP}      =  0;  #Display seqs with separate ruler track?
$opt_{seqIsCirc}        =  0;  #Expected input sequences are circular?
$opt_{fieldSeparator}  =  "\t";  #Field separator used to compute @head.
$opt_{validSeqChars}   =  '-a-zA-Z';  #Chars possible in sequence, others are ignored.
our $printORFs_arg  =  undef;  #When defined, indicates length of ORFs to print.

$opt_{htmlP}      =  0;  #Should markup be done in html?
$opt_{htmlTags}   = '';  #If so, with what open/close tags.

our $seqIsCirc;  #Sequence is from a circular molecule?  Inialized by command line, but okay to modify in record script.


my @ARGVSpec = (
    'e|script=s'              =>  \$opt_{mainScript},
    'b|begin=s'               =>  \$opt_{begScript},
    'O|circular'              =>  sub{  $opt_{seqIsCirc}= $seqIsCirc= 1  },
    'c|count'                 =>  \$opt_{countP},
    'C|colorize'              =>  \$opt_{colorizeP},
    'd|dna'                   =>  \$opt_{dnaP},
#   'e|script=s'              =>  \$mainScript_arg_,
    'z|end=s'                 =>  \$opt_{endScript},
    'F|field-separator=s'     =>  \$opt_{fieldSeparator},
    'g|grep'                  =>  \$opt_{grepP},
    'H|html'                  =>  \$opt_{htmlP},
    'html-tags=s'             =>  sub{  @opt_{qw(htmlP htmlTags)} = (1,$_[1])  },
    'l|line-length=i'         =>  \$opt_{seqLineLen},
    'm|masking'               =>  \$masking_flag_,
    'M|module=s@'             =>  \$opt_{modules},
    'n|generate-standalone'   =>  \$opt_{generateStandaloneP},
#       'h|help      Reserved for possible future addition.
    'o|output-file=s'         =>  \$opt_{outputFilename},
    'ORF:s'                   =>  \$printORFs_arg,
    'p|print-each-record'     =>  \$opt_{printEachRecSS_P},
    'P|print-each-record-ds'  =>  \$opt_{printEachRecDS_P},
    '1|print-seq-on-1-line'   =>  \$opt_{printSeqOn1LineP},
    'q|grep-quiet'            =>  \$opt_{quietP},
    'Q|valid-seq-chars=s'     =>  \$opt_{validSeqChars},
    'r|ruler-inline'          =>  \$opt_{rulerInlineP},
    'R|ruler-track'           =>  \$opt_{rulerTrackP},
#        w|warn      Reserved for possible future addition.
    's|sort'                  =>  \$opt_{sortP},
#        v|verbose   Reserved for possible future addition.
    );


#  Return literal to include in standalone, suitable for passing to parseCommandLineOptions.
sub optionsForStandalone{
    my $retVal;

    $retVal .=  '-C '   if $opt_{colorizeP};
    $retVal .=  '-d '   if $opt_{dnaP};
    $retVal .=  "-F $opt_{fieldSeparator} "    if $opt_{fieldSeparator} ne "\t";
    $retVal .=  '-O '   if $opt_{seqIsCirc};
    $retVal .=  '-p '   if $opt_{printEachRecSS_P};
    $retVal .=  '-P '   if $opt_{printEachRecDS_P};
    $retVal .=  '-1 '   if $opt_{printSeqOn1LineP};
    $retVal .=  '-r '   if $opt_{rulerInlineP};
    $retVal .=  '-R '   if $opt_{rulerTrackP};
    $retVal .=  "-l $opt_{seqLineLen} "    if $opt_{seqLineLen};

    $retVal .=  '-H '                            if $opt_{htmlP};
    $retVal .=  "--html-tags $opt_{htmlTags} "   if $opt_{htmlTags};

    $retVal .=  '--ORF '             if         $printORFs_arg;
    $retVal .=  "$printORFs_arg "    if defined $printORFs_arg;

    $retVal   or   return '';
    return  "[qw( $retVal)]";
}



my $htmlP;  #speed optimization copy of $opt_{htmlP}
my %color;  #$color{RESIDUE} --> COLOR


sub parseCommandLineOptions{
    my $argv =  shift;

    my $getOptionsRetval  =  Getopt::Long::GetOptionsFromArray( $argv, @ARGVSpec, @_ );
    $getOptionsRetval  or  return 0;

    $htmlP = $opt_{htmlP};
    init_htmlOpenClozStrs()   if $htmlP;

    %color  =  (
        $opt_{dnaP}?
        qw(a blue c yellow g green t red n gray)
        :
        qw(h red k red_bold r red_bold  l magenta i magenta v magenta m magenta_bold
       n blue q blue b blue z blue d blue_bold e blue_bold  c cyan p cyan_bold
       f green y green w green_bold  a yellow s yellow t yellow g yellow_bold  x gray)
        );

    return $getOptionsRetval;
}


#  Print $message to STDERR.
#  if STDOUT is a terminal:  print usage summary to it;  otherwise print $message to it.
sub _dieWithUsage{
    my $message  =  "Command line parsing error; $_[0]";
    say STDERR  $message;
    say(    -t STDOUT?   "\nFor more info try:  fastapl --usage"   :   "fastapl: $message"    );
    exit 64;

}#  _dieWithUsage( $message )

our $seq;           #Sequence as scalar.
our @seq;           #Sequence as array, one element per residue.
our $head;          #Head line (line starting with '>').
our @head;          #Head line after id, split by $opt_{fieldSeparator} into array.
our %head;          #Holds key,value pairs when $head looks like key=value
our $id;            #Record id (non-white-space string following the '>').
our $comment = '';  #Record comment lines as single string (newlines included).
our $commentAtStartOfStream = '';  #To hold comment lines at start of stream.  Emptied after one record printed by pr().

our %codonAA  =  qw(
aaa K aac N aag K aat N aau N aan X aca T acc T acg T act T acu T acn T aga R agc S agg R agt S agu S agn X ata I aua I atc I auc I atg M aug M att I auu I atn X aun X ana X anc X ang X ant X anu X ann X
caa Q cac H cag Q cat H cau H can X cca P ccc P ccg P cct P ccu P ccn P cga R cgc R cgg R cgt R cgu R cgn R cta L cua L ctc L cuc L ctg L cug L ctt L cuu L ctn L cun L cna X cnc X cng X cnt X cnu X cnn X
gaa E gac D gag E gat D gau D gan X gca A gcc A gcg A gct A gcu A gcn A gga G ggc G ggg G ggt G ggu G ggn G gta V gua V gtc V guc V gtg V gug V gtt V guu V gtn V gun V gna X gnc X gng X gnt X gnu X gnn X
taa * tac Y tag * tat Y tau Y tan X tca S tcc S tcg S tct S tcu S tcn S tga * tgc C tgg W tgt C tgu C tgn X tta L tua L ttc F tuc F ttg L tug L ttt F tuu F ttn X       tna X tnc X tng X tnt X tnu X tnn X
uaa Y uac Y uag * uau Y uan X uca S ucc S ucg S ucu S ucn S uga * ugc C ugg W ugu C ugn X uua L uuc F uug L uuu F uun X una X unc X ung X unu X unu X unn X
naa X nac X nag X nat X nau X nan X nca X ncc X ncg X nct X ncu X ncn X nga X ngc X ngg X ngt X ngu X ngn X nta X nua X ntc X nuc X ntg X nug X ntt X nuu X ntn X nun X nna X nnc X nng X nnt X nnu X nnn X
);


our( $seq1, $seq2 );  #For sorting

my $idCopy;  #ID of current record as read in or as it was in the most recent call of pr().

our $numSeqLines_  =  0;  #Number of sequence lines.

our $maxSeqLineLenSeen_ = 0;


#  vecs to hold values in [0,7].  O no AA, 1-3 AA frame.  +4 iff initial M of ORF.
my( $ORFtrackFrameFwd,  $ORFtrackFrameRev)  =  ('','');

my @ORFtrackAny  =  (0,0);   #Is ORFtrackFrameFwd or ORFtrackFrameRev used at all?
my $ORF_colorize;   #Color ORF track?



sub seq($);     #  seq(i) returns ith character in $seq.
sub fin();      #  Returns length($seq) - 1
sub rc;         #  Reverse complement in place.
sub len();      #  Returns length($seq)
sub codonAA($); #  Returns one-letter amino acid code for residue corresponding to codon $.
sub ORFs;       #  List of ORFs.

sub ps;    #  ps( $beg, $end ) prints ($beg, $end) substring of $seq.
sub Ps;    #  Ps( $beg, $end ) prints ($beg, $end) substring of $seq as double stranded DNA.


# ──────────  Internal Functions  ──────────

sub processHead(){
    chomp $head;
    ($idCopy)  =  ($id)  =  ( $head =~ /^ \s* (\S*) /x );
    @head  =   split $opt_{fieldSeparator}, $head;
    %head  =   map {/([^=]+)=([^=]+)/} @head;
}

#  Final index of $seq.
sub fin(){   length($seq) - 1   }

#  Sequence length.
sub len(){    length $seq   }


sub fin1(){   length($seq1) - 1   }
sub fin2(){   length($seq2) - 1   }
sub len1(){   length $seq1        }
sub len2(){   length $seq2        }


#  Return $ith base in sequence.
sub seq($){
    my $i  =  $_[0];
    die  "In function seq(), index ($i) out of range"    if  $i > fin;

    return  substr $seq, $i, 1;
}



#  Return sequence segment $up,$dw as string
sub seg{
    my( $u, $d )  =  @_;

    @_ == 2   or   die 'seq expected two args but got', join ', '=>@_;

    my $L = len;

    if(  $d < $u  ){
        if(  0 <= $d  ){                          #   0 < $d < $u
            ;    $L >= $u  or  goto DIE_RANGE;    return rc   substr  $seq,    $d,  $u-$d
        }else{#    $d < 0  $d < $u
            if(  $u <= 0  ){                      #  $d < $u ≦  0
                -$L <= $d  or  goto DIE_RANGE;    return rc   substr  $seq, $L+$d,  $u-$d
            }
            $seqIsCirc     or  goto DIE_WRAPD;    #  $d <  0 ≦ $u
            $L >= $u-$d    or  goto DIE_RANGE;    return rc   substr( $seq, $L+$d, -$d )  .  substr $seq, 0, $u
        }
    }else{#  $u ≦ $d
        if(  0 <= $u  ){                          #   0 ≦ $u ≦ $d
            ;    $L >= $d  or  goto DIE_RANGE;    return      substr  $seq,    $u,  $d-$u
        }else{
            if(  $d <= 0  ){                      #  $u ≦ $d ≦  0
                -$L <= $u  or  goto DIE_RANGE;    return      substr  $seq, $L+$u,  $d-$u
            }
            $seqIsCirc     or  goto DIE_WRAPD;    #  $u <  0 < $d
            $L >=  $d-$u   or  goto DIE_RANGE;    return      substr( $seq, $L+$u, -$u )  .  substr $seq, 0, $d
        }
    }

  DIE_RANGE:
    die  "In seg(); coords ($u,$d) are incompatible with \$seq length $L\n";

  DIE_WRAPD:
    die  "In seg(); coords ($u,$d) cross zero, but \$seq is not circular\n";

}#END: seg


#  Internal function for use when pattern matching on circular segments.
#  Like seg() but returns seg concatenated against itself.
sub segseg{
    my( $up, $dw )  =  @_;

    @_ == 2   or   die 'seq expected two args but got', join ', '=>@_;

    my $L = len;

    if(  $dw < $up  ){
        return rc  $seq.$seq   if $dw >= 0;
        return rc  $seq.$seq   if $up >= 0;
        return rc  substr $seq, $L+$dw  .  $seq  .  substr $seq, 0, $up;
    }else{#  $up ≦ $dw
        return     $seq.$seq    if $up >= 0;
        return     $seq.$seq    if $dw <= 0;
        return     substr $seq, $L+$up  .  $seq  .  substr $seq, 0, $dw;
    }

}#END segseq


#  Convert a ($up, $dw) seqSeg to a pair of ranges,
#  suitable for use in memory efficient for loops.
#
#  Returns list of one or two triples.
#
#    $rev  =  ($up > $dw)?;
#    @range = seqSeg_to_directedRanges( $up, $dw );
#    while(  my( $iDéb, $iFin, $printAsNeg )  =  splice @range, 0, 3  ){
#      for my $i ($iDéb..$iFin){
#        $i *= -1  if $rev;
#        # current char is substr $i, 0, 1
#        $pos  =   $i  -  len * $printAsNeg;


#  Differs from seqSeg_to_numRanges in the ORDER of the i's.
#  This function returns a range suitable for spelling out a substring
#  of the DNA on the stipulated strand and possibly crossing the end of $seq and
#  wrapping back into the beg of $seq in the case of circular DNA.
#
#  The third returned argument is always zero (FALSE) for linear $seq.
#  For circular sequences it is one (TRUE) for the negative part of wrapped seqSeg.
#
#  Six example cases in order of treatment in code below.
#
#                           0 1 2 3 4 5 6 7
#         → a b c d e f g h|a b c d e f g h
#           A B C D E F G H|A B C D E F G H ←
#  (+5,+2)  -7-6-5-4-3-2-1-0 1 [ 3 4←] 6 7 8 → -4..-2              = EDC
#  (-2,-5)  -7-6-[-4-3←]-1-0 1 2 3 4 5 6 7 8 → -5..-3              = FED
#  (+3,-2)  -7-6-5-4-3-[-1-0 1 2←] 4 5 6 7 8 → -2..-0 F, -7..-6 T  = CBAHG
#
#  (+2,+5)  -7-6-5-4-3-2-1 0 1 [→3 4 ] 6 7 8 →  2..+4              = def
#  (-4,-1)  -7-6-5-[→3-2-] 0 1 2 3 4 5 6 7 8 →  4..+6              = efg
#  (-2,+3)  -7-6-5-4-3-[→1 0 1 2 ] 4 5 6 7 8 →  6..+7 T,  0..+2 F  = ghabc
sub seqSeg_to_directedRanges{
    my( $u, $d, $msgPrefix )  =  @_;
    my $L  =  len();

    if(  $d < $u  ){
        if(  0 <= $d  ){                      #   0 < $d < $u
            +$L >= $u  or  goto DIE_RANGE;    return    1-$u,     -$d, 0
        }
        if(  0 >= $u  ){                      #  $d < $u ≦ 0
            -$L <= $d  or  goto DIE_RANGE;    return 1-$u-$L,  -$L-$d, 0
        }
        $seqIsCirc     or  goto DIE_WRAPD;    #  $d <  0 < $u
        $L >=  $u-$d   or  goto DIE_RANGE;    return    1-$u,       0, 0      ,1-$L, -$L-$d, 1
    }{#else $u ≦ $d
        if(  0 <= $u  ){                      #   0 ≦ $u ≦ $d
            +$L >= $d  or  goto DIE_RANGE;    return      $u,    $d-1, 0
        }
        if(  0 >= $d  ){                      #  $u ≦ $d ≦ 0
            -$L <= $u  or  goto DIE_RANGE;    return   $L+$u, $L+$d-1, 0
        }
        $seqIsCirc     or  goto DIE_WRAPD;    #  $u <  0 < $d
        $L >=  $d-$u   or  goto DIE_RANGE;    return   $L+$u,    $L-1, 1      ,0,      $d-1, 0
    }

  DIE_RANGE:
    die  "$msgPrefix seqseg ($u,$d) is incompatible with \$seq length $L\n";

  DIE_WRAPD:
    die  "$msgPrefix seqseg ($u,$d) crosses zero, but \$seq is not circular\n";
}#END: seqSeg_to_directedRanges



#  seqSeg_to_numRanges( POS1, POS2 )
#
#  Return range (sometimes two for circular sequences)
#  for seqSeg (POS1, POS2) or equivalently (POS2, POS1)
#  suitable for use in memory efficient for loops.
#
#    @range = seqSeg_to_numRanges( $up, $dw );
#    while(  ($a, $b)  =  splice @range, 0, 2  ){
#      for my $i ($a..$b){
#        do something with $i...
#
#  Differs from seqSeg_to_directedRanges in the ORDER of the i's.
#  This function gives the i's in ascending numerical order.
#
#                 0 1 2 3 4 5 6 7
#             g h|a b c d e f g h
#  (-2,+3)  -B-1 0 1 2 E 4 5 6 7 0  -->  0..2, 6..7  -->  abcgh
#  (+3,-2)  -E-1 0 1 2 B 4 5 6 7 0  -->  0..2, 6..7  -->  abcgh
#
#  DOES NOT DO RANGE CHECKING!
sub seqSeg_to_numRanges{
    my( $p, $G )  =  $_[0] < $_[1]?  @_[0,1]  :  @_[1,0];   # petite, Grand

    my $l  =  len();

    return     $p,    $G-1    if $p >  0;
    return  $l+$p, $l+$G-1    if $G <= 0;

    #  $p < 0 < $G  only valid for circular sequence
    return      0, $G-1,     $l+$p, $l-1
}#END: seqSeg_to_numRanges



#  Shorten sequence and quality to maximum length $len, setting it to
#  the forward strand seqSeg ($up, $dw).
#  Return length of sequence after trimming.
sub trim{
    my( $up, $dw )  =  ( @_, len );

    die  "trim called with up ($up) > dw ($dw)"   if $up > $dw;

    $up  +=  len    if  $up < 0;
    $dw  +=  len    if  $dw < 0;

    my $lenToKeep  =  $dw - $up + 1;
    $seq  =   substr $seq, $up, $lenToKeep;

    @seq  =  splice @seq, $up, $lenToKeep    if $::needSeqArray_;

    return $lenToKeep;
}#END:  trim( $up, $dw )



# Internal function for rc().
# Set $_[0] to its reverse complement.
sub _rc{
    defined $_[0]   or   confess  '_rc passed undefined arg';
    $_[0]  =  reverse $_[0];

    if(  length $_[0]  !=  $_[0] =~ tr/ACGTUNSWRYKMBDHVacgtunswrykmbdhv/TGCAANSWYRMKVHDBtgcaanswyrmkvhdb/  ){
        (my $offendingPart)  =  $_[0] =~ /([^ACGTUNSWRYKMBDHVacgtunswrykmbdhv]+)/;
        die  "While attempting to reverse complement, encountered chars [$offendingPart]\n.";
    }
}


#  Reverse Complement sequence.
#  Modifies in place when called in void context, otherwise returns reverse complement(s).
#  When called with no args, set $seq to its reverse complement.
#  Otherwise it should be called with a list of sequences.
#  The conversion covers ambiguous letters and partially covers RNA by mapping 'U' --> 'A', but 'A' is always mapped to 'T'.
sub rc{
    if(  defined wantarray  ){
        if(  wantarray  ){
            if(  @_  ){
                my @retVal  =  @_;
                _rc $_   for @retVal;
                return @retVal;
            }
        }else{# scalar context
            my $retVal  =  @_? $_[0] : $seq;
            _rc $retVal;
            return $retVal;
        }
    }

    # ELSE  ─────  Void context  ─────
    if(  @_  ){
        _rc( ref $_? $$_ : $_ )  for @_;
    }
    else{  #No arguments passed to rc().
        die  'rc() without arguments called when sorting'   if $opt_{sortP};

        _rc $seq;

        #  Synch @seq up with $seq.
        @seq  =  split //, $seq    if $::needSeqArray_;
    }
}


{# Scope for expandRegex
#  When $opt_{dnaP} given, expandRegex( $regex ) returns an expanded version of $regex
#  Examples: 'SATYW' --> '[GC]AT[CT][AT]'
#            'a<MA>' --> 'aA[TU]GGC.'
#  However does not substitute inside of existing brackets.  e.g. in 'A[TW]G' would remain 'A[TW]G'

    my %nt_IUPAC = qw(R [AG]  Y [CT]  S [CG] W [AT]  K [GT]  M [AC]  B [CGT]  D [AGT]  H [ACT]  V [ACG]  N [ACGT]  r [ag]  y [ct]  s [cg] w [at]  k [gt]  m [ac]  b [cgt]  d [agt]  h [act]  v [acg]  n [acgt]);

    my %AAcodonRegex = qw(A GC.  B [AG]A[CTU]  C [TU]G[CTU]  D GA[CTU]  E GA[AG]  F [TU][TU][CTU]  G GG.  H CA[CTU]  I A[TU][ACTU]  K AA[AG]  L (?:C[TU].|[TU][TU][AG])  M A[TU]G  N AA[CTU]  P CC.  Q CA[AG]  R (?:CG.|AG[AG])  S (?:[TU]C.|AG[CTU])  T AC.  V G[TU].  W [TU]GG  Y [TU]A[CTU]  Z [CG]A[AG]);

    sub expandRegex{
        $opt_{dnaP}   or   return  $_[0];

        my $origRegex  =  $_[0];
        my @regexPart  =    split   /(  < [^>]* >  )/x,   $origRegex;

        my $newRegexStr;
        for my $regexPart (@regexPart){
            if(   (my $pat)  =  $regexPart =~  /< ([^>]*) >/x   ){
                die  "Not able to process codon substitution into string '<$pat>' containing brackets"
                    if  $pat =~   /  \[  |  ]  /x;
                my $reversed  =  reverse $pat;
                while(  my $c = chop $reversed  ){
                    $newRegexStr .=  $AAcodonRegex{$c} // $c;
                }
            }
            else{
                my $netCount  =  0;
                my $reversed  =  reverse $regexPart;
                while(  my $c = chop $reversed  ){
                    $netCount++   if $c eq '[';
                    $netCount--   if $c eq ']';
                    $newRegexStr .=  (!$netCount  and  $_ = $nt_IUPAC{$c})?  $_  :  $c;
                }
                $netCount == 0   or   die  "Unbalanced [] brackets in regex '$origRegex'.";
            }
        }

        return $newRegexStr;
    }
}#END: expandRegex scope.


my %expandedRegex;   #Cache for expanded regular expressions.

#  Returns FastaplSeqsegs object containing list of segments in $seq matching $regex
#  Optional trailing control chars ［%_=］ in $regex can be used to contol what type of matches are included.
#  Empty FastaplSeqsegs object is returned if $regex is empty.
sub matches{
    @_   or   die  'matches() expected one args but got none';

    my $rawRegex  =  shift   or   die  'matches() expected at least one regex arg but got none';
    my $matches  =  new FastaplSeqsegs( $rawRegex );

    my $regex  =  $expandedRegex{$rawRegex} //= expandRegex $rawRegex;

    my $overlapOK       = 0;  #Default: do not return overlapping matches
    my $lookOnStrand    = 0;  #Default: just match on same strand
    my $recordAllOnFwd  = 0;  #Default: record reverse strand matches as on the reverse strand.
    while(  substr( $regex, -1, 1 ) =~ /[_=%1]/  ){
        my $c  =  chop $regex;
        $overlapOK      = 1   if $c eq '%';  #Include overlapping matches
        $lookOnStrand   = 1   if $c eq '_';  #Just match on opposite strand
        $lookOnStrand   = 2   if $c eq '=';  #Match on both strands
        $recordAllOnFwd = 1   if $c eq '1';  #Record all matches as if on the forward strand.
    }


    @_   or   @_ = (0,len);
    my( $seqSegs, @coord1_ )  =  @_;

    ref( $seqSegs )  =~ /FastaplSeqsegs/
        or  $seqSegs  =  FastaplSeqsegs->new( 'to search in', $seqSegs, @coord1_ );

    my $L = len;

    my $nextPair  =  $seqSegs->itr();
    while(   my( $searchUp, $searchDw )  =  $nextPair->()   ){
        my $revP  =  ($searchDw < $searchUp);
        if(  !$seqIsCirc   or   $L != abs $searchUp - $searchDw   ){
            my $seg  =  seg $searchUp, $searchDw;
            if(  $lookOnStrand != 1  ){
                while(  $seg =~ /$regex/g  ){
                    if( $revP ){  $matches->push(  $searchUp - $-[0],  $searchUp - $+[0]  )}
                    else       {  $matches->push(  $searchUp + $-[0],  $searchUp + $+[0]  )}
                    pos $seg  =  $-[0] + 1    if $overlapOK;
                }
            }
            if(  $lookOnStrand != 0  ){
                rc $seg;
                while(  $seg =~ /$regex/g  ){
                    if( $revP ){  $matches->push(  $searchUp + $-[0],  $searchUp + $+[0]  )}
                    else       {  $matches->push(  $searchUp - $-[0],  $searchUp - $+[0]  )}
                    pos $seg  =  $-[0] + 1    if $overlapOK;
                }
            }
        }
        else{#  Current seg is circular.
            my $segseg  =  segseg( $searchUp, $searchDw );

            if(  $lookOnStrand != 1  ){   #Search on same strand as seg.
                while(  $segseg =~ /$regex/g  ){
                    my( $matchUp, $matchDw )  =  ( $-[0], $+[0] );
                    last   if  $matchUp >= $L;
                    if(   $matchDw - $matchUp  >  $L   ){
                        next  if  substr( $segseg, $matchUp, $L )  !~  /\A$regex/;
                        $matchDw  =  $matchUp + $+[0];
                    }
                    if(  $revP  ){
                        my $offset  =  $searchUp;
                        $offset += $L   if $offset - $matchDw < -$L;
                        $matches->push(  $offset - $matchUp, $offset - $matchDw  );
                    }else{
                        my $offset  =  $searchUp;
                        $offset -= $L   if $offset + $matchDw > +$L;
                        $matches->push(  $offset + $matchUp, $offset + $matchDw  );
                    }
                    pos $segseg  =  $matchUp + 1   if $overlapOK;
                }
            }

            if(  $lookOnStrand != 0  ){   #Search on opposite strand of seg.
                rc $segseg;
                while(  $segseg =~ /$regex/g  ){
                    my( $matchUp, $matchDw )  =  ( $-[0], $+[0] );
                    last   if  $matchUp >= $L;
                    if(   $matchDw - $matchUp  >  $L   ){
                        next  if  substr( $segseg, $matchUp, $L )  !~  /\A$regex/;
                        $matchDw  =  $matchUp + $+[0];
                    }
                    my $offset  =  $searchUp;
                    if(  $revP  ){
                        $offset += $L   if $offset - $matchDw < -$L;
                        $matches->push(  $offset - $matchDw, $offset - $matchUp   );
                    }
                    else{
                        $offset += $L   if $offset - $matchDw < -$L;
                        $matches->push(  $offset - $matchUp, $offset - $matchDw  );
                    }
                    pos $segseg  =  $matchUp + 1   if $overlapOK;
                }
            }
        }#END if/else seg cicular
    }#END while ($searchUp,$searchDw)

    $matches->allToFwd()   if $recordAllOnFwd;
    return $matches;
}#END matches



#  ━━━━━━━━━━━━━━━━━━━━━━━━━  BEGIN translation related functions    ━━━━━━━━━━━━━━━━━━━━━━━━━
#  ORFspec form is [a][c][NUM][m][NUM]
sub setORFtrack{
    my $ORFspec  =  $printORFs_arg;

    my $needAcceptor  =  $ORFspec =~  /a/i;
    $ORF_colorize     =  $ORFspec =~ s/c//;

    my @minWeight = (100,80);  #Minimum weight for non-M, and M starting ORFs respectively.

    @minWeight  =  defined $2?  ($1, $2)  :  ($1, $1)    if $ORFspec =~  /a?  (\d+)  (?: m (\d+) )?/x;

    _setORFtrack( 0, $needAcceptor, @minWeight );

    rc $seq;

    _setORFtrack( 1, $needAcceptor, @minWeight );

    #  Adjust coordinates in ORFtrackFrameRev to count forward, instead of reverse.
    for my $i (0..fin/2){
        my $tmp  =  vec( $ORFtrackFrameRev, $i, 4 );
        vec( $ORFtrackFrameRev, $i, 4 )  =  vec( $ORFtrackFrameRev, fin-$i, 4 );
        vec( $ORFtrackFrameRev, fin-$i, 4 )  =  $tmp;
    }

    rc $seq;
}#END: setORFtrack


#  @minWeight has two weights: one for ORFs starting with M (MORFS) and the second one for other ORFs.
#  These balance the choice between a long ORF not starting M and its longest MORF suffix.
#  If a maximal ORF has a MORF suffix, the MORF is chosen unless the weight of ORF minus its MORF suffix
#  exceeds the threshold $minWeight[0]
sub _setORFtrack{
    my( $strand, $needAcceptor, @minWeight )  =  @_;
    my $ORFtrackFrameRf  =  $strand? \$ORFtrackFrameRev : \$ORFtrackFrameFwd;

    my $begCodonMid;   #Index into $seq of the middle of the first   codon in the current ORF.
    my $begModonMid;   #Index into $seq of the middle of the first M codon in the current ORF.
    my $curCodonMid;   #Middle of current codon to add to the current ORF.
    my @ORFcurWeight;  #The number of non-masked codons in current ORF.  @curWeight[0,1] correspond to weight before and from first 'M' onward, respectively.
    my $aa;            #To hold amino acid of the most recently examined codon.
    my $Mseen;         #1 iff an M has already been seen in current ORF.

    my $len  =  length $seq;
    my $acceptorRgx  =   $masking_flag_?  qr{[CT]AG}  :  qr{[ct]ag}i;
  FRAME:
    for my $frame (+1,+2,+3){
        @ORFcurWeight  =  (0,0);
        $begCodonMid = $begModonMid  =  undef;
        my $inORF = $Mseen =  0;
        my $curCodon;
        my $numNtsWrapped  =  undef;
        for(   $curCodonMid  =  $frame;
               $curCodonMid - 1  <  $len;
               $curCodonMid += 3
            ){

            $numNtsWrapped   =   2  +  $curCodonMid -  $len;
            if(  $numNtsWrapped >= 0  ){
                next   unless  $seqIsCirc;
                next   if  defined $begCodonMid  &&  $begCodonMid > $numNtsWrapped+1;
                $curCodon   =   substr(  $seq,  $numNtsWrapped - 3  )   .   substr(  $seq, 0, $numNtsWrapped  );
            }else{
                $curCodon  =  substr( $seq, $curCodonMid-1, 3 );
            }
            $aa   =   codonAA $curCodon;


            # ──────────  Start ORF if appropriate  ──────────
            if(  $aa eq 'M'  ){
                $inORF = 1;
                $Mseen = 1;
                $begCodonMid  //=  $curCodonMid;
                $begModonMid  //=  $curCodonMid;
            }
            else{
                if(  !$inORF  ){
                    if(  $aa ne '*'  or  $minWeight[0] == 0  ){
                        if(  $needAcceptor  ){
                            my $flank;
                            my $leftFlankSize  =  min( $curCodonMid-1, 5 );
                            if(  $seqIsCirc  &&  $leftFlankSize < 5  ){
                                $flank   =   substr( $seq, $leftFlankSize - 5 )  .  substr( $seq, 0, 5 - $leftFlankSize );
                            }else{
                                $flank   =   substr $seq, $curCodonMid-1-$leftFlankSize, $leftFlankSize;
                            }
                            if(  $flank =~ m/$acceptorRgx/  ){
                                $inORF  =  1;
                                $begCodonMid  //=  $curCodonMid;
                            }
                        }
                        else{#  ORF can start anywhere.
                            $inORF = 1;
                            $begCodonMid  //=  $curCodonMid;
                        }
                    }
                }
            }

            if(  $inORF  ){
                ++$ORFcurWeight[ $Mseen ]    if  $aa =~ /[A-WYZ]/;
                if(  $aa eq '*'  ){
                    $begCodonMid = undef          if  $minWeight[1] >  $ORFcurWeight[1]  &&  $minWeight[0] > $ORFcurWeight[0] + $ORFcurWeight[1];
                    $begCodonMid = $begModonMid   if  $minWeight[1] <= $ORFcurWeight[1]  &&  $minWeight[0] > $ORFcurWeight[0];
                    if(  defined $begCodonMid  ){
                        my $canExtendLeftward  =  0;
                        if(  $seqIsCirc  &&  $begCodonMid == $frame  &&  $curCodonMid+1 < $len  ){
                            my $prevCodon   =   substr( $seq, $frame-4 )  .  substr( $seq, 0, $frame-1 );
                            $canExtendLeftward  =  (codonAA $prevCodon ne '*');
                        }
                        unless(  $canExtendLeftward  ){
                            $ORFtrackAny[$strand]  =  1;
                            for(  my $i = $begCodonMid;   $i <= $curCodonMid;  $i += 3  ){
                                vec( $$ORFtrackFrameRf, ($i<$len? $i:0) , 4 )  =  $frame;
                            }
                            vec( $$ORFtrackFrameRf, 0, 4 )  =  $frame    if  $numNtsWrapped == 2;
                            vec( $$ORFtrackFrameRf, ($begCodonMid<$len? $begCodonMid:0), 4 )  +=  4;
                        }
                    }
                    @ORFcurWeight  =  (0,0);
                    $begCodonMid = $begModonMid  =  undef;
                    $inORF = $Mseen =  0;
                }

            }
        }#END:  for $seq starting positions in current frame.
        my $endCodonMid  =  $curCodonMid;

        if(  $inORF  &&  $seqIsCirc  ){
            #──────────  Finish off wrapped ORF  ──────────
            for(   $curCodonMid  =  $numNtsWrapped + 1;
                   $curCodonMid  <  $begCodonMid - 2;
                   $curCodonMid += 3
                ){

                $aa   =   codonAA  substr $seq,$curCodonMid-1,3;

                last  if  $aa eq '*';
                ++$ORFcurWeight[ $Mseen ]   if  $aa =~/[A-WYZ]/;

            }#END:  for $curCodonMid < $begCodonMid

            $begCodonMid = undef          if  $ORFcurWeight[1] <  $minWeight[1]  &&  $minWeight[0] > $ORFcurWeight[0] + $ORFcurWeight[1];
            $begCodonMid = $begModonMid   if  $ORFcurWeight[1] >= $minWeight[1]  &&  $minWeight[0] > $ORFcurWeight[0];
            if(  defined $begCodonMid  ){
                $ORFtrackAny[$strand]  =  1;
                for(  my $i = $begCodonMid;  $i <= $endCodonMid;  $i += 3  ){
                    vec( $$ORFtrackFrameRf, ($i<$len? $i:0), 4 )  =  $frame;
                }
                vec( $$ORFtrackFrameRf, ($begCodonMid<$len? $begCodonMid:0), 4 )  +=  4;
                for(   my $i = $numNtsWrapped + 1;  $i <= $curCodonMid;  $i += 3   ){
                    vec( $$ORFtrackFrameRf, $i, 4 )  =  $frame;
                }
            }
        }#END:  If circular
        else{#  end of sequence in non-circular mode.
            if(  $inORF  ){
                $begCodonMid = undef          if  $ORFcurWeight[1] <  $minWeight[1]  &&  $minWeight[0] > $ORFcurWeight[0] + $ORFcurWeight[1];
                $begCodonMid = $begModonMid   if  $ORFcurWeight[1] >= $minWeight[1]  &&  $minWeight[0] > $ORFcurWeight[0];
                if(  defined $begCodonMid  ){
                    $ORFtrackAny[$strand]  =  1;
                    for(  my $i = $begCodonMid;  $i < $len;  $i += 3  ){
                        vec( $$ORFtrackFrameRf, $i, 4 )  =  $frame;
                    }
                    vec( $$ORFtrackFrameRf, $begCodonMid, 4 )  +=  4;
                }
            }
        }
    }#Next frame
}#END: _setORFtrack


#
#  ORFs(  seqRef, ORFspec  )
#
#  seqRef defaults to \$seq
#
#  ORFspec contains an integer minWeight and an optional M or A.
#  An M means only consider ORFs starting with M.
#  An A means consider ORFs preceded by an AG splicing acceptor or starting with M
#  minWeight is the minimum length of an ORF, except that
#  when the masking flag is given amino acids from codons
#  overlapping with masked nucleotides are not included in the ORF length
#
#  Return list of ORFs in nucleotide sequence $nSeq matching the ORFspec criteria
#  Each item of this list is a quadruple [AASequence, begin, end, frame].
#  Where,
#    ・ [begin, end) are the nucleotide sequence coordinates of the ORF
#    ・ frame is the reading frame relative to the start of the $$seqRef for forward frames
#       and the end of $$seqRef for reverse strand frames.
#    ・ AASequence holds the ORF as an amino acid sequence.
#
#  Three cases of ORF termination.
#    ・ A stop codon or 'nnn' codon is encountered:
#         - the stop codon is included in [begin, end)
#         - amino_acid_sequence ends in '*'
#    ・ The ORF runs off the end of the nucleotide sequence:
#         - [begin, end) includes no stop codon.
#         - '\' is appended to AASequence, but not reflected in [begin, end)
#    ・ (in circularMode) the ORF wraps around to where it started
#         - [begin, end) includes the part of the ORF such that no codons overlap at all
#         - '@' is appended to AASequence, but not reflected in [begin, end)
#
#  Can be called with 0 to 2 arguments.
#      The nucleotide sequence reference defaults to $seq.
#      The minimum weight defaults to one.  Weight is like length, but masked codons are not added in.
#      setting the minimum weight to zero is a special case, in which stop codons are read through
#      and only running of the sequence stops an ORF.  This was done for convenience when calling
#      ORFs from setORFtrack.
#
#  Usage:
#    ORFs()             -->  ORFs( \$seq, 'm50' )
#    ORFs( $r )         -->  ORFs( \$r,   'm50' )
#    ORFs( $r, '' )     -->  ORFs( \$r,   'm50' )
#    ORFs( 'a70'  )     -->  ORFs( \$seq, 'a70' )
#    ORFs( $r, 'a70' )  -->  ORFs( \$r,   'a70' )
#    ORFs( 'a70', $r )  ERROR.
sub ORFs{
    my $err  =  'Error in sub ORFs;';
    @_ < 3   or   die  "$err expected at most two args, but got (@_)";
    my( $ntSeqRf, $ORFspec )  =  (\$seq, 'm50');

    if(  @_ == 1  ){
        if(  $_[0]  =~  /\d+/   ){   $ORFspec = $_[0]   }
        else                     {   $ntSeqRf = $_[0]   }
    }
    if(  @_ == 2  ){
        $ntSeqRf  =  $_[0];
        $ORFspec  =  $_[1]   if length $_[0];
    }


    $ORFspec  =~  /^[mM]?[aA]?\d+$/   or   die   "$err invalid ORF specifier '$ORFspec', should match [ma]?\\d";

    my $needM         =  $ORFspec =~ /m/i;
    my $needAcceptor  =  $ORFspec =~ /a/i;
    $needM && $needAcceptor   and  die  "$err ORF specificier '$ORFspec' stipulated needing both an initial methionine and an AG splicing acceptor";

    my( $minWeight )  =  $ORFspec =~ /(\d+)/;
    $minWeight  //=  50;

    my $len  =  length $$ntSeqRf;
    return ()   if  $len < 3  or  $len < $minWeight;

    my @forwardORF  =  _ORFs( $ntSeqRf, $minWeight, $needM, $needAcceptor );

    #  Temporarily reverse complement sequence,
    #  and compute forward strand ORFs on that.
    rc $$ntSeqRf;
    my @reverseORF  =  _ORFs( $ntSeqRf, $minWeight, $needM, $needAcceptor );

    #  Adjust coordinates and frame.
    for my $ORF (@reverseORF){
        @$ORF[0..2]  =  ($$ORF[0] < 0)
            ?  (    -$$ORF[0],     -$$ORF[1], -$$ORF[2])
            :  ($len-$$ORF[0], $len-$$ORF[1], -$$ORF[2]);
    }

    rc $$ntSeqRf;

    return( @forwardORF, @reverseORF );
}#END: ORFs().


#  Internal function doing the work of ORFrags,
#  which checks arguments and provides default values.
#  In the future might think of a mode which allows {AUG,CUG,UUG} as start codons.
sub _ORFs{
    my $ntSeqRf      =  $_[0];  #Reference to nucleotide sequence.
    my $minWeight    =  $_[1];  #Minimum number of unmasked codons needed to report an ORF.
    my $needM        =  $_[2];  #Starting Methionine codon required or not
    my $needAcceptor =  $_[3];  #'ag' splicing acceptor site (immediately upstream

    my @retVal  =  ();

    my $aaSeq;         #Current amino acid sequence.
    my $begCodonBeg;   #In $$ntSeqRf, Start of current ORF.
    my $curCodonBeg;   #In $$ntSeqRf, Current position for examination to extend the current ORF.
    my $ORFcurWeight;  #The number of non-masked codons in current ORF.
    my $aa;            #To hold amino acid of the most recently examined ORF.

    my $len  =  length $$ntSeqRf;
    my @ORF1idxInFrame;  #to hold the idx into retVal of the first ORF in each frame, which in circular case may need to be extended leftward as a wrapped ORF
    my $acceptorRgx  =  $masking_flag_?  qr/[CT]AG/  :  qr/[ct]ag/i;
  FRAME:
    for my $frame (+1,+2,+3){
        $aaSeq  =  '';
        $ORFcurWeight  =  0;
        my $inORF  =  !($needM || $needAcceptor);   #Currently in ORF?
        for(   $begCodonBeg = $curCodonBeg  =  $frame - 1;
               $curCodonBeg +2  <  $len;
               $curCodonBeg += 3
            ){
            $aa   =   codonAA  substr $$ntSeqRf,$curCodonBeg,3;

            # ──────────  Start ORF if appropriate  ──────────
            if(  !$inORF  ){
                if(  $aa eq 'M'  ){
                    $inORF  =  1;
                    $begCodonBeg  =  $curCodonBeg;
                }
                elsif(  $needAcceptor  ){
                    my $flank;
                    my $leftFlankSize  =  min( $curCodonBeg, 5 );
                    if(  $seqIsCirc  &&  $leftFlankSize < 5  ){
                        $flank   =   substr( $$ntSeqRf, $leftFlankSize - 5 )  .  substr( $$ntSeqRf, 0, 5 - $leftFlankSize );
                    }else{
                        $flank   =   substr $$ntSeqRf, $curCodonBeg-$leftFlankSize, $leftFlankSize;
                    }
                    if(  $flank =~ m/$acceptorRgx/  ){
                        $inORF  =  1;
                        $begCodonBeg  =  $curCodonBeg;
                    }
                }
            }

            if(  $inORF  ){
                $aaSeq  .=  $aa;
                $ORFcurWeight++    if $aa =~ /[A-Z]/;

                if(  $aa eq '*'  &&  $minWeight  ){
                    if(  $ORFcurWeight >= $minWeight  ){
                        push @retVal, [$begCodonBeg, $curCodonBeg+3, $frame, $aaSeq];
                        $ORF1idxInFrame[$frame]  //=  [$#retVal, $begCodonBeg];
                    }
                    $aaSeq  =  '';   $ORFcurWeight  =  0;
                    $begCodonBeg  =  $curCodonBeg + 3;
                    $inORF  =  !$needM;
                }
            }
        }#END:  for $ntSeqRf starting positions in current frame.

        if(   $inORF  &&  $seqIsCirc  &&  $len > 5   ){
            my $begCodonBeg_wrapped  =  $begCodonBeg - $len;   #Convert to wrapped segment coordinate.

            #Frame of the wrapped ORF as viewed from the sequence start.
            my $wrappedFrame  =  $frame - (len % 3);
            $wrappedFrame  +=  3 * ($wrappedFrame < 1);

            #──────────  Process wrapped codon  ──────────
            my $numNucsRght  =  $len - $curCodonBeg;
            my $numNucsLeft  =  0;

            if(  $numNucsRght  ){
                $numNucsLeft  =  3 - $numNucsRght;

                my $wrappedCodon   =   substr( $$ntSeqRf, -$numNucsRght )  .  substr( $$ntSeqRf, 0, $numNucsLeft );
                $aaSeq    .=    $aa  =  codonAA $wrappedCodon;
                $ORFcurWeight++    if $aa =~ /[A-Z]/;

                if(  $aa eq '*'  &&  $minWeight  ){
                    push @retVal, [$begCodonBeg_wrapped, $numNucsLeft, $frame, $aaSeq]   if  $ORFcurWeight >= $minWeight;
                    next FRAME;
                }
            }#END:  if(  $numNucsRght  )


            #──────────  Finish off wrapped ORF  ──────────
            for(   $curCodonBeg  =  $numNucsLeft;
                   $curCodonBeg  <  $begCodonBeg - 2;
                   $curCodonBeg += 3
                ){

                $aaSeq     .=     $aa   =   codonAA  substr $$ntSeqRf,$curCodonBeg,3;
                $ORFcurWeight++    if $aa =~ /[A-Z]/;

                if(  $aa eq '*'  &&  $minWeight  ){
                    if(  $ORFcurWeight >= $minWeight  ){
                        #  If necessary, replace ORF1 in retVal for complete wrapped ORF.
                        if(  defined $ORF1idxInFrame[$wrappedFrame]  &&
                             $ORF1idxInFrame[$wrappedFrame]->[1] < $curCodonBeg  ){
                            $retVal[$ORF1idxInFrame[$wrappedFrame]->[0]]  =  [$begCodonBeg_wrapped, $curCodonBeg+3, $frame, $aaSeq];
                        }
                        else{
                            push @retVal, [$begCodonBeg_wrapped, $curCodonBeg, $frame, $aaSeq];
                        }
                    }
                    next FRAME;
                }
            }#END:  for $curCodonBeg < $begCodonBeg

            if(  $aaSeq  ){
                $aaSeq .= '@';
                if(  $ORFcurWeight >= $minWeight  ){
                    #  If necessary, replace ORF1 in retVal for complete wrapped ORF.
                    if(  defined $ORF1idxInFrame[$wrappedFrame]  &&
                        $ORF1idxInFrame[$wrappedFrame]->[1] < $curCodonBeg  ){
                        $retVal[$ORF1idxInFrame[$wrappedFrame]->[0]]  =  [$begCodonBeg_wrapped, $curCodonBeg, $frame, $aaSeq];
                    }
                    else{
                        push @retVal, [$begCodonBeg_wrapped, $curCodonBeg, $frame, $aaSeq];
                    }
                }
            }
        }#END:  If circular
        else{
            if(   $aaSeq  &&  $ORFcurWeight >= $minWeight   ){
                $aaSeq  .=  '\\';
                push @retVal, [$begCodonBeg, $curCodonBeg, $frame, $aaSeq];
            }
        }
    }#END: _ORFfrags(), for my $frame ('+1','+2','+3')


    return @retVal;

}#END: _ORFs( $ntSeqRf, $minWeight )






#  Return one-letter amino acid code for residue corresponding to $codon.
sub codonAA($){
    my $codon  =  lc $_[0];
    my $retVal  =  $codonAA{$codon};
    defined $retVal   or   confess  "no amino acid defined for codon '$_[0]'";
    $retVal =~ tr/[A-Z]/[a-z]/    if  $masking_flag_  &&  $codon ne $_[0];
    return $retVal;
}


#  Returns translated version of (nucleotide) sequence.
#  $frame defaults to +1, sequence reference $$s defaults to $seq.
#  frame, sequence or both can be omitted from the argument list, but when both are given the order
#  must be (sequence, frame).
#  Usage:
#    translate();          #  Returns translation of $seq in reading frame +1
#    translate( -2 );      #  Returns translation of $seq in reading frame -2
#    translate( $s );      #  Returns translation of $s   in reading frame +1
#    translate( $s, -2 );  #  Returns translation of $s   in reading frame +1
#    translate( -2, $s );  #  ERROR.
sub translate{
    my $err  =  'Error in sub translate;';

    @_ < 3   or   die  "$err expected at most two args, but got (@_)";
    if   (  @_ == 0  ){
        return translateFrame( \$seq, 1 );
    }
    elsif(  @_ == 1  ){
        if(  $_[0]  =~  /^[-+]?[0123]$/  ){   #  Looks like frame specifier.
            my $frame  =  $_[0];
            die  "$err zero is not a valid frame, use {1,2,3,-1,-2,-3}"   if  $frame == 0;
            return  translateFrame(  \$seq,  $frame  );
        }
        else{   #  $_[0] should be a sequence.
            return  translateFrame(  \$_[0],  1  );
        }
    }
    else{   #2 arguments given.
        my ($s, $frame)  =  @_;
        $frame  =~  /^[-+]?[123]$/   or   die  "$err invalid frame specifier '$frame'";
        return  translateFrame(  \$s,  $frame  );
    }
} # END: translate().


#  Return a 6-element list holding the translations of nucleotide.
#  sequence $nucSeqRf (defaulting to $seq) in each of the 6 possible frames.
sub translations{
    die  "Error in sub translations; called with more than one argument (@_)"   if  @_ > 1;

    my $nucSeqRf  =   @_?  \$_[0]  :  \$seq;   #  Reference to nucleotide sequence.
    my @retVal;

    push @retVal, translateFrame( $nucSeqRf, $_)   for (1,2,3,-1,-2,-3);  #  Translate each frame.

    return @retVal;
}


#  Internal function, return translation of nucleotide sequence $$nucSeqRf in frame $frame.
#  $frame ∈ {+1,+2,+3,-1,-2,-3}
sub translateFrame{
    my( $nucSeqRf, $frame )  =  @_;

    my $codon;
    my $retVal  =  '';
    my $nucSeqLen  =  length $$nucSeqRf;

    rc $$nucSeqRf  if $frame < 0;

    my $i;   # index into $$nucSeqRf, value also used after loop exit.
    for(   $i = abs($frame) - 1;  $i < $nucSeqLen-2;  $i += 3   ){
        $retVal .=  codonAA( substr $$nucSeqRf, $i, 3 )
    }

    if(    $seqIsCirc    and    abs($frame) - 1  +  $nucSeqLen - $i  >  2    ){
        $codon    =    substr( $$nucSeqRf, $i )   .   substr( $$nucSeqRf, 0, abs($frame)-1 );
        $retVal .=  codonAA $codon;
    }

    rc $$nucSeqRf  if $frame < 0;

    return $retVal;
}#END: translateFrame()



#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  Markup code mostly below  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

my $termType_   = 'xterm';  #Type of terminal to assume when marking up output.

#Use like:  printin( 'green,ul', 'leprechaun' )
sub printin{
    $htmlP?  &fastaplUtils::printinHtml  :  &fastaplUtils::printinXterm;
}


sub  bold;     #  bold( COORDS ) marks up COORDS part of $seq as bold.
sub  iv;       #    iv( COORDS ) marks up COORDS part of $seq as inverse video.
sub  ul;       #    ul( COORDS ) marks up COORDS part of $seq as underline.

#  Functions to allow color words to be used like barewords
sub  black  {  'black'    }
sub  blue   {  'blue'     }
sub  default{  'default'  }
sub  gray   {  'gray'     }
sub  green  {  'green'    }
sub  magenta{  'magenta'  }
sub  red    {  'red'      }
sub  yellow {  'yellow'   }

#  The reason why the following markup vecs come paired with boolean "any" variables
#  is that I do not know of an effecient way to test if a vec has no non-zero elements.
#  It seems that tests like:  if $blFwd,  if length($blFwd),  can return false even when $blFwd is not all zero
#  Execution speed is the reason I did not package the {vec, any} pairs in a class.
my $blFwd  = '';   #  byte vector to hold forward strand boolean vec for blink.
my $blRev  = '';   #  byte vector to hold reverse strand boolean vec for blink.
my $blAny  =  0;   #  Is blink used at all?

my $boFwd  = '';   #  byte vector to hold forward strand boolean vec for bold.
my $boRev  = '';   #  byte vector to hold reverse strand boolean vec for bold.
my $boAny  =  0;   #  Is bold used at all?

my $ivFwd  = '';   #  byte vector to hold forward strand boolean vec for iv
my $ivRev  = '';   #  byte vector to hold reverse strand boolean vec for iv
my $ivAny  =  0;   #  Is iv used at all?

my $ulFwd  = '';   #  byte vector to hold forward strand boolean vec for ul.
my $ulRev  = '';   #  byte vector to hold reverse strand boolean vec for ul.
my $ulAny  =  0;   #  Is ul used at all?

my @blOpen  = (0,0);   #  $blOpen[STRAND] true iff blink currently open on STRAND.
my @boOpen  = (0,0);   #  $boOpen[STRAND] true iff bold  currently open on STRAND.
my @ivOpen  = (0,0);   #  $ivOpen[STRAND] true iff iv    currently open on STRAND.
my @ulOpen  = (0,0);   #  $ulOpen[STRAND] true iff ul    currently open on STRAND.

sub blink{
    $boAny  =  1;
    return setBoolVec( 'blink', \$blFwd, \$blRev, @_ );
}

sub bold{
    $boAny  =  1;
    return setBoolVec( 'bold',  \$boFwd, \$boRev, @_ );
}

sub iv{
    $ivAny  =  1;
    return setBoolVec( 'iv',    \$ivFwd, \$ivRev, @_ );
}

sub ul{
    $ulAny  =  1;
    return setBoolVec( 'ul',    \$ulFwd, \$ulRev, @_ );
}


sub setBoolVec{
    my( $caller, $fwdVec, $revVec, $seqSegs, @coord1_ )  =  @_;

    ref( $seqSegs )  =~ /FastaplSeqsegs/
        or  $seqSegs  =  FastaplSeqsegs->new( $caller, $seqSegs, @coord1_ );

    my $next  =  $seqSegs->rangeItr();
    while(   my( $revP, $d, $e )  =  $next->()   ){
        if( $revP ){   vec( $$revVec, $_, 1 )  =  1   for  $d..$e  }
        else       {   vec( $$fwdVec, $_, 1 )  =  1   for  $d..$e  }
    }

    return $seqSegs;
}#END:  sub setBoolVec.


#  Downcase closed interval seq regions given by arguments.
sub downcase{  changecase( 'downcase', @_ )  }
sub   upcase{  changecase(   'upcase', @_ )  }

#  Change case of sequence segments.
#  Not strand specific.
sub changecase{
    my( $caller, $seqSegs, @coord1_ )  =  @_;

    ref( $seqSegs )  =~ /FastaplSeqsegs/
        or  $seqSegs  =  FastaplSeqsegs->new( $caller, $seqSegs, @coord1_ );

    my $next  =  $seqSegs->rangeItr();
    while(   (undef, $a, $b)  =  $next->()   ){
        if(  $caller =~ /^d/  ){   substr $seq, $a, $b-$a+1,  lc substr $seq, $a, $b-$a+1   }
        else                   {   substr $seq, $a, $b-$a+1,  uc substr $seq, $a, $b-$a+1   }
    }

    return $seqSegs;
}#END: changecase()



my $coAny  =  0;

my @colorXcurNib  =  (0,0);   #  $colorXcurNib[STRAND] Most recently output xterm color nibble on STRAND.

#  For html output.  $colorRcur[STRAND] holds current red color value of STRAND.
#  A value of 256 indicates currently in a non-colored state.
my @colorRcur = my @colorGcur = my @colorBcur =  (256,256);


#  nibble vector for coloring with xterm accordings to @nibbleToColorCode.
my $colorXfwd = '';
my $colorXrev = '';

my $coloredFwd = '';  #Bit vector:  vec( $coloredFwd, $i, 1 )  true if forward strand is colored at position $i.
my $coloredRev = '';  #Bit vector:  vec( $coloredRev, $i, 1 )  true if reverse strand is colored at position $i.
my $colorRfwd = '';   #Byte vector to hold forward strand color vec for red.
my $colorGfwd = '';   #Byte vector to hold forward strand color vec for green.
my $colorBfwd = '';   #Byte vector to hold forward strand color vec for blue.
my $colorRrev = '';   #Byte vector to hold forward strand color vec for red.
my $colorGrev = '';   #Byte vector to hold forward strand color vec for green.
my $colorBrev = '';   #Byte vector to hold forward strand color vec for blue.


sub color{   $coAny = 1;  $htmlP?  &colorHtml : &colorXterm   }


my %colorName_RGB =  qw(  black 000000  red CD0000  green 00CD00  yellow CDCD00  blue 0000CD  magenta CD00CD  cyan 00CDCD  gray E5E5E5  red_bold FF0000  green_bold 00FF00  yellow_bold FFFF00  blue_bold 5C5CFF  magenta_bold FF00FF  cyan_bold 00FFFF  white FFFFFF  );

sub colorHtml{
    my( $color, $seqSegs, @coord1_ )  =  @_;

    ref( $seqSegs )  =~ /FastaplSeqsegs/
        or  $seqSegs  =  FastaplSeqsegs->new( 'color', $seqSegs, @coord1_ );

    $color =~ /^#?[0-9A-F]{6}$/ || $colorName_RGB{$color}   or   die  "Color '$color' not defined";

    my $colorHex  =   (  $color =~ /^#?([0-9A-F]{6})$/?  $1  :  $colorName_RGB{$color}  );

    my( $R,$G,$B )  =   map  {hex}  $colorHex =~ /^#?([0-9A-F][0-9A-F])([0-9A-F][0-9A-F])([0-9A-F][0-9A-F])$/;

    my $next  =  $seqSegs->rangeItr();
    while(   my( $revP, $d, $e )  =  $next->()   ){
        if(  $revP  ){
            for( $d..$e ){
                vec( $coloredRev, $_, 1 )  =  1;
                vec( $colorRrev,  $_, 8 )  =  $R;
                vec( $colorGrev,  $_, 8 )  =  $G;
                vec( $colorBrev,  $_, 8 )  =  $B;
            }
        }else{
            for( $d..$e ){
                vec( $coloredFwd, $_, 1 )  =  1;
                vec( $colorRfwd,  $_, 8 )  =  $R;
                vec( $colorGfwd,  $_, 8 )  =  $G;
                vec( $colorBfwd,  $_, 8 )  =  $B;
            }
        }
    }
    return $seqSegs;
}#END:  sub colorHtml


my @numToHex_  =   map  {sprintf '%02X', $_ }  0..255;
my %colorNameToNibble_  =  (  black=>1, red=>2, green=>3, yellow=>4, blue=>5, magenta=>6, cyan=>7, gray=>8, default=>0  );
my @nibbleToColorCode_ =  ( "\e[39m",   #0 default
                            "\e[30m",   #1 black
                            "\e[31m",   #2 red
                            "\e[32m",   #3 green
                            "\e[33m",   #4 yellow
                            "\e[34m",   #5 blue
                            "\e[35m",   #6 magenta
                            "\e[36m",   #7 cyan
                            "\e[37m",   #8 gray
);

sub colorXterm{
    my( $color, $seqSegs, @coord1_ )  =  (lc $_[0], @_[1..$#_]);

    ref( $seqSegs )  =~ /FastaplSeqsegs/
        or  $seqSegs  =  FastaplSeqsegs->new( 'color', $seqSegs, @coord1_ );

    $color   or   die  'sub color expected color argument';
    exists $colorNameToNibble_{$color}
        or   die  "when using xterm, sub color expected a color in {black,red,green,yellow,blue,magenta,cyan,gray}, but got '$color'";
    my $colorCode  =  $colorNameToNibble_{$color};

    my $next  =  $seqSegs->rangeItr();
    while(   my( $revP, $d, $e )  =  $next->()   ){
        if( $revP ){   vec( $colorXrev, $_, 4 )  =  $colorCode   for  $d..$e   }
        else        {   vec( $colorXfwd, $_, 4 )  =  $colorCode   for  $d..$e   }
    }
    return $seqSegs;
}#END:  sub colorXterm.


sub _colorizeXterm_1Seg{
    my( $seqSegs )  =  @_;

    my $next  =  $seqSegs->rangeItr();
    while(   my( $revP, $d, $e )  =  $next->()   ){
        if(  $revP  ){
            for my $i ($d..$e){
                my $c  =  lc  substr $seq, $i, 1;
                $c =~ tr/ACGTUNSWRYKMBDHVacgtunswrykmbdhv/TGCAANSWYRMKVHDBtgcaanswyrmkvhdb/;
                my $colorSpec  =  $color{$c}   or   die  "color not defined for character seq[$i] = '$c'";
                my ( $color, $makeBold )   =  split '_', $colorSpec;
                vec( $colorXrev, $i, 4 )  =  $colorNameToNibble_{$color};
                vec( $boRev    , $i, 1 )  =  1    if $makeBold;
            }
        }else{
            for my $i ($d..$e){
                my $c  =  lc  substr $seq, $i, 1;
                my $colorSpec  =  $color{$c}   or   die  "color not defined for character seq[$i] = '$c'";
                my ( $color, $makeBold )   =  split '_', $colorSpec;
                vec( $colorXfwd, $i, 4 )  =  $colorNameToNibble_{$color};
                vec( $boFwd    , $i, 1 )  =  1    if $makeBold;
            }
        }
    }
}#END: _colorizeXterm_1Seg

my %nameToR =  (  black=>  0,  red=>205,  green=>  0,  yellow=>205,  blue=>  0,  magenta=>205,
                  cyan=>  0,  gray=>229,  red_bold=>255,  green_bold=>  0,  yellow_bold=>255,
                  blue_bold=> 92,  magenta_bold=>255,  cyan_bold=>  0,  white=>255   );
my %nameToG =  (  black=>  0,  red=>  0,  green=>205,  yellow=>205,  blue=>  0,  magenta=>  0,
                  cyan=>205,  gray=>229,  red_bold=>  0,  green_bold=>255,  yellow_bold=>255,
                  blue_bold=> 92,  magenta_bold=>  0,  cyan_bold=>255,  white=>255    );
my %nameToB =  (  black=>  0,  red=>  0,  green=>  0,  yellow=>  0,  blue=>205,  magenta=>250,
                  cyan=>205,  gray=>229,  red_bold=>  0,  green_bold=>  0,  yellow_bold=>  0,
                  blue_bold=>255,  magenta_bold=>255,  cyan_bold=>255,  white=>255    );

sub _colorizeHtml_1Seg{
    my( $seqSegs )  =  @_;

    my $next  =  $seqSegs->rangeItr();
    while(   my( $revP, $d, $e )  =  $next->()   ){
        if(  $revP  ){
            for my $i ($a..$b){
                my $c   =   lc  substr $seq, $i, 1;
                $c =~ tr/ACGTUNSWRYKMBDHVacgtunswrykmbdhv/TGCAANSWYRMKVHDBtgcaanswyrmkvhdb/;
                my $colorName  =  $color{$c}   or   die  "color not defined for character seq[$i] = '$c'";
                if(  $colorName eq 'default'  ){
                    vec( $coloredRev, $i, 1 )  =  0;
                }
                else{
                    vec( $coloredRev, $i, 1 )  =  1;
                    vec( $colorRrev,  $i, 8 )  =  $nameToR{$colorName};
                    vec( $colorGrev,  $i, 8 )  =  $nameToG{$colorName};
                    vec( $colorBrev,  $i, 8 )  =  $nameToB{$colorName};
                }
            }
        }else{  #Reverse strand
            for my $i ($a..$b){
                my $c   =   lc  substr $seq, $i, 1;
                my $colorName  =  $color{$c}   or   die  "color not defined for character seq[$i] = '$c'";
                if(  $colorName eq 'default'  ){
                    vec( $coloredFwd, $i, 1 )  =  0;
                }
                else{
                    vec( $coloredFwd, $i, 1 )  =  1;
                    vec( $colorRfwd,  $i, 8 )  =  $nameToR{$colorName};
                    vec( $colorGfwd,  $i, 8 )  =  $nameToG{$colorName};
                    vec( $colorBfwd,  $i, 8 )  =  $nameToB{$colorName};
                }
            }
        }
    }
}#END: _colorizeHtml_1Seg{


sub colorize{
    my( $seqSegs, @coord1_ )  =   @_?  @_  :  (0, len);

    ref( $seqSegs )  =~ /FastaplSeqsegs/
        or  $seqSegs  =  FastaplSeqsegs->new( 'color', $seqSegs, @coord1_ );

    $coAny = 1;

    $htmlP?  _colorizeHtml_1Seg $seqSegs  :  _colorizeXterm_1Seg $seqSegs;

    return $seqSegs;
}



my @frame_colorRGB   =  qw(0 CD0000 CDCD00 0000CD);
my @frame_colorCode  =    (39,   31,   33,   34  );


sub ORFtrackAny{
    my $strand  =  shift;  #True for reverse strand, otherwise forward strand.
    # Remaining @_[0,1] should be coordinantes of seqSeg, in either (up,dw) or (dw,up) order.

    my @range  =  seqSeg_to_numRanges @_;
    while(  ($a, $b)  =  splice @range, 0, 2  ){
        if(  $strand  ){    vec( $ORFtrackFrameRev, $_, 4 )  &&  return 1   for $a..$b    }
        else           {    vec( $ORFtrackFrameFwd, $_, 4 )  &&  return 1   for $a..$b    }
    }
    return 0;
}


sub ORFtrack{  shift @_?  &ORFtrackRev : &ORFtrackFwd  }


sub ORFtrackFwd{
    my $i  = $_[0];   $i += len   if $i < 0;

    my $f  =  vec $ORFtrackFrameFwd, $i, 4
        or  return ' ';                            #<---  SUBROUTINE EXIT POINT

    my $codon;
    if(  $i > 0  &&  $i < fin  ){
        $codon   =  substr $seq, $i-1, 3;
    }else{
        $seqIsCirc   or   return ' ';              #<---  SUBROUTINE EXIT POINT
        $codon =  ($i == 0)
            ?  substr( $seq, -1 )  .  substr( $seq, 0, 2 )
            :  substr( $seq, -2 )  .  substr( $seq, 0, 1 );
    }

    my $lcCodon  =  lc $codon;
    my $aa  =  $codonAA{$lcCodon};
    $aa =~ tr/[A-Z]/[a-z]/    if $masking_flag_  &&  $codon ne $lcCodon;

    return $aa    unless $ORF_colorize;            #<---  SUBROUTINE EXIT POINT


    if(  $f > 3  ){  #If starting codon
        $f  -=  4;
        return(  $htmlP
                 ?  qq(<U><FONT color="#$frame_colorRGB[$f]">${aa}</FONT></U>)
                 :  "\e[4m\e[$frame_colorCode[$f]m${aa}\e[39m\e[24m"   );
    }
    else{
        return(  $htmlP
                 ?  qq(<FONT color="#$frame_colorRGB[$f]">${aa}</FONT>)
                 :  "\e[$frame_colorCode[$f]m$aa\e[39m"   );
    }
}#END: ORFtrackFwd


sub ORFtrackRev{
    my $i  = $_[0];   $i += len   if $i < 0;

    my $f   =   vec  $ORFtrackFrameRev, $i, 4
        or  return ' ';                            #<---  SUBROUTINE EXIT POINT

    my $codon;
    if(  $i > 0  &&  $i < fin  ){
        $codon  =  substr $seq, $i-1, 3;
    }else{
        $seqIsCirc   or   return ' ';              #<---  SUBROUTINE EXIT POINT
        $codon =  ($i == 0)
            ?  substr( $seq, -1 )  .  substr( $seq, $i, 2 )
            :  substr( $seq, -2 )  .  substr( $seq, $i, 1 );
    }
    $codon  =  reverse $codon;
    $codon =~ tr/ACGTUNSWRYKMBDHVacgtunswrykmbdhv/TGCAANSWYRMKVHDBtgcaanswyrmkvhdb/;
    my $lcCodon  =  lc $codon;
    my $aa  =  $codonAA{$lcCodon};
    $aa =~ tr/[A-Z]/[a-z]/    if  $masking_flag_  &&  $codon ne $lcCodon;

    return $aa    unless $ORF_colorize;            #<---  SUBROUTINE EXIT POINT


    if(  $f > 3  ){  #If starting codon
        $f  -=  4;
        return(  $htmlP
                 ?   qq(<U><FONT color="#$frame_colorRGB[$f]">${aa}</FONT></U>)
                 :   "\e[4m\e[$frame_colorCode[$f]m${aa}\e[39m\e[24m"   );
    }else{
        return(  $htmlP
                 ?  qq(<FONT color="#$frame_colorRGB[$f]">${aa}</FONT>)
                 :  "\e[$frame_colorCode[$f]m${aa}\e[39m"   );
    }
}#END: ORFtrackRev




#  Internal subroutines
sub dieHtmlOptionMalformed($){
    my $message  =  $_[0];
    say2  "fastapl, error in --html-tags argument '$opt_{htmlTags}' $message";
    exit 64;
}


my $htmlOpenStr  =    '<PRE>';   #String to use when opening html output (e.g. <HTML><BODY><PRE>).
my $htmlClozStr  =   '</PRE>';   #String to use when closing html output (e.g. </PRE></BODY></HTML>).

#  Set $htmlOpenStr, $htmlOpenStr and $termType based on $opt{htmlTags}
sub init_htmlOpenClozStrs{

    $termType_  =  'html';

    my @directive   =   sort  uniq  split ',' => $opt_{htmlTags};

    my $body_handled_once  =  0;
    for(  @directive  ){
        if(  /^body-inverse$/i  ){
            dieHtmlOptionMalformed  'html BODY tag stipulated twice'   if $body_handled_once++;
            $htmlOpenStr  .=    '<BODY bgcolor="#000000" text="#FFFFFF">';
            $htmlClozStr   =   '</BODY>'  .  $htmlClozStr;
        }elsif(  /^body$/i      ){
            dieHtmlOptionMalformed  'html BODY tag stipulated twice'   if $body_handled_once++;
            $htmlOpenStr   =    '<BODY>'  .  $htmlOpenStr;
            $htmlClozStr  .=   '</BODY>';
        }elsif(  /^html$/i      ){
                $htmlOpenStr   =    '<HTML>'  .  $htmlOpenStr;
                $htmlClozStr  .=   '</HTML>';
        }elsif(  /^no-pre$/i    ){
                $htmlOpenStr  =~   s|<PRE>||;
                $htmlClozStr  =~  s|</PRE>||;
        }
        else{
            dieHtmlOptionMalformed  "unknown html directive '$_' ∉ {html,body,body-inverse,no-pre}";
        }
    }
}



#  HTML_do( open|cloz )
#  Output opening or closing tags when appropriate.
sub HTMLdo{
    state $htmlOpened  =  0;

    $htmlOpened++ || print $htmlOpenStr   if  $_[0] eq 'open';
    $htmlOpened   && print $htmlClozStr   if  $_[0] eq 'cloz';
}


sub anyMarkup(){  $blAny || $boAny || $ivAny || $ulAny || $coAny  }


sub clearMarkup(){
    $blAny = $boAny  = $ivAny = $ulAny = $coAny = $ORFtrackAny[0] = $ORFtrackAny[1] = 0;
    $blFwd = $blRev = $boFwd = $boRev = $ivFwd = $ivRev = $ulFwd = $ulRev = $ORFtrackFrameFwd = $ORFtrackFrameRev = '';
    $coloredFwd = $coloredRev = $colorRfwd = $colorGfwd = $colorBfwd = $colorRrev = $colorGrev = $colorBrev = $colorXfwd = $colorXrev = '';
}



#  ──────────────────────────────  Print Seq Functions  ──────────────────────────────
#  Print $seq interval ($up, $dw) with appropriate length lines.


#  ─────  Vars shared by print sequence functions  ─────
my( $ps_up, $ps_dw );  #Current seqSeg to print.

my( $ps_start, $ps_past);  #SUBTLY DIFFERENT THAN $ps_up, $ps_dw.  Defined as:
#  ($ps_start, $ps_past) = ($ps_up,   $ps_dw )    when $ps_up ≦ $ps_dw,
#  ($ps_start, $ps_past) = ($ps_dw-1, $ps_up-1)   when $ps_dw > $ps_dw,
#  so that:   for(  $i = $ps_start;  $i != $ps_past;  $i += $ps_dir  )
#  works for both directions.

my $csq;           #Complemented, but not reversed, copy of $seq.
my $ps_ds;         #Print both strands.
my $ps_dir;        #Direction of printing: +1 if printing along forward strand, otherwise -1
my $ps_len;        #Length of seqSeg to print.
my $ps_lineLen;    #Length of lines to print.
my @ps_range;      #To hold seqSeg_to_directedRanges( $up, $dw );
my $ps_bot;        #Index of bottom sequence: 0 if reverse strand should be printed at top, otherwise 1.
my $ps_top;        #Index of  top   sequence: 1 if reverse strand should be printed at top, otherwise 0.
my $ps_topSeq;     #Reference to sequence to be printed at top.
my $ps_botSeq;     #Reference to sequence to be printed at bottom.
my $ps_1stSegP;    #True iff currently printing first segment.  Used to add blank lines between segs when printing both strands.
my $ps_numResPerLine;  #Number of residues per line.



sub openTagsNotfall{   &openTags     if anyMarkup   }
sub clozTagsNotfall{   &clozTags     if anyMarkup   }
sub clozAllOpenTagsNotfall{   &clozAllOpenTags  if anyMarkup   }



#  print_openTags( s, i )
#  Open tags for closed tags with true value at ith position in strand S.
#  Completely redundant to, but faster than, simply printing the string returned by openTags.
sub print_openTags{
    # $htmlP?  openTagsHtml : openTagsXterm
    my( $s, $i )  =  @_;  #s for strand.
    $i += len   if $i < 0;

    if(  $htmlP  ){
        if(  $s  ){
            #  Reverse Strand.
            do{  $ulOpen[$s] = 1; print '<U>'    }    if  !$ulOpen[$s]  &&  vec( $ulRev, $i, 1 );
            do{  $boOpen[$s] = 1; print '<B>'    }    if  !$boOpen[$s]  &&  vec( $boRev, $i, 1 );
            do{  $blOpen[$s] = 1; print '<BLINK>'}    if  !$blOpen[$s]  &&  vec( $blRev, $i, 1 );
            my( $R, $G, $B )  =  (  vec( $colorRrev, $i, 8 ),
                                    vec( $colorGrev, $i, 8 ),
                                    vec( $colorBrev, $i, 8 )  );
            if(  $R != $colorRcur[$s]  or  $G != $colorGcur[$s]  or  $B != $colorBcur[$s]  ){
                if(  vec( $coloredRev, $i, 1 )   ){
                    print '<FONT color="', $numToHex_[ $R ], $numToHex_[ $G ], $numToHex_[ $B ], '">';
                    ($colorRcur[$s], $colorGcur[$s], $colorBcur[$s])  =  ($R, $G, $B);
                }
            }
        }else{  #Forward Strand.
            do{  $ulOpen[$s] = 1; print '<U>'    }    if  !$ulOpen[$s]  &&  vec( $ulFwd, $i, 1 );
            do{  $boOpen[$s] = 1; print '<B>'    }    if  !$boOpen[$s]  &&  vec( $boFwd, $i, 1 );
            do{  $blOpen[$s] = 1; print '<BLINK>'}    if  !$blOpen[$s]  &&  vec( $blFwd, $i, 1 );
            my( $R, $G, $B )  =  (  vec( $colorRfwd, $i, 8 ),
                                    vec( $colorGfwd, $i, 8 ),
                                    vec( $colorBfwd, $i, 8 )  );
            if(  $R != $colorRcur[$s]  or  $G != $colorGcur[$s]  or  $B != $colorBcur[$s]  ){
                if(  vec( $coloredFwd, $i, 1 )   ){
                    print '<FONT color="', $numToHex_[ $R ], $numToHex_[ $G ], $numToHex_[ $B ], '">';
                    ($colorRcur[$s], $colorGcur[$s], $colorBcur[$s])  =  ($R, $G, $B);
                }
            }
        }
    }else{   #  xterm
        my $colorNib;
        if(  $s  ){
            #  Reverse Strand.
            do{  $ulOpen[$s] = 1; print "\e[4m"  }    if  !$ulOpen[$s]  &&  vec( $ulRev, $i, 1 );
            do{  $boOpen[$s] = 1; print "\e[1m"  }    if  !$boOpen[$s]  &&  vec( $boRev, $i, 1 );
            do{  $blOpen[$s] = 1; print "\e[5m"  }    if  !$blOpen[$s]  &&  vec( $blRev, $i, 1 );
            do{  $ivOpen[$s] = 1; print "\e[7m"  }    if  !$ivOpen[$s]  &&  vec( $ivRev, $i, 1 );
            $colorNib  =  vec( $colorXrev, $i, 4 );
        }else{
            #  Forward Strand.
            do{  $ulOpen[$s] = 1; print "\e[4m"  }    if  !$ulOpen[$s]  &&  vec( $ulFwd, $i, 1 );
            do{  $boOpen[$s] = 1; print "\e[1m"  }    if  !$boOpen[$s]  &&  vec( $boFwd, $i, 1 );
            do{  $blOpen[$s] = 1; print "\e[5m"  }    if  !$blOpen[$s]  &&  vec( $blFwd, $i, 1 );
            do{  $ivOpen[$s] = 1; print "\e[7m"  }    if  !$ivOpen[$s]  &&  vec( $ivFwd, $i, 1 );
            $colorNib  =  vec( $colorXfwd, $i, 4 );
        }
        if(   $colorNib != $colorXcurNib[$s]   ){
            print  $nibbleToColorCode_[$colorNib];
            $colorXcurNib[$s]  =  $colorNib;
        }
    }
}#END: print_openTags.


#  openTags( s, i )
#  Return open tags string for closed tags with true value at ith position in strand S.
sub openTags{
    # $htmlP?  openTagsHtml : openTagsXterm
    my( $s, $i )  =  @_;  #s for strand.
    $i += len   if $i < 0;

    my $retVal  =  '';
    if(  $htmlP  ){
        if(  $s  ){
            #  Reverse Strand.
            do{  $ulOpen[1] = 1; $retVal .= '<U>'    }    if  !$ulOpen[1]  &&  vec( $ulRev, $i, 1 );
            do{  $boOpen[1] = 1; $retVal .= '<B>'    }    if  !$boOpen[1]  &&  vec( $boRev, $i, 1 );
            do{  $blOpen[1] = 1; $retVal .= '<BLINK>'}    if  !$blOpen[1]  &&  vec( $blRev, $i, 1 );
            my( $R, $G, $B )  =  (  vec( $colorRrev, $i, 8 ),
                                    vec( $colorGrev, $i, 8 ),
                                    vec( $colorBrev, $i, 8 )  );
            if(  $R != $colorRcur[1]  or  $G != $colorGcur[1]  or  $B != $colorBcur[1]  ){
                if(  vec( $coloredRev, $i, 1 )   ){
                    $retVal  .=  '<FONT color="' . $numToHex_[$R] . $numToHex_[$G] . $numToHex_[$B] . '">';
                    ($colorRcur[1], $colorGcur[1], $colorBcur[1])  =  ($R, $G, $B);
                }
            }
        }
        else{
            #  Forward Strand.
            do{  $ulOpen[0] = 1; $retVal .= '<U>'    }    if  !$ulOpen[0]  &&  vec( $ulFwd, $i, 1 );
            do{  $boOpen[0] = 1; $retVal .= '<B>'    }    if  !$boOpen[0]  &&  vec( $boFwd, $i, 1 );
            do{  $blOpen[0] = 1; $retVal .= '<BLINK>'}    if  !$blOpen[0]  &&  vec( $blFwd, $i, 1 );
            my( $R, $G, $B )  =  (  vec( $colorRfwd, $i, 8 ),
                                    vec( $colorGfwd, $i, 8 ),
                                    vec( $colorBfwd, $i, 8 )  );
            if(  $R != $colorRcur[0]  or  $G != $colorGcur[0]  or  $B != $colorBcur[0]  ){
                if(  vec( $coloredFwd, $i, 1 )   ){
                    $retVal  .=  '<FONT color="' . $numToHex_[$R] . $numToHex_[$G] . $numToHex_[$B] . '">';
                    ($colorRcur[0], $colorGcur[0], $colorBcur[0])  =  ($R, $G, $B);
                }
            }
        }
    }else{   #  xterm
        my $colorNib;
        if(  $s  ){
            #  Reverse strand.
            do{  $ulOpen[1] = 1; $retVal .= "\e[4m"  }    if  !$ulOpen[1]  &&  vec( $ulRev, $i, 1 );
            do{  $boOpen[1] = 1; $retVal .= "\e[1m"  }    if  !$boOpen[1]  &&  vec( $boRev, $i, 1 );
            do{  $blOpen[1] = 1; $retVal .= "\e[5m"  }    if  !$blOpen[1]  &&  vec( $blRev, $i, 1 );
            do{  $ivOpen[1] = 1; $retVal .= "\e[7m"  }    if  !$ivOpen[1]  &&  vec( $ivRev, $i, 1 );
            $colorNib  =  vec( $colorXrev, $i, 4 );
        }else{
            #  Forward strand.
            do{  $ulOpen[0] = 1; $retVal .= "\e[4m"  }    if  !$ulOpen[0]  &&  vec( $ulFwd, $i, 1 );
            do{  $boOpen[0] = 1; $retVal .= "\e[1m"  }    if  !$boOpen[0]  &&  vec( $boFwd, $i, 1 );
            do{  $blOpen[0] = 1; $retVal .= "\e[5m"  }    if  !$blOpen[0]  &&  vec( $blFwd, $i, 1 );
            do{  $ivOpen[0] = 1; $retVal .= "\e[7m"  }    if  !$ivOpen[0]  &&  vec( $ivFwd, $i, 1 );
            $colorNib  =  vec( $colorXfwd, $i, 4 );
        }
        if(   $colorNib != $colorXcurNib[$s]   ){
            $retVal  .=  $nibbleToColorCode_[$colorNib];
            $colorXcurNib[$s]  =  $colorNib;
        }
    }
    return $retVal;
}#END:  openTags.


#  print_clozTags( s, i )
#  Print close tags for open tags with false value at ith position in strand S.
#  Completely redundant to, but significantly faster than, simply printing the string returned by clozTags.
sub print_clozTags{
    my( $s, $i )  =  @_;  #s for strand.
    $i += len   if $i < 0;

    my $n  =  $i+$ps_dir;  #n for next.
    if(  $seqIsCirc  ){
        $n = 0      if  $n == len;
        $n = len-1  if  $n < 0;
    }
    #Edge case of $n = len or $n = -1 for linear seq is not a problem.  Although out of range, vec silently returns 0 in that case.

    # $htmlP?  clozTagsHtml : clozTagsXterm;
    if(  $htmlP  ){
        if(  $s  ){#  Reverse strand.
            if(   $colorRcur[1] < 256   and
                  (!vec $coloredRev, $n, 1   or
                   vec $colorRrev, $i, 8  !=  vec $colorRrev, $n, 8   or
                   vec $colorGrev, $i, 8  !=  vec $colorGrev, $n, 8   or
                   vec $colorBrev, $i, 8  !=  vec $colorBrev, $n, 8)  ){
                print '</FONT>';   $colorRcur[1] = 256;
            }
            do{  $blOpen[1] = 0; print '</BLINK>'}    if  $blOpen[1] && !vec( $blRev, $n, 1 );
            do{  $boOpen[1] = 0; print '</B>'    }    if  $boOpen[1] && !vec( $boRev, $n, 1 );
            do{  $ulOpen[1] = 0; print '</U>'    }    if  $ulOpen[1] && !vec( $ulRev, $n, 1 );
        }
        else{#  Forward strand.
            if(   $colorRcur[0] < 256   and
                  (!vec $coloredFwd, $n, 1   or
                   vec $colorRfwd, $i, 8  !=  vec $colorRfwd, $n, 8   or
                   vec $colorGfwd, $i, 8  !=  vec $colorGfwd, $n, 8   or
                   vec $colorBfwd, $i, 8  !=  vec $colorBfwd, $n, 8)  ){
                print '</FONT>';   $colorRcur[0] = 256;
            }
            do{  $blOpen[0] = 0; print '</BLINK>'}    if  $blOpen[0] && !vec( $blFwd, $n, 1 );
            do{  $boOpen[0] = 0; print '</B>'    }    if  $boOpen[0] && !vec( $boFwd, $n, 1 );
            do{  $ulOpen[0] = 0; print '</U>'    }    if  $ulOpen[0] && !vec( $ulFwd, $n, 1 );
        }
    }
    else{   #  xterm.
        #  Implementation note: This part would be simpler if the turn off bold code (21) was better supported.
        #  Because it is poorly supported I use the reset code (0) instead, which also wipes out the underline state.
        #  Color tags do not need to be closed here because openTags opens the default color when appropriate.
        if(  $s  ){#  Reverse strand.
            if(   $boOpen[1]  &&  !vec( $boRev, $n, 1 )   ){
                print "\e[0m";
                $boOpen[1] = $ulOpen[1] = $colorXcurNib[1] = 0;
            }
            do{  $ivOpen[1] = 0;  print "\e[27m"  }    if  $ivOpen[1] && !vec( $ivRev, $n, 1 );
            do{  $blOpen[1] = 0;  print "\e[25m"  }    if  $blOpen[1] && !vec( $blRev, $n, 1 );
            do{  $ulOpen[1] = 0;  print "\e[24m"  }    if  $ulOpen[1] && !vec( $ulRev, $n, 1 );
        }else{#  Forward strand.
            if(   $boOpen[0]  &&  !vec( $boFwd, $n, 1 )   ){
                print "\e[0m";
                $boOpen[0] = $ulOpen[0] = $colorXcurNib[0] = 0;
            }
            do{  $ivOpen[0] = 0;  print "\e[27m"  }    if  $ivOpen[0] && !vec( $ivFwd, $n, 1 );
            do{  $blOpen[0] = 0;  print "\e[25m"  }    if  $blOpen[0] && !vec( $blFwd, $n, 1 );
            do{  $ulOpen[0] = 0;  print "\e[24m"  }    if  $ulOpen[0] && !vec( $ulFwd, $n, 1 );
        }
    }
}#END: print_clozTags()


#  clozTags( s, i )
#  Return close tags string for open tags with false value at ith position in strand S.
sub clozTags{
    my( $s, $i )  =  @_;  #s for strand.
    $i += len   if $i < 0;

    my $n  =  $i+$ps_dir;  #n for next.

    if(  $seqIsCirc  ){
        $n = 0      if  $n == len;
        $n = len-1  if  $n < 0;
    }
    #Edge case of $n = len or $n = -1 for linear seq is not a problem.  Although out of range, vec silently returns 0 in that case.

    my $retVal  =  '';
    # $htmlP?  clozTagsHtml : clozTagsXterm;
    if(  $htmlP  ){
        if(  $s  ){#  Reverse strand.
            if(   $colorRcur[1] < 256   and
                  (!vec $coloredRev, $n, 1   or
                   vec $colorRrev, $i, 8  !=  vec $colorRrev, $n, 8   or
                   vec $colorGrev, $i, 8  !=  vec $colorGrev, $n, 8   or
                   vec $colorBrev, $i, 8  !=  vec $colorBrev, $n, 8)  ){
                $retVal .= '</FONT>';   $colorRcur[1] = 256;
            }
            do{  $blOpen[1] = 0; $retVal .= '</BLINK>'}    if  $blOpen[1] && !vec( $blRev, $n, 1 );
            do{  $boOpen[1] = 0; $retVal .= '</B>'    }    if  $boOpen[1] && !vec( $boRev, $n, 1 );
            do{  $ulOpen[1] = 0; $retVal .= '</U>'    }    if  $ulOpen[1] && !vec( $ulRev, $n, 1 );
        }
        else{#  Forward strand.
            if(   $colorRcur[0] < 256   and
                  (!vec $coloredFwd, $n, 1   or
                   vec $colorRfwd, $i, 8  !=  vec $colorRfwd, $n, 8   or
                   vec $colorGfwd, $i, 8  !=  vec $colorGfwd, $n, 8   or
                   vec $colorBfwd, $i, 8  !=  vec $colorBfwd, $n, 8)  ){
                $retVal .= '</FONT>';   $colorRcur[0] = 256;
            }
            do{  $blOpen[0] = 0; $retVal .= '</BLINK>'}    if  $blOpen[0] && !vec( $blFwd, $n, 1 );
            do{  $boOpen[0] = 0; $retVal .= '</B>'    }    if  $boOpen[0] && !vec( $boFwd, $n, 1 );
            do{  $ulOpen[0] = 0; $retVal .= '</U>'    }    if  $ulOpen[0] && !vec( $ulFwd, $n, 1 );
        }
    }
    else{   #xterm.
        #  Implementation note: This part would be simpler if the turn off bold code (21) was better supported.
        #  Because it is poorly supported I use the reset code (0) instead, which also wipes out the underline state.
        #  Color tags do not need to be closed here because openTags opens the default color when appropriate.
        if(  $s  ){#  Reverse strand.
            if(   $boOpen[1]  &&  !vec( $boRev, $n, 1 )   ){
                $retVal .= "\e[0m";
                $boOpen[1] = $ulOpen[1] = $colorXcurNib[1] = 0;
            }
            do{  $ivOpen[1] = 0;  $retVal .= "\e[27m" }    if  $ivOpen[1] && !vec( $ivRev, $n, 1 );
            do{  $blOpen[1] = 0;  $retVal .= "\e[25m" }    if  $blOpen[1] && !vec( $blRev, $n, 1 );
            do{  $ulOpen[1] = 0;  $retVal .= "\e[24m" }    if  $ulOpen[1] && !vec( $ulRev, $n, 1 );
        }else{#  Forward strand.
            if(   $boOpen[0]  &&  !vec( $boFwd, $n, 1 )   ){
                $retVal .= "\e[0m";
                $boOpen[0] = $ulOpen[0] = $colorXcurNib[0] = 0;
            }
            do{  $ivOpen[0] = 0;  $retVal .= "\e[27m" }    if  $ivOpen[0] && !vec( $ivFwd, $n, 1 );
            do{  $blOpen[0] = 0;  $retVal .= "\e[25m" }    if  $blOpen[0] && !vec( $blFwd, $n, 1 );
            do{  $ulOpen[0] = 0;  $retVal .= "\e[24m" }    if  $ulOpen[0] && !vec( $ulFwd, $n, 1 );
        }
    }
    return $retVal;
}#END: clozTags()



#  Unconditionally close all tags on STRAND.
#  Completely redundant to, but significantly faster than, simply printing the string returned by clozAllOpenTags.
sub print_clozAllOpenTags{
    my( $s )  =  @_;  #strand
    if(  $htmlP  ){
        do{  $blOpen[$s] = 0; print '</BLINK>'}   if $blOpen[$s];
        do{  $boOpen[$s] = 0; print '</B>'    }   if $boOpen[$s];
        do{  $ulOpen[$s] = 0; print '</U>'    }   if $ulOpen[$s];
        if(  $colorRcur[$s] < 256  ){
            print '</FONT>';   $colorRcur[$s]  =  256;
        }
    }else{
        do{  $ivOpen[$s] = 0; print "\e[27m"  }   if $ivOpen[$s];
        do{  $blOpen[$s] = 0; print "\e[25m"  }   if $blOpen[$s];
        if(  $boOpen[$s] || $ulOpen[$s]  ){
            print "\e[0m";
            $ulOpen[$s] = $boOpen[$s] = 0;
        }
        if(  $colorXcurNib[$s]  ){
            print "\e[39m";
            $colorXcurNib[$s] = 0;
        }
    }
}#END: print_clozAllOpenTags



#  Return string to close all open tags on STRAND.
sub clozAllOpenTags($){
    my( $s )  =  @_;  #strand
    my $retVal  =  '';
    if(  $htmlP  ){
        do{  $blOpen[$s] = 0; $retVal .= '</BLINK>'}   if $blOpen[$s];
        do{  $boOpen[$s] = 0; $retVal .= '</B>'    }   if $boOpen[$s];
        do{  $ulOpen[$s] = 0; $retVal .= '</U>'    }   if $ulOpen[$s];
        if(  $colorRcur[$s] < 256  ){
            $retVal .= '</FONT>';   $colorRcur[$s]  =  256;
        }
    }else{
        do{  $ivOpen[$s] = 0; $retVal .= "\e[27m"  }   if $ivOpen[$s];
        do{  $blOpen[$s] = 0; $retVal .= "\e[25m"  }   if $blOpen[$s];
        if(  $boOpen[$s] || $ulOpen[$s]  ){
            $retVal .= "\e[0m";
            $ulOpen[$s] = $boOpen[$s] = 0;
        }
        if(  $colorXcurNib[$s]  ){
            $retVal .= "\e[39m";
            $colorXcurNib[$s] = 0;
        }
    }
    return $retVal;
}#END: clozAllOpenTags



#  Print current fasta record to STDOUT.
sub pr(){
    state $initialized  =  0;

    die  'pr() called in non-void context'   if defined wantarray;

    if(  !$initialized++  ){
        HTMLdo 'open'   if $htmlP;
        print  $commentAtStartOfStream;
    }

    if(  $id  &&  $id ne $idCopy  ){
        $head  =~   s/^(\s*) \Q$idCopy\E (?!\S)  /$1$id/x;
        $idCopy  =  $id;
    }

    print  ">$head\n$comment";

    $opt_{printEachRecDS_P}?  Ps()  :  ps();
}#END: pr()



#  ─────  Declare auxiliary functions to subs ps,Ps  ─────
sub ps3();


#  Print sequence length string looking like:     .    1    .    2…
#  Based on globals $ps_start, $ps_past, $ps_dir and $ps_dw
#  Or just the starting position if the line is too short.
sub sayRulerTrack{
    use integer;   #Use integer division in this function.

    my ($start, $past)  =  (@_, $ps_start, $ps_past);

    my $lineLen  =  abs( $start - $past );
    my $rulerTrack   =   ' '  x  $lineLen;

    my $i  =  $start;
    if(  $ps_dir == 1  ){
        ++$i  while(  $i % 5  );  #Proceed to multiple of 5.
        my $notDivBy10  =  $i % 10;
        while(  $i < $past  ){
            if( $notDivBy10 ){  substr $rulerTrack, $i-$start, 1, '.'  }
            else{
                my $rem00  =  $i % 100;
                if(  $rem00 ){
                    substr $rulerTrack, $i-$start, 1, $rem00/10;
                }else{
                    my $div00  =  $i / 100;
                    substr $rulerTrack, $i-$start, length $div00, $div00;
                }
            }
            $notDivBy10  =  !$notDivBy10;
            $i += 5;
        }
        substr $rulerTrack, 0, length($start)+1, $start.' ';
    }
    else{
        --$i  while(  $i % 5  );  #Proceed to multiple of 5.
        my $notDivBy10  =  $i % 10;
        while(  $i > $past  ){
            if( $notDivBy10 ){  substr $rulerTrack, $start-$i, 1, '.'  }
            else{
                my $rem00  =  $i % 100;
                if(  $rem00 ){
                    substr $rulerTrack, $start-$i, 1, $rem00/10;
                }else{
                    my $div00  =  $i / 100;
                    substr $rulerTrack, $start-$i, length $div00, $div00;
                }
            }
            $notDivBy10  =  !$notDivBy10;
            $i -= 5;
        }
        substr $rulerTrack, 0, length($start)+1, $start.'-';
    }

    say  $rulerTrack;
}#END:  sayRulerTrack


#  ps( SEQSEGS )  Prints SEQSEGS.
#  SEQSEGS can take the form of a list of paired positions or a reference to FastaplSeqsegs object
#  Defaults:  ps(    ) -->  ps( 0,  len )
#             ps(  3 ) -->  ps( 3,  len )
#             ps( -3 ) -->  ps( 3, -len )
sub ps{
    $ps_ds  =  0;
    $csq  =  '';
    &ps2;
}


#  Like ps, but prints sequence as double stranded DNA.
sub Ps{
    HTMLdo 'open'   if $htmlP;
    $ps_ds  =  1;
    ($csq = $seq)  =~  tr/ACGTUNSWRYKMBDHVacgtunswrykmbdhv/TGCAANSWYRMKVHDBtgcaanswyrmkvhdb/;
    &ps2;
}



#  Pass $seq intervals in @seqSegs onto ps3
sub ps2{
    @_   or   @_ = (0,len);

    my $ref  =  ref $_[0];
    if(  !$ref  ){
        my @seg  =  @_;

        #  Fill in default value of $dw as needed.
        push @seg,  $seg[0] < 0? -&len : len    if @seg == 1;

        $ps_1stSegP = 1;   ps3  if    ($ps_up, $ps_dw)  =  splice @seg, 0, 2;
        $ps_1stSegP = 0;   ps3  while ($ps_up, $ps_dw)  =  splice @seg, 0, 2;
    }
    else{
        $ref eq 'FastaplSeqsegs'   or   die  $ps_ds? 'Ps:' : 'ps:', "passed unknown ref type '$ref'";
        my $nextPair  =  $_[0]->itr();
        $ps_1stSegP = 1;   ps3  if    ($ps_up, $ps_dw)  =  $nextPair->();
        $ps_1stSegP = 0;   ps3  while ($ps_up, $ps_dw)  =  $nextPair->();
    }
}


#  Print $seq interval [UP, DW) in appropriate format.
#  Does range checking.
sub ps3(){

    return  if $ps_up == $ps_dw;   #Printing empty seg is noop.
    @ps_range  =   seqSeg_to_directedRanges  $ps_up,  $ps_dw,  qw(ps Ps)[$ps_ds];


    #  from $ps_up and  $ps_dw
    #  set $ps_start, $ps_past, $ps_top, $ps_bot, $ps_dir, $ps_len, $ps_topSeq, $ps_botSeq and $csq;
    if(  $ps_dw < $ps_up  ){
        $ps_start  =  $ps_up-1;   $ps_past  =  $ps_dw-1;
        $ps_top = 1;  $ps_bot = 0;  $ps_dir = -1;
        $ps_len  =  $ps_up - $ps_dw;
        $ps_topSeq = \$csq;
        $ps_botSeq = \$seq    if $ps_ds;
        $csq   or   ($csq = $seq)  =~  tr/ACGTUNSWRYKMBDHVacgtunswrykmbdhv/TGCAANSWYRMKVHDBtgcaanswyrmkvhdb/;
    }else{
        $ps_start  =  $ps_up;   $ps_past  =  $ps_dw;
        $ps_top = 0;  $ps_bot = 1;  $ps_dir = +1;
        $ps_len  =  $ps_dw - $ps_up;
        $ps_topSeq = \$seq;
        $ps_botSeq = \$csq    if $ps_ds;
    }

    die  "...Omitting long ($_ chars) sequence.\n"
      if -t STDOUT  and   100_000_000  <  ($_ = abs($ps_dw - $ps_up));

    colorize  $ps_up, $ps_dw,  $ps_ds? ($ps_dw, $ps_up) : ()    if $opt_{colorizeP};

    if(  $opt_{printSeqOn1LineP}  ||  ($numSeqLines_ == 1) && !$opt_{seqLineLen}  ){
        if(  !$ps_1stSegP  ){    say ''  if $ps_ds   }#  Blank line between segs.
        goto &printSeq_1line_inlineRuler   if $opt_{rulerInlineP};
        goto &printSeq_1line;
    }

    if(  !$ps_1stSegP  ){    say '';  say ''  if $ps_ds   }#  Blank line(s) between segs.
    $ps_numResPerLine  =   $opt_{seqLineLen}  ||  $maxSeqLineLenSeen_;
    goto &printSeq_multiline_inlineRuler   if $opt_{rulerInlineP};
    goto &printSeq_multiline;

}#END: _ps()



#  Print $seq on one line.
sub printSeq_1line(){

    sayRulerTrack   if $opt_{rulerTrackP};

    if(  $ORFtrackAny[$ps_top]  ){
        for(  my $i = $ps_start;  $i != $ps_past;  $i += $ps_dir  ){
            print ORFtrack $ps_top, $i;
        }
        say '';
    }

    if(  anyMarkup  ){
        {
            for(  my $i = $ps_start;  $i != $ps_past;  $i += $ps_dir  ){
                print_openTags $ps_top, $i;
                print  substr $$ps_topSeq, $i, 1;
                print_clozTags $ps_top, $i;
            }
            print_clozAllOpenTags $ps_top;
            say  '';
        }
        if(  $ps_ds  ){
            for(  my $i = $ps_start;  $i != $ps_past;  $i += $ps_dir  ){
                print_openTags $ps_bot, $i;
                print  substr $$ps_botSeq, $i, 1;
                print_clozTags $ps_bot, $i;
            }
            print_clozAllOpenTags $ps_bot;
            say  '';
        }
    }else{  #No markup
        my $s  =   substr  $$ps_topSeq,  min($ps_up,$ps_dw),  $ps_len;
        say  $ps_top?  scalar reverse $s  :  $s;
        if(  $ps_ds  ){
            $s  =   substr  $$ps_botSeq,  min($ps_up,$ps_dw),  $ps_len;
            say  $ps_top?  scalar reverse $s  :  $s;
        }
    }

    if(  $ps_ds  &&  $ORFtrackAny[$ps_bot]  ){
        for(  my $i = $ps_start;  $i != $ps_past;  $i += $ps_dir  ){
            print ORFtrack $ps_bot, $i;
        }
        say '';
    }


}#END:  printSeq_1line



#  ─────────────────────────  Ruler Functions  ─────────────────────────
my $ps_numWidth;  #Printing width of inline ruler position.

sub set_ps_numWidth{
    my $absMax  =  max( abs($ps_up), abs($ps_dw) ) / 100.0;
    for(  $ps_numWidth = 1;  $absMax >= 10;  $absMax /= 10  ){ ++$ps_numWidth }
}

#  Return inline rule string appropriate for just before from-0 position $i
sub ruleBefore($){
    return  ''                                             if  $_[0] % 10   !=  0;
    return  sprintf " %${ps_numWidth}d ", int $_[0]/100    if  $_[0] % 100  ==  0;
    return  '  '                                           if  $_[0] %  50  ==  0;
    return  ' ';
}

#  Return inline rule string appropriate for just before from-0 position $i
sub ruleBeforeSpacer($){
    return  ''                                             if  $_[0] % 10   !=  0;
    return  ' ' x ($ps_numWidth+2)                         if  $_[0] % 100  ==  0;
    return  '  '                                           if  $_[0] %  50  ==  0;
    return  ' ';
}


#  Return 1 iff fewer than two numbers would be printed inline in the seqseg ps_beg, ps_end.
sub ps_inlineNumsOneOrZero(){
    my $count = 0;
    for(  my $i = $ps_start;  $i != $ps_past;  $i += $ps_dir  ){
        ++$count    if  $i % 100 == 0;
        return 0    if  $count > 1;
    }
    return 1;
}


#  Print $seq on one line with inline ruler.
#  $ds ∈ {0,1}; 1 means print as double stranded DNA.
sub printSeq_1line_inlineRuler(){
    my( $minPos, $maxPos )  =  minmax $ps_up, $ps_dw;

    set_ps_numWidth;
    my( $i, $printedNumCount );

    my $shouldPrintNumAtEnd  =  ps_inlineNumsOneOrZero;

    my $r;

    if(  anyMarkup  ){
        {
            if(   $ORFtrackAny[$ps_top]  &&  ORFtrackAny( $ps_top, $ps_up, $ps_dw )  ){
                $i = $ps_start;
                print ORFtrack $ps_top, $i;

                for(  $i += $ps_dir;  $i != $ps_past;  $i += $ps_dir  ){
                    print  ruleBeforeSpacer $i,  ORFtrack $ps_top, $i;
                }
                say  '';
            }
            {#  Top sequence track.
                $i = $ps_start;
                print_openTags $ps_top, $i;
                print  substr $$ps_topSeq, $i, 1;
                print_clozTags $ps_top, $i;

                for(  $i += $ps_dir;  $i != $ps_past;  $i += $ps_dir  ){
                    $r  =  ruleBefore $i;
                    print clozAllOpenTags $ps_top,  $r   if $r;
                    print_openTags $ps_top, $i;
                    print  substr $$ps_topSeq, $i, 1;
                    print_clozTags $ps_top, $i;
                }

                print_clozAllOpenTags $ps_top;

                print  $ps_top? ' >' : '  ',  $ps_dw     if   $shouldPrintNumAtEnd;
                say  '';
            }
        }
        if(  $ps_ds  ){
            {#  Bottom sequence track.
                $i = $ps_start;
                print_openTags $ps_bot, $i;
                print  substr $$ps_botSeq, $i, 1;
                print_clozTags $ps_bot, $i;

                for(  $i += $ps_dir;  $i != $ps_past;  $i += $ps_dir  ){
                    $r  =  ruleBefore $i;
                    print clozAllOpenTags $ps_bot,  $r   if $r;
                    print_openTags $ps_bot, $i;
                    print  substr $$ps_botSeq, $i, 1;
                    print_clozTags $ps_bot, $i;
                }

                print_clozAllOpenTags $ps_bot;

                print  ' ' x (1 + length $ps_dw)    if   $shouldPrintNumAtEnd;
                say  '';
            }
            if(   $ORFtrackAny[$ps_bot]  &&  ORFtrackAny( $ps_bot, $ps_up, $ps_dw )  ){
                $i = $ps_start;
                print ORFtrack $ps_bot, $i;

                for(  $i += $ps_dir;  $i != $ps_past;  $i += $ps_dir  ){
                    print  ruleBeforeSpacer $i,  ORFtrack $ps_bot, $i;
                }
                say  '';
            }
        }
    }
    else{#  No Markup.
        {
            if(   $ORFtrackAny[$ps_top]  &&  ORFtrackAny( $ps_top, $ps_up, $ps_dw )  ){
                $i = $ps_start;
                print ORFtrack $ps_top, $i;

                for(  $i += $ps_dir;  $i != $ps_past;  $i += $ps_dir  ){
                    print  ruleBeforeSpacer $i,  ORFtrack $ps_top, $i;
                }
                say  '';
            }
            {#  Top sequence track.
                $i = $ps_start;
                print  substr $$ps_topSeq, $i, 1;

                for(  $i += $ps_dir;  $i != $ps_past;  $i += $ps_dir  ){
                    print  ruleBefore $i;
                    print  substr $$ps_topSeq, $i, 1;
                }

                print  $ps_top? ' >' : '  ',  $ps_dw     if   $shouldPrintNumAtEnd;
                say  '';
            }
        }
        if(  $ps_ds  ){
            {#  Bottom sequence track.
                $i = $ps_start;
                print  substr $$ps_botSeq, $i, 1;

                for(  $i += $ps_dir;  $i != $ps_past;  $i += $ps_dir  ){
                    print  ruleBefore $i;
                    print  substr $$ps_botSeq, $i, 1;
                }

                print  ' ' x (1 + length $ps_dw)    if   $shouldPrintNumAtEnd;
                say  '';
            }

            if(   $ORFtrackAny[$ps_bot]  &&  ORFtrackAny( $ps_bot, $ps_up, $ps_dw )  ){
                $i = $ps_start;
                print ORFtrack $ps_bot, $i;

                for(  $i += $ps_start;  $i != $ps_past;  $i += $ps_dir  ){
                    print  ruleBeforeSpacer $i,  ORFtrack $ps_bot, $i;
                }
                say  '';
            }

        }
    }
}#END: printSeq_1line_inlineRuler()


#  Print $seq, wrapping lines at $lineLen characters.
#  $ds={0,1}; 1 means print as double stranded DNA.
#  See printSeq_multiline_inlineRuler for similar function but with inline ruler.
sub printSeq_multiline{
    my $numLineBlocks  =  ceil $ps_len/$ps_numResPerLine;

    my $start = $ps_start;  #To hold start of current line.

    if(  $ps_dir == 1  ){
        for my $lineBlockNum(  1 .. $numLineBlocks  ){
            my $past = min $start+$ps_numResPerLine, $ps_dw;

            sayRulerTrack( $start, $past)  if $opt_{rulerTrackP};

            if(  ORFtrackAny $ps_top, $start, $past-1  ){
                say  ''    if  !$ps_ds;
                print ORFtrack $ps_top, $_    for  $start .. $past-1;
                say  '';
            }

            if(  anyMarkup  ){
                {
                    for ($start..$past-1){
                        my $i = $_;   $i += len  if $i < 0;
                        print   openTags( 0, $i ),  substr( $seq, $i, 1 ),  clozTags( 0, $i );
                    }
                    say  clozAllOpenTags 0;
                }
                if(  $ps_ds  ){
                    for ($start..$past-1){
                        my $i = $_;   $i += len  if $i < 0;
                        print   openTags( 1, $i ),  substr( $csq, $i, 1 ),  clozTags( 1, $i );
                    }
                    say  clozAllOpenTags 1;
                }
            }
            else{  #No markup.
                if(  $start < 0  ){
                    my $beg =  len + $start;
                    if(  $past > 0  ){
                        say  substr( $seq, $beg ), substr $seq, 0, $past;
                        say  substr( $csq, $beg ), substr $csq, 0, $past   if $ps_ds;
                    }else{
                        say  substr $seq, $beg, $past-$start;
                        say  substr $csq, $beg, $past-$start    if $ps_ds;
                    }
                }else{
                    say  substr $seq, $start, $past-$start;
                    say  substr $csq, $start, $past-$start    if $ps_ds;
                }
            }

            if(  $ps_ds  ){
                if(  ORFtrackAny $ps_bot, $start, $past-1  ){
                    print ORFtrack $ps_bot, $_    for  $start .. $past-1;
                    say  '';
                }
                say  ''    if $lineBlockNum < $numLineBlocks;
            }

            $start += $ps_numResPerLine;
        }#END:  for each line block.
    }
    else{#  $ps_dir == -1
        for my $lineBlockNum(  1 .. $numLineBlocks  ){
            my $past  =  max $start-$ps_numResPerLine, $ps_dw-1;

            sayRulerTrack( $start, $past)  if $opt_{rulerTrackP};

            if(  ORFtrackAny $ps_top, $past+1, $start  ){
                say  ''    if  !$ps_ds;
                print ORFtrack $ps_top, $_    for  reverse( $past+1 .. $start);
                say  '';
            }

            if(  anyMarkup  ){
                {
                    for (reverse $past+1..$start){
                        my $i = $_;   $i += len  if $i < 0;
                        print   openTags( 1, $i ),  substr( $csq, $i, 1 ),  clozTags( 1, $i );
                    }
                    say  clozAllOpenTags 1;
                }
                if(  $ps_ds  ){
                    for (reverse $past+1..$start){
                        my $i = $_;   $i += len  if $i < 0;
                        print   openTags( 0, $i ),  substr( $seq, $i, 1 ),  clozTags( 0, $i );
                    }
                    say  clozAllOpenTags 0;
                }
            }
            else{  #No markup.
                if(  $past < -1  ){
                    my $end  =  len + $past+1;
                    if(  $start >= 0  ){
                        say  scalar reverse  substr($csq, $end) . substr $csq, 0, $start+1;
                        say  scalar reverse  substr($csq, $end) . substr $seq, 0, $start+1    if $ps_ds;
                    }else{
                        say  scalar reverse  substr $csq, $end, $start-$past;
                        say  scalar reverse  substr $seq, $end, $start-$past    if $ps_ds;
                    }
                }else{
                    say  scalar reverse  substr $csq, $past+1, $start-$past;
                    say  scalar reverse  substr $seq, $past+1, $start-$past    if $ps_ds;
                }
            }

            if(  $ps_ds  ){
                if(  ORFtrackAny $ps_bot, $past+1, $start  ){
                    print ORFtrack $ps_bot, $_    for  reverse( $past+1 .. $start);
                    say  '';
                }
                say  ''    if $lineBlockNum < $numLineBlocks;
            }

            $start -= $ps_numResPerLine;
        }#END:  for each line block.

    }#END:  if/else  $ps_dir == 1

}#END:  printSeq_multiline



#  Print $seq with inline ruler, wrapping lines at no more than $lineLen characters.
#  If $lineLen is long enough, sequences are printed as vertically aligned blocks of 10.
#  Otherwise lines of $lineLen are output, with no vertical alignment of blocks.
#  $ps_ds ∈ {0,1}  0: print as single stranded DNA,  1: print as double stranded DNA.
#  Assumes $ps_up ≠ $ps_dw.
sub printSeq_multiline_inlineRuler{
    goto &printSeq_multiline_inlineRuler_shortUnalignedBlocks   if $ps_numResPerLine < 30;
    goto &printSeq_multiline_inlineRuler_longUnalignedBlocks    if $ps_numResPerLine > 999;

    my $numBlocksPerLine  =   ceil  $ps_numResPerLine / 10;

    my $b5spc  =  $ps_numResPerLine > 99 ?  '  '  :  ' ';  #spacer after block 5,10,...
    my $b5len  =  length $b5spc;
    my $lineLen  =  11 * $numBlocksPerLine + 1 + floor(  ($numBlocksPerLine-1) / 5 ) * ($b5len - 1);

    my $numBasesPerLine  =  $ps_numResPerLine;

    my $L  =  len;  #To save time.

    #  Usually padding needs to be added to the start of the first line to make blocks of 10 bases align.
    my $leadBlockPadSize  =   $ps_top?  (10 - $ps_start%10) % 10  :  $ps_start%10;

    my $numBasesInPartialLine   =   min(  $numBasesPerLine - $leadBlockPadSize,   $ps_len  )   % $ps_numResPerLine;

    my $numWidth  =   max  6,  1+ length len;   #Widest coordinate is -len
    my $numFieldWidth  =   $numWidth + 2;
    my $leadBlockPadWithNumSize  =  $leadBlockPadSize + $numFieldWidth;

    #  $lineStartPos holds the circular position corresponding to first column in current line.
    #  In the initial partial line, this position be in the lead block.
    my $lineStartPos  =  $ps_start  -  $ps_dir * $leadBlockPadSize;
    my $linePastPos   =  $ps_start  +  $ps_dir * $numBasesInPartialLine;
    my $lineBegPos    =  $lineStartPos + $ps_top;
    my $lineEndPos    =  $linePastPos  + $ps_top;

    my $numSpacesOutput;  #Number of spaces printed so far on current line.
    my $pos  =  $ps_start;
    my $i;

    if(  $numBasesInPartialLine  ){
        # ─────  Print Initial Partial Line  ─────
        {
            if(   $ORFtrackAny[$ps_top]  &&  ORFtrackAny $ps_top, $ps_up, $lineEndPos   ){
                print  ' ' x $leadBlockPadWithNumSize;
                do{
                    print(  ($pos-$lineStartPos) %50 ?   ' '  :  $b5spc  )    if  $pos%10==0  &&  $pos != $lineStartPos;
                    $i  =  $pos < 0 ?  $L+$pos  :  $pos;
                    print ORFtrack $ps_top, $i;
                    $pos  +=  $ps_dir;
                }while  $pos != $linePastPos;

                say  '';
                $pos  =  $ps_start;
            }
            printf  "%${numWidth}d  ", $lineStartPos;
            print  ' ' x $leadBlockPadSize;   $numSpacesOutput = $leadBlockPadWithNumSize;
            do{
                if(  $pos%10==0  &&  $pos != $lineStartPos  ){
                    print_clozAllOpenTags $ps_top;
                    if(  ($pos-$lineStartPos) %50  ){  print ' ';     $numSpacesOutput++          }
                    else                            {  print $b5spc;  $numSpacesOutput += $b5len  }
                }
                $i  =  $pos < 0 ?  $L+$pos  :  $pos;
                print  openTags( $ps_top, $i ),  substr( $$ps_topSeq, $i, 1 ),  clozTags( $ps_top, $i );
                $pos  +=  $ps_dir;
            }while  $pos != $linePastPos;

            print_clozAllOpenTags $ps_top;
            print  ' ' x  ($lineLen - $numSpacesOutput - abs($pos-$ps_start));
            printf  "  %${numWidth}d\n",  $lineEndPos;
        }
        if(  $ps_ds  ){
            $pos  =  $ps_start;
            print  ' ' x $leadBlockPadWithNumSize;

            do{
                if(  $pos%10==0  &&  $pos != $lineStartPos  ){
                    print_clozAllOpenTags $ps_bot;
                    print(  ($pos-$lineStartPos) %50 ?   ' '  :  $b5spc  );
                }
                $i  =  $pos < 0 ?  $L+$pos  :  $pos;
                print  openTags( $ps_bot, $i ),  substr( $$ps_botSeq, $i, 1 ),  clozTags( $ps_bot, $i );
                $pos  +=  $ps_dir;   $i  =  $pos < 0 ?  $L+$pos  :  $pos;
            }while  $pos != $linePastPos;

            print_clozAllOpenTags $ps_bot;
            say  '';
            if(   $ORFtrackAny[$ps_bot]  &&  ORFtrackAny $ps_bot, $lineEndPos, $ps_up   ){
                print  ' ' x $leadBlockPadWithNumSize;
                $pos  =  $ps_start;
                do{
                    $i  =  $pos < 0 ?  $L+$pos  :  $pos;
                    print(  ($pos-$lineStartPos) %50 ?   ' '  :  $b5spc  )    if  $pos%10==0  &&  $pos != $lineStartPos;
                    print ORFtrack $ps_bot, $i;

                    $pos  +=  $ps_dir;
                }while  $pos != $linePastPos;

                say  '';
            }#END  if ORFtrackAny $ps_bot, $lineEndPos, $ps_up
            say  ''    if  $linePastPos != $ps_past;
        }#END  if $ps_ds
    }#END  print partial first line.

    return   if $linePastPos == $ps_past;          #<---  SUBROUTINE EXIT POINT


    #...                        sub printSeq_multiline_inlineRuler continues...

    my $numRemainingBases  =  $ps_len - $numBasesInPartialLine;
    my $numFullLines  =  int( $numRemainingBases / $numBasesPerLine );


    # ─────  Print all full lines (lines with $numBasesPerLine bases to print).
    for (1..$numFullLines){
        $lineStartPos  =  $pos;
        $lineBegPos  =  $lineStartPos + $ps_top;
        $lineEndPos  =   $pos  +  $ps_dir * $numBasesPerLine  +  $ps_top;
        {#  Forward strand.
            if(  $ORFtrackAny[$ps_top]   &&  ORFtrackAny $ps_top, $lineBegPos, $lineEndPos   ){
                print  ' ' x $numFieldWidth;
                for my $blockNum (0..$numBlocksPerLine-1){
                    print  $blockNum%5 ?  ' '  :  $b5spc    if $blockNum;
                    for (0..9){
                        $i  =  $pos < 0 ?  $L+$pos  :  $pos;
                        print ORFtrack $ps_top, $i;
                        $pos  +=  $ps_dir;
                    }
                }
                say  '';
                $pos  =  $lineStartPos;
            }

            printf  "%${numWidth}d  ", $lineStartPos;

            for my $blockNum (0..$numBlocksPerLine-1){
                print  clozAllOpenTags $ps_top,  $blockNum%5 ?  ' '  :  $b5spc   if $blockNum;
                for (0..9){
                    $i  =  $pos < 0 ?  $L+$pos  :  $pos;
                    print  openTags( $ps_top, $i ),  substr( $$ps_topSeq, $i, 1 ),  clozTags( $ps_top, $i );
                    $pos  +=  $ps_dir;
                }
            }
            print  clozAllOpenTags $ps_top;
            printf "  %${numWidth}d\n", $pos;
        }
        if(  $ps_ds  ){
            $pos  =  $lineStartPos;

            print  ' ' x  $numFieldWidth;
            for my $blockNum (0..$numBlocksPerLine-1){
                print  clozAllOpenTags $ps_bot,  $blockNum%5 ?  ' '  :  $b5spc   if $blockNum;
                for (0..9){
                    $i  =  $pos < 0 ?  $L+$pos  :  $pos;
                    print  openTags( $ps_bot, $i ),  substr( $$ps_botSeq, $i, 1 ),  clozTags( $ps_bot, $i );
                    $pos  +=  $ps_dir;
                }
            }
           say  '';
            if(   $ORFtrackAny[$ps_bot]   &&  ORFtrackAny $ps_bot, $lineEndPos, $lineBegPos   ){
                $pos  =  $lineStartPos;
                print  ' ' x $numFieldWidth;
                for my $blockNum (0..$numBlocksPerLine-1){
                    print  $blockNum%5 ?  ' '  :  $b5spc    if $blockNum;
                    for (0..9){
                        $i  =  $pos < 0 ?  $L+$pos  :  $pos;
                        print ORFtrack $ps_bot, $i;
                        $pos  +=  $ps_dir;
                    }
                }
                say  '';
            }
            say  ''    if  $pos != $ps_past;
        }#END  if $ps_ds
    }#END  for all full lines.

    return    if  $pos == $ps_past;                #<---  SUBROUTINE EXIT POINT


    # Print Final Partial Line. sub printSeq_multiline_inlineRuler continues...
    $lineStartPos  =  $pos;
    $lineBegPos  =  $pos + $ps_top;
    {#  Forward strand.
        if(  $ORFtrackAny[$ps_top]   &&  ORFtrackAny $ps_top, $lineBegPos, $ps_dw   ){
            print  ' ' x $numFieldWidth;
            do{
                print(   ($pos - $lineStartPos) %50 ?  ' '  :  $b5spc   )    if  $pos%10 == 0  &&  $pos != $lineStartPos;
                $i  =  $pos < 0 ?  $L+$pos  :  $pos;
                print  ORFtrack $ps_top, $i;
                $pos  +=  $ps_dir;
            }while(  $pos != $ps_past  );
            say  '';
            $pos  =  $lineStartPos;
        }
        printf  "%${numWidth}d  ", $lineStartPos;
        do{
            if(   $pos%10 == 0  &&  $pos != $lineStartPos   ){
                print  clozAllOpenTags $ps_top;
                print(   ($pos - $lineStartPos) %50 ?  ' '  :  $b5spc   );
            }
            $i  =  $pos < 0 ?  $L+$pos  :  $pos;
            print  openTags( $ps_top, $i ),  substr( $$ps_topSeq, $i, 1 ),  clozTags( $ps_top, $i );
            $pos  +=  $ps_dir;
        }while(  $pos != $ps_past  );
        print  clozAllOpenTags $ps_top;
        say  '';
    }
    if(  $ps_ds  ){
         $pos  =  $lineStartPos;
         print  ' ' x $numFieldWidth;
         do{
             if(   $pos%10 == 0  &&  $pos != $lineStartPos   ){
                 print  clozAllOpenTags $ps_bot;
                 print(   ($pos - $lineStartPos) %50 ?  ' '  :  $b5spc   );
             }
             $i  =  $pos < 0 ?  $L+$pos  :  $pos;
             print  openTags( $ps_bot, $i ),  substr( $$ps_botSeq, $i, 1 ),  clozTags( $ps_bot, $i );
             $pos  +=  $ps_dir;
         }while(  $pos != $ps_past  );
         print  clozAllOpenTags $ps_bot;
         say  '';
         if(  $ORFtrackAny[$ps_bot]   &&  ORFtrackAny $ps_bot, $ps_dw, $lineBegPos   ){
             $pos  =  $lineStartPos;
             print  ' ' x $numFieldWidth;
             do{
                 $i  =  $pos < 0 ?  $L+$pos  :  $pos;
                 print(   ($pos - $lineStartPos) %50 ?  ' '  :  $b5spc   )    if  $pos%10 == 0  &&  $pos != $lineStartPos;
                 print  ORFtrack $ps_bot, $i;
                 $pos  +=  $ps_dir;
             }while(  $pos != $ps_past  );
             say  '';
         }
    }

}#END  printSeq_multiline_inlineRuler()



#  ──────────────────────────────  Print Seq Functions   ──────────────────────────────
#  Print $seq with inline ruler, wrapping lines at no more than $lineLen characters.
#  $lineLen lines are output, with no vertical alignment of blocks.
#  According to $ps_ds,  0: print as single stranded DNA,  1: print as double stranded DNA.
#
#  This format is for line lengths so small that is does not look good,
#  and is only supported in the spirit of completeness.
sub printSeq_multiline_inlineRuler_shortUnalignedBlocks{

    my $lineLen  =  $ps_numResPerLine;  #local alias for readability.

    my( $line, @seqTop, @seqBot, @ORFtop, @ORFbot );

    my $numStr;

    if(  $ps_top  ||  $ps_up != 0  ){
        $numStr  =  ($ps_up<0? '':'+') . "$ps_up ";
        push @seqTop,  split( //, $numStr );
        push @seqBot,  (' ') x length $numStr   if $ps_ds;
        push @ORFtop,  (' ') x length $numStr;
        push @ORFbot,  (' ') x length $numStr   if $ps_ds;
    }

    #  Loop over each sequence position to be printed in the correct order to print them.
    while(   my( $iDéb,$iFin,$printAsNeg )  =  splice @ps_range, 0, 3   ){
        for my $i(  $iDéb..$iFin  ){
            $i *= -1   if $ps_top;    my $pos =  $i - len*$printAsNeg;

            $numStr = '';
            $numStr = $pos/10    if  !($pos%10)  &&  $pos != $ps_up  &&  $pos != $ps_dw;

            {
                #  Print top sequence.
                if(  $ORFtrackAny[$ps_top]   ){
                    push @ORFtop,   (' ')  x  (2 + length $numStr)    if  $numStr ne '';
                    push @ORFtop,  ORFtrack $ps_top, $i;
                    if(  @ORFtop >= $lineLen  ){
                        $_ =   join  '',  splice( @ORFtop, 0, $lineLen );   say    if /[a-zA-Z]/;
                    }
                }
                push @seqTop,   clozAllOpenTags($ps_top).' ',  split( //, $numStr),  ' '    if  $numStr ne '';
                push @seqTop,   openTags( $ps_top, $i )  .  substr( $$ps_topSeq, $i, 1 )  .  clozTags( $ps_top, $i );
                if(  @seqTop >= $lineLen  ){
                    say   join(  '',  splice @seqTop, 0, $lineLen  ),  clozAllOpenTags $ps_top;
                }
            }
            if(  $ps_ds  ){
                #  Print bottom sequence.
                my $atEndOfLine  =  0;
                push @seqBot,   clozAllOpenTags($ps_bot).' ',  (' ') x (1 + length $numStr)    if  $numStr ne '';
                push @seqBot,   openTags( $ps_bot, $i )  .  substr( $$ps_botSeq, $i, 1 )  .  clozTags( $ps_bot, $i );
                if(  @seqBot >= $lineLen  ){
                    say   join(  '',  splice @seqBot, 0, $lineLen  ),  clozAllOpenTags $ps_bot;
                    $atEndOfLine  =  1;
                }
                if(  $ORFtrackAny[$ps_bot]  ){
                    push @ORFbot,   (' ')  x  (2 + length $numStr)   if $numStr ne '';
                    push @ORFbot,  ORFtrack( $ps_bot, $i );
                    if(  @ORFbot >= $lineLen  ){
                        $line   =   join  '',  splice( @ORFbot, 0, $lineLen );
                        say $line    if  $line =~ /[a-zA-Z]/;
                    }
                }
                say  ''   if $atEndOfLine;
            }
        }#END  each $i
    }

    #  Add final position number.
    if(  !$ps_top  ||  $ps_dw != len  ){
        $numStr  =   ($ps_dw<0?  ' ' : ' +')  .  $ps_dw;
        push @seqTop,  clozAllOpenTags($ps_top),  split( //, $numStr );
        push @seqBot,  clozAllOpenTags($ps_bot),  (' ') x length $numStr   if $ps_ds;
        push @ORFtop,  (' ') x length $numStr;
        push @ORFbot,  (' ') x length $numStr   if $ps_ds;
    }

    #  Print partial line(s) at end.
    while(  @seqTop  ){
        $line   =   join  '',  splice( @ORFtop, 0, $lineLen );
        say $line    if  $line =~ /[a-zA-Z]/;
        say   join(  '',  splice @seqTop, 0, $lineLen  ),  clozAllOpenTags(  $ps_top)   if  @seqTop;
        say   join(  '',  splice @seqBot, 0, $lineLen  ),  clozAllOpenTags(1-$ps_top)   if  @seqBot;
        $line   =   join  '',  splice( @ORFbot, 0, $lineLen );
        say $line    if  $line =~ /[a-zA-Z]/;
        say  ''      if  $ps_ds  &&  @seqTop
    }

}#END:  printSeq_multiline_inlineRuler_shortUnalignedBlocks


#  Print $seq with binary megabase ruler, wrapping lines at $lineLen characters.
#  $lineLen lines are output, with no vertical alignment of blocks, but instead in horizontal blocks.
#  Because with such long lines, the font must be set too small to read for normal monitors,
#  printing in this format does not make much sense unless markup is used.
#  Currently, $ORFtrackFrameFwd,$ORFtrackFrameRev are ignored.
sub printSeq_multiline_inlineRuler_longUnalignedBlocks(){

    my $lineLen  =  $ps_numResPerLine;  #local alias for readability.

    my( $line, @seqTop, @seqBot, @ORFtop, @ORFbot );


    while(  my( $iDéb, $iFin, $printAsNeg )  =  splice @ps_range, 0, 3  ){
        for my $i ( $iDéb..$iFin ){
            $i *= -1   if $ps_top;    my $pos =  $i - len*$printAsNeg;

            if(  $pos % 100_000 == 0  ){
                if(  $pos % 1_000_000 == 0  ){
                    #  Flush seqTop, seqBot and print position in ascii art binary.
                    if(  @seqTop >= $lineLen  ){
                        {
                            $line  =  join(  '',  splice @ORFtop, 0, $lineLen  );
                            say   $ps_ds? "\n$_" : $_    if /[a-zA-Z]/;
                            say   join(  '',  splice @seqTop, 0, $lineLen  ),  clozAllOpenTags $ps_top;
                        }
                        if(   $ps_ds  ){
                            say   join(  '',  splice @seqBot, 0, $lineLen  ),  clozAllOpenTags $ps_bot;
                            $_ =  join(  '',  splice @ORFbot, 0, $lineLen  );   say    if /[a-zA-Z]/;
                            say  '';
                        }
                    }
                    my  $posInMB  =  int( abs $pos / 1_000_000 );
                    my  @posInMBbinary   =   split(  //, sprintf( '%012b', $posInMB )  );
                    say  '';
                    for (1..2){
                        print  $pos < 0?  ('Z'x50) : (' 'x50);
                        print  ' ' x 50;   #Print minus sign if needed.
                        print(  $_?  'Z' x 50 :  ' ' x 50   )    for  @posInMBbinary;
                        # print(  ($_ eq '1')? ('Z'x50) : (' 'x50)  )    for  @posInMBbinary;
                        say  '';
                    }
                    print(  ($_ eq '1')? ('Z'x50) : (' 'x50)  )    for  qw(0 0 1 1 1 1 0 0 0 0 1 1 1 1);
                    say  "\n\n";
                }
                else{#  $pos divisible by 100_000 but not 1_000_000
                    my @blockSep  =  $pos % 500_000
                        ?   (' ') x 10
                        :   split   //,   '   ' . sprintf( '%04d', $pos/1_000_000 ) . '   ';
                    {
                        push  @ORFtop,  (' ') x 10;
                        push  @seqTop,  @blockSep;
                    }
                    if(  $ps_ds  ){
                        push  @seqBot,  @blockSep;
                        push  @ORFbot,  (' ') x 10;
                    }
                }
            }

            {
                #  Print top sequence.
                if(  $ORFtrackAny[$ps_top]  ){
                    push  @ORFtop,  ORFtrack $ps_top, $i;
                    if(  @ORFtop >= $lineLen  ){
                        $_  =  join(  '',  splice @ORFtop, 0, $lineLen  );
                        say  $ps_ds?  $_ : "\n$_"    if /[a-zA-Z]/;
                    }
                }
                push  @seqTop,  openTags( $ps_top, $i) . substr( $$ps_topSeq, $i, 1 ) . clozTags( $ps_top, $i );
                if(  @seqTop >= $lineLen  ){
                    say   join(  '',  splice @seqTop, 0, $lineLen  ),  clozAllOpenTags $ps_top;
                }
            }
            if(  $ps_ds  ){
                #  Print bottom sequence.
                push  @seqBot,  openTags( $ps_bot, $i) . substr( $$ps_botSeq, $i, 1 ) . clozTags( $ps_bot, $i );
                my $atEndOfLine  =  0;
                if(  @seqBot >= $lineLen  ){
                    say   join(  '',  splice @seqBot, 0, $lineLen  ),  clozAllOpenTags $ps_bot;
                    $atEndOfLine  =  1;
                }
                if(  $ORFtrackAny[$ps_bot]  ){
                    push  @ORFbot,  ORFtrack $ps_bot, $i;
                    if(  @ORFbot >= $lineLen  ){
                        $_ =  join  '',  splice @ORFbot, 0, $lineLen;
                        say    if /[a-zA-Z]/;
                    }
                }
                say  ''  if $atEndOfLine;
            }

        }
    }#END: for each $i

    #  Print partial line(s) at end.
    while(  @seqTop  ){
        {
            $_ =  join(  '',  splice @ORFtop, 0, $lineLen  );   say   if /[a-zA-Z]/;
            say   join(  '',  splice @seqTop, 0, $lineLen  ),  clozAllOpenTags $ps_top;
        }
        if(  $ps_ds  ){
            say   join(  '',  splice @seqBot, 0, $lineLen  ),  clozAllOpenTags $ps_bot;
            $_ =  join(  '',  splice @ORFbot, 0, $lineLen  );   say   if /[a-zA-Z]/;
            say  '';
        }
    }

}#END:  printSeq_multiline_inlineRuler_longUnalignedBlocks



#  Check input filesnames passed as @_, dying with error if problem found.
#  Achtung!  Adds '-' to @_ when empty.
sub checkInputFilenames{

    my $numDashes  =  0;
    for my $filename (@_){
        if(  $filename eq '-'  ){
            _dieWithUsage  'STDIN given as input twice'                                 if $numDashes++;
            _dieWithUsage  'STDIN stipulated for input, but is attached to a terminal'  if -t *STDIN;
        }
        else{  # filename ne '-'
            _looksLikeFastaFile( $filename )   or  _dieWithUsage  "file '$filename' is not a fasta file";
        }
    }

}#END: checkInputFilenames( filenames... )



#  Returns runs of $runLen or more occurrences of the same character as FastaplSeqsegs object.
sub runs($){
    @_ == 1   or   die  'lcRuns() expected one argument but got', 0+@_, ' args.';
    my $runLen  =  $_[0];

    $runLen > 1   or   die  "lcRuns() expected runLen to be a positive integer, but got '$runLen'";

    my $retVal  =  FastaplSeqsegs->new( 'runs' );

    my $i  =  0;
    my $j  =  $i + 1;
    my $d  =  lc substr $seq, $i, 1;
    my $len  =  len;

    for(  ;  $j < $len;  $i = $j++  ){
        my $c  =  $d;
        ++$j    while   $j < $len   &&   $c  eq  ($d = lc substr $seq, $j, 1);
        #  $j should now be the smallest index > $i, such that seq $j ≠ seq $i.
        $retVal->push( $i,$j )    if  $j-$i > $runLen;
    }
    return $retVal;
}#END: runs($)



#  Return hash holding any KEY=VALUE pairs present in @_
#  @_ defaults to @head.
sub asHash{

    @_   or   @_  =  @head;

    my @keyValuePairs;

    #                             KEY     =    VALUE
    push  @keyValuePairs,  /^  ( [^=]+ )  =  ( [^=]+ )  $/x    for  @_;

    return @keyValuePairs;
} # END: asHash()



#  kmers( $k, [$regex] )  Returns associative list:  kmer₁, freq₁, ...
#  of the fixed lengthed substrings (e.g. k-mers) in $seq and their frequencies.
#
#  Optionally $regex can be used to select which k-mers to include.
#  So for example [acgt] would exclude k-mers with n's or upcased ACGT
sub kmers{
    defined wantarray   or   die  'useless call of kmers in void context';

    #  else return kmers as list.
    0 < @_ and @_ < 3   or   die  'kmers expected one or two arguments';   #perl v5.32 could use 1 < @_ < 3
    my %freq;
    my $k = shift;

    $k  =~  /^[1-9]\d*$/   or   die  "kmers expected positive integer argument but got $k";

    my $regex= shift;

    if(  $regex  ){
        #  Reality check regex.  Once I passed $seq by mistake when it held a chromosome.
        my $regexLen= length $regex;
        $regexLen < 999   or   die  "kmers passed regex of length $regexLen, seems too long to be want you want.";

        for  my $i  (0..fin+1-$k){
            my $kmer =  substr $seq, $i, $k;
            ++$freq{$kmer}   if  $kmer=~m/$regex/;
        }
    }
    else{
        ++$freq{substr $seq, $_, $k}   for  0..fin+1-$k;
    }

    _kmerFreqs_add_wrapped_fragment( \%freq, $k , $regex )   if $seqIsCirc;

    return  %freq;
}#END: kmers()


#  _kmerFreqs_add_wrapped_fragment( \%freq, $k, [$regex] )
#  For circular DNA molecules.
#  Add to counts of k-mers which straddle the end of $seq, to the frequencies in hash %freq
sub _kmerFreqs_add_wrapped_fragment{
    my( $freqRf, $k, $regex )=  @_;

    if( $seqIsCirc ){
        my $wrapped_frag=   substr ($seq, len+1-$k)  .  substr ($seq, 0, $k-1);
        for  my $i  (0..$k-2){
            my $kmer=  substr $wrapped_frag, $i, $k;
            ++$freqRf->{$kmer}   if !$regex or $kmer=~m/$regex/;
        }
    }
}


#  addKmers( %freq, $k, [$regex] ) adds counts of length $k kmers in $seq to %h
#  if $regex given, only kmers matching it are added.
#  addKmers( %freq, $k )  is equivalent to addTo( %h, kmers $k ), but saves memory (0~40%) for large $k.
#
#  This function is for use in fastapl user scripts.
#  Call syntax is addKmers( %freq, $k )
#  Note the prototype \%@ here.
sub addKmers( \%@ ){
    die  'call of addKmers in non-void context'   if  defined wantarray;
    1 < @_ and @_ < 4   or   die  'addKmers expected two or three arguments';   #v5.32 could use 1 < @_ < 4
    my( $freqRf, $k )=  @_[0,1];

    ref $freqRf eq 'HASH'   or   die  'addKmers expected first argument to be a hash reference';

    $k  =~  /^[1-9]\d*$/   or   die  "addKmers expected positive integer argument but got $k";

    my $regex=  qr/$_[2]/  if defined $_[2];

    if(  $regex  ){
        my $kmer;
        for  my $i  (0..fin+1-$k){
            $kmer =  substr $seq, $i, $k;
            ++$freqRf->{$kmer}   if  $kmer=~m/$regex/;
        }
    }
    else{
            ++$freqRf->{substr $seq, $_, $k}   for 0..fin+1-$k;
    }

    _kmerFreqs_add_wrapped_fragment( $freqRf, $k , $regex )   if $seqIsCirc;

}#END: addKmers(\%@)



#  Print input parsing error and exit.
sub sayInputError{
    say2  "fastapl; Input parsing error; $_[0]";
    exit -1;
}


#  open $opt_{outputFilename} and redirect STDOUT to it.
sub redirectSTDOUT{
    -t STDOUT
        or   _dieWithUsage  "not sure whether to output to file '$opt_{outputFilename}' or STDOUT";

    my $outStream  =  *STDOUT;   #STDOUT is now redirected to $outStream;
    open(   $outStream, '>', $opt_{outputFilename}   )
        or   _dieWithUsage  "could not open output file '$opt_{outputFilename}', $!";
}


#  Return true iff filename $_[0] appears to be a fasta file.
#  Determined heuristically by looking for a line starting with '>'.
sub _looksLikeFastaFile{
    open  my $file, '<', $_[0]
        or  _dieWithUsage  "could not open fasta file '$_[0]', $!";

    while(  <$file>  ){
        return 1    if  /^>/;
        return 0    if  /^[^#]/;
    }
    return 0;
}


sub widen{
    my $self = pop;

    return  $self->widen( @_ );
}




1;
