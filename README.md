# fastapl -- Manipulate sequence data from the command line in a perl one-liner style.


## Getting started
`% export PERL5LIB=$PERL5LIB:THIS_DIRECTORY/modules/perl5`
`% fastapl5 --help`

## Running regression tests
`% perl t/runRegressionTests.pl


## Summary
fastapl facilitates manipulation of fasta format sequence data via snippets of perl code.  The sister program fastqpl is similar, but for fastq format data.  fastapl is basically just wraps user provided code to execute in Perl, but it also provides some convenient functions such as reverse complementation (rc) and translation from nucleotide to amino acid sequence (translate, translations).  Although fastapl is a command line tool, it can be used to mark sequence regions of interest in bold or color in an xterm terminal or in html format.  For example, color red, matches "CACA=1" displays matches to CACA in red.  fastapl functions like these are designed with the understanding that DNA has two strands and sometimes is circular, so that those functions flexibly handle various use cases.

### Examples

#### fastapl Examples
* Truncate sequences to length 39.<BR>
`% fastapl5 -pe 'trim 0,38'  proteins.fasta`

* Reverse complement DNA sequences.<BR>
`% fastapl5 -p  -e '$seq = reverse $seq; $seq =~ tr/acgtACGT/tgcaTGCA/'  shortDNA.fasta`

* Compute fraction of sequences matching regex L$.<BR>
`% fastapl5 -e '++$tot; ++$m if( $seq =~ /L$/ );'  -f 'say "$m/$tot"'  proteins.fasta`

* Print records of sequences not starting with 'M'.<BR>
`% fastapl5 -ge '$seq !~ /^M/'  proteins.fasta`

* Compute total sequence length from two input files.<BR>
`% fastapl5 -e '$tot += length $seq'  --end 'print "total seq. length = $tot\n"' proteins1.fasta proteins2.fasta`

* Compute residue composition of collection of sequences.<BR>
`% fastapl5 -e '++$tot, $a{$_}+=100 for @seq'  --end 'printf "%s: %5.2f\n", $_, $a{$_}/$tot  for sort keys %a'	proteins.fasta`

* Same but skip initial residue of each sequence.<BR>
`% fastapl5 -e '++$tot, $a{$_}+=100 for @seq[1..fin]'  --end 'printf "%s: %5.2f\n", $_, $a{$_}/$tot  for sort keys %a'  proteins.fasta`

* Compute residue composition of each sequence.<BR>
`% fastapl5 -e '%a=(); $a{$_}+=100 for @seq;  print $id;  printf " $_:%4.2f", $a{$_}/len  for sort keys %a;  say ""'  proteins.fasta`

* Randomly shuffle sequences.<BR>
`% fastapl5 -b 'srand 7'  -pe '$seq = join "", shuffle @seq'  proteins.fasta`

* Trim poly-A tail allowing for some miscalled bases.<BR>
`% fastapl5 -p  -e '$s = 0; $p = 0; $m = -9; for $i (reverse 0..$#seq){ $s += $seq[$i] eq "a" ? 1 : -2; if( $s > $m ){ $m = $s; $p = $i }} $m < 8 or substr( $seq, $p ) = ""'  DNAwithPolyAtail.fasta`

* Sort records by id.<BR>
`% fastapl5 --sort -e '$id1 cmp $id2'  proteins.fasta`

* Sort records by sequence length.<BR>
`% fastapl5 -se 'len1 <=> len2'  proteins.fasta`

* Reformat records so that sequence line lengths are 100 (at most).<BR>
`% fastapl5 -pl 100  proteins.fasta`

* Print first 3 records.<BR>
`% fastapl5 -ge '$c++ < 3'  proteins.fasta`

* Select records matching given ids.<BR>
`% fastapl5 -b '%keep = qw(Q12154|GET3_YEAST 1 P24583|KPC1_YEAST 1 P38179|ALG3_YEAST 1)' -ge '$keep{$id}'  proteins.fasta`

* Select records matching an id in ids.txt, ids.txt holds one id per line.<BR>
`% fastapl5 -b '%keep = map {chop;$_,1} <STDIN>' -ge '$keep{$id}'  proteins.fasta  <  ids.txt`

* Reformat records so that record sequences are printed on one line.<BR>
`% fastapl5 -p -1  proteins.fasta`

* Extract records with first field after ID matching /mito/.<BR>
`% fastapl5 -ge '$head[1] =~ /mito/'  proteins.fasta`

* Print forward and reverse complementary strands with appropriate change to ids.<BR>
`% fastapl5 -pe '$id .= "_Crick"; pr; rc; $id =~ s/Crick$/Watson/'  shortDNA.fasta`

* Convert DNA sense strand to RNA transcript sequence.<BR>
`% fastapl5 -pe '$seq =~ tr/tT/uU/'  shortDNA.fasta`

* Dump 6 possible frame translations of DNA sequence.<BR>
`% fastapl5 -e 'say $id;  say "  $_"  for translations'  DNAwithPolyAtail.fasta`

* Add comments listing ORF fragments of length 29 or more to DNA sequence records.<BR>
`% fastapl5 -pe '$comment .= sprintf "# %2d %s\n", @$_[2,3]  for ORFs 29'  DNAwithPolyAtail.fasta`

* Add comments listing ORFs of length 29 or more to DNA sequence records.<BR>
`% fastapl5 -pe '$comment .= sprintf "# %2d %s\n", @$_[2,3]  for ORFs "m29"'  DNAwithPolyAtail.fasta`

* Upcase sequences.<BR>
`% fastapl5 -pe '$seq = uc;'  DNAwithPolyAtail.fasta`

* Print 2-mer frequences for each seq.<BR>
`% fastapl5 -e 'say hstr kmers 2'  mitoGenomes.fa`

* Print 2-mer frequences for each seq in compact format.<BR>
`% fastapl5 -e 'say  "$id\t",  hstrf ": ", kmers 2'  mitoGenomes.fa`

* Print cumulative 2-mer frequences of sequences, one 2-mer per line.<BR>
`% fastapl5 -e 'addKmers %c, 2' -f 'say hstr %c'  mitoGenomes.fa`

* Print table of homopolymer runs of length 20 or more.<BR>
`% fastapl5 -e '$m = runs 20;  say  join "\t",  $a, $b, $b-$a, seq $a  while $m->ab'  yeast_chrII.fa`

* Color matches to 'CACA' on either strand red, also underline the matches on the reverse strand; show all of these on the forward strand.<BR>
-H -pe 'color red, matches "CACA=1"; ul matches "CACA_1"'  shortDNA.fasta

* Print sequences in format of 100 bases per line.<BR>
`% fastapl5 -pl 100  proteins.fasta`

* Count the number of records with first field 'E.R.'.<BR>
`% fastapl5 -ce '$head[1] eq "E.R."'  proteins.fasta`

* Print information encode in header as key=value pairs.<BR>
`% fastapl5 -e 'say "  $id"; say hstr asHash'  proteinsHashHeaders.fasta`

* Print aromatic residues in bold with rule.<BR>
`% fastapl5 -r -pe 'bold matches "[YWF]"'  proteins.fasta`

* Print sequence as double stranded DNA.<BR>
`% fastapl5 -P  shortDNA.fasta`

* Colorize final four residues in xterm format.<BR>
`% fastapl5 -pe 'colorize -4,0'  proteins.fasta`

* Print sequence as double stranded with inline ruler and ORF colored ORF track in xterm format.<BR>
`% fastapl5 --ORF a110m80c  -Pr  -l100  InfluenzaA_H1N1_Iwate200903.fa`

* Print bases 100 bases starting at position 5000 (counting from zero) in each record with ruler track.<BR>
`% fastapl5 -Re 'say "$id 5000,5100"; ps 5000,5100'  mitoGenomes.fa`

* Print entries with residues colorized.<BR>
`% fastapl5 -C  proteins.fasta`

* Print entries with residues colorized as DNA.<BR>
`% fastapl5 -Cd  shortDNA.fasta`

* Dump final field of each head, treating header line as semicolon separated fields.<BR>
`% fastapl5 '-F;' -e 'say $head[-1]'  semicolonHeadFields.fasta`

* For circular genome, print the reverse strand sequence seqment of length 10 centered on the origin.<BR>
`% fastapl5 -Oe 'say seg 5, -5'  mitoGenomes.fa`

* Mark TATAT on forward strand red, allowing matches to overlap.<BR>
`% fastapl5 -pe 'color red, matches "TATAT%"'  m.fa`

* Mark matches to 'GAATAT' on either strand in inverse video, for circular DNA.<BR>
`% fastapl5 -OPe 'iv matches "GAATAT="'  mitoGenomes.fa`

* Downcase $seq, then upcase matches to 'tata[at]a'.  Note the -d flag.<BR>
`% fastapl5 -dpe '$seq = lc; upcase matches "tatawa"'  yeast_chrII.fa`

* Print regex to case insensitively match an adenine followed by six bases coding for amino acids sequence 'MA'.<BR>
`% fastapl5 -db 'say expandRegex "(?i)a<MA>"'`

* Count the number of start codons.  Works for RNA or DNA sequences, upper or lower case.<BR>
`% fastapl5 -e '$c=()=/(?i)(A[UT]G)/g; say "$id $c"'  yeast_chrII.fa`


#### fastqpl Examples
* Print records with id containing substring "4401B".<BR>
`fastqpl5 -g  -e  '$id =~ /4401B/'  longreads_original_sanger.fastq`

* For each record: print (id, p), where p is the probability that all bases are correct.<BR>
`fastqpl5 -e '$p = 1;  $p *= (1 - $_) for( errProb() ); print "$p\n"'  longreads_original_sanger.fastq`

* Extract fasta file from fastq input.<BR>
`fastqpl5 -e 'print ">$head\n$seq\n"'  longreads_original_sanger.fastq`

* Sort records in order of descending average Phred score.<BR>
`fastqpl5 -se '(sum Phred1) / length $seq1  <=>  (sum Phred2) / length $seq2'  longreads_original_sanger.fastq`

* Trim off first 4 bases of each read.<BR>
`fastqpl5 -pe 'trim 4'   longreads_original_sanger.fastq`

* Another way to trim off final 4 bases of each read.<BR>
`fastqpl5 -pe 'trim 0,-5'   longreads_original_sanger.fastq`

* Dump sequence with Phred scores.  Output looks like: 'FSRRS4401BE7HA  t:37 c:37 a:37 g:35 T:35 T:35...'.<BR>
`fastqpl5 -e '@p = Phred; @v = pairwise {"$a:$b"} @seq, @p; say "$id  @v"'   longreads_original_sanger.fastq`

* Trim off any bases with error prob > 0.1 at end of read.  Discard record if none found.<BR>
`fastqpl5 -ge  'trim  0,  lastidx {$_ < 0.1} errProb;  $seq'   longreads_original_sanger.fastq`

* Add fake poly-A tail to each read.  Each added base having an increasing error prob.<BR>
`fastqpl5 -b 'srand 5'  -pe  '$len = rand 30; $p = .05; @p = map {$p *= 1.1} 1..$len; $seq .= "a" x $len; $quality .= prob_to_qualStr @p'  longreads_original_sanger.fastq`

* Keep the (shortest) prefix with the highest total adjusted (15 subtracted from each score) Phred scores.  Similar to BWA trimming.<BR>
`fastqpl5 -ge '$sum = 0;  ($i,$max) = maxidx apply {$_ = $sum += $_-15} Phred;  trim 0,$i;  $max > 0'  longreads_with_dummy.fastq`

* Keep the segment with the highest total adjusted (15 subtracted from each score) Phred scores.  Similar to Solexa trimming.<BR>
`fastqpl5 -ge '($beg,$end,$max) = maxSeg apply {$_ -= 15} Phred;  trim $beg,$end;  $max'  longreads_with_dummy.fastq`



### Conversion to fasta format
fastapl only accepts fasta format files.  If you have a text file of sequences, one per line, you can convert them to fasta format with perl like this:<br>
`% perl -pE 'say ">seq$."' seqs.text > seqs.fasta`

### Dependencies
fastapl, fastqpl depend on several standard modules which may or may not already be installed on your system.  For example `List::MoreUtils` and `indirect`.  If you get an immediate error message when you first try to run fastapl, please read the error message to see if it is complaining about not being able to find a module.

### Caveats.
Runs on a linux box under several version of perl (e.g. v5.18.2).  Not tested elsewhere.


### See Also
perltab


### Author
Paul Horton.  Copyright 2010,...,2017.
