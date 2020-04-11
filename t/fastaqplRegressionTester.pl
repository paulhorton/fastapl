#!/usr/bin/perl -w
#  Author: Paul Horton
#  Copyright (C) 2011, 2017 Paul Horton, All rights reserved.
#  Created: 20110108
#  Updated: 20170805
#
#  regression tester for fastapl or fastqpl one-liners
#
#  See pod below for details.

use strict;
use Getopt::Long qw(:config posix_default bundling permute);
use Pod::Usage;
use File::Basename;
use File::Spec;
use List::MoreUtils 'any';
use Cwd;
use Time::HiRes;
use feature qw(say state);


sub nextLine($);    #Return next line from $file, skipping empty and comment lines.
sub printError($);  #Print $message to STDERR.  And to STDOUT as well if it is bound to a terminal.


#  ──────────  Declarations of ARGV related global variables.  ──────────
my $fastaqplPath;  #Pathname to fastapl or fastqpl file to be tested.

my $oneLinerFile    =  undef;  #File holding one-liners to test.
my $oneLinerFilename    = '';  #Name of file holding one-liners to test.
my $oneLinerDirPathname = '';  #pathname of directory holding $oneLinerFile.

my @oneLinersToRun        =  ();  #To hold one-liners to run, if stipulated on command line.

my %cksumDone;  #$cksumDone{FILE} TRUE iff cksum done on FILE when reading the one-liner header.

my $generatePrograms_flag=   0;  #If given, also test programs generated with fastapl -n option.
my $printExpected_flag=      0;  #Print expected output file pathname.
my $opt_abortOnErrorP=   undef;  #Abort immediately upon error? 
my $absolutePathnames_flag=  0;  #Print pathnames as absolute pathnames?
my $justPrint_flag        =  0;  #Print but do not execute commands?
my $verbose_flag          =  0;  #Print commands before executing?


#  ──────────  Declarations of other global variables.  ──────────
my $countString = '';        #To hold string formatted version of current command number.
my $expectedOutputPathname;  #To hold pathname of expected output.
my $curCommandGrepish;       #Current command is grepish {-g|-c|-q}?

$| = 1;  #Do not buffer output.


# ━━━━━━━━━━━━━━━  Parse Command Line  ━━━━━━━━━━━━━━━

{   # ────────────  BEG: Process command line options  ────────────
    my %docFlag = ();
    my @commandLineOptionSpecs = (
        'help'                =>  \$docFlag{help},
        'man'                 =>  \$docFlag{man},
        'options'             =>  \$docFlag{options},
        'usage'               =>  \$docFlag{usage},
        'a|abortOnErrorP'     =>  \$opt_abortOnErrorP,
        'A|absolute-pathnames'=>  \$absolutePathnames_flag,
        'g|generate-programs' =>  \$generatePrograms_flag,
        'e|print-expected'    =>  \$printExpected_flag,
        'n|just-print'        =>  \$justPrint_flag,
        'v|verbose'           =>  \$verbose_flag,
        );

    my $getOptionsRetval  =  GetOptions( @commandLineOptionSpecs );

    $getOptionsRetval  or  pod2usage();

    $verbose_flag |= $justPrint_flag;

    # ──────────  Process flags related to usage or documentation  ──────────
    my $isDumbTerminal   =   exists($ENV{TERM})  &&  $ENV{TERM} eq 'dumb';

    $docFlag{man}    and  pod2usage( -verbose   => 2,
                                     -noperldoc => $isDumbTerminal );
    $docFlag{help}   and  pod2usage( -verbose => 1 );
    $docFlag{usage}  and  pod2usage( -verbose => 0 );


    my @optionSummary = (
        'a|abortOnErrorP       Abort immediately upon error',
        'A|absolute-pathnames  Print pathnames as absolute pathnames',
        'd|direct-only         Only run fastapql directly, do not try standalone',
        'e|print-expected      Print expected output file pathname',
        'n|just-print          Do not run commands, just print them'
        );         

    if(  $docFlag{options}  ){
        say " $_"  for @optionSummary;
        exit 1;
    }

}   # ────────────  END: Process command line options  ────────────



{ # ───────────────  BEG: Check arguments  ───────────────
    
    printError  'not enough arguments'  if @ARGV < 2;
    printError  'too many arguments'    if @ARGV > 3;

    my $curArg  =  shift @ARGV;

    if(   my ($beg,$end)  =  $curArg =~  /^  (\d+)  (?: \.\. (\d+)? )?  $/x   ){
        defined $end   or   $end  =  $curArg =~ /\.\./?  99  : $beg;
        ($beg < 100)  or  printError  'maximum one-liner number is 99, but got "$beg"';
        ($end < 100)  or  printError  'maximum one-liner number is 99, but got "$end"';
        @oneLinersToRun  =  ($beg .. $end);
        $curArg  =  shift @ARGV;
    }

    $fastaqplPath  =  $curArg;

    $fastaqplPath =~ /fast[aq]pl[56]?$/
        or   printError  "2nd argument should give a pathname ending in 'fastapl' or 'fastqpl', but got '$curArg'";

    -e $fastaqplPath   or   die  "pathname '$fastaqplPath' given for fastaqpl program does not exist.\n";
    convertToAbsPath( $fastaqplPath );

    $curArg  =  shift @ARGV;
    printError  'too many arguments?'  if @ARGV;
    $oneLinerFilename  =  $curArg;

    printError  "one-liner file '$oneLinerFilename' should not be a directory"  if -d $oneLinerFilename;

    ( open  $oneLinerFile, '<', $oneLinerFilename )
        or  printError  "could not open file '$oneLinerFilename', $!";

    $oneLinerDirPathname  =  dirname $oneLinerFilename;
    chdir  $oneLinerDirPathname;

} # ───────────────  END: Check arguments  ───────────────


# ───────────────  Define some vars  ───────────────
my $count=  0;
my $oneLinerNumber;

my $oneLinerNumberLine;
my $oneLinerCoreLine;


# ━━━━━━━━━━━━━━━━━━━━  Parse one-liner file header  ━━━━━━━━━━━━━━━━━━━━

# ──────────  Set input, expectedOutput partial pathnames  ──────────
my $inputDir         =  '';  #To prepend to  input pathnames as given in one-liners file
my $expectedOutputDir=  '';  #To prepend to output pathnames as given in one-liners file
{  
    my $curLine  =  nextLine $oneLinerFile;

    ($inputDir)  =  $curLine =~ /\$INPUT_DIR\s*=\s*(\S+)\s*$/
        or   die  "Error; when reading one-liner file '$oneLinerFilename', expected '\$INPUT_DIR =' assignment line but got: '$curLine'\n";

    $curLine  =  nextLine $oneLinerFile;

    ($expectedOutputDir)  =  $curLine =~ /\$OUTPUT_DIR\s*=\s*(\S+)\s*$/
        or   die  "Error; when reading one-liner file '$oneLinerFilename', expected '\$OUTPUT_DIR =' assignment line but got: '$curLine'\n";

    -d $expectedOutputDir  or
        die  "Error; when reading one-liner file '$oneLinerFilename', expected '$expectedOutputDir' to be a directory, but it is not\n";
    
# ──────────  Check inputfile cksums  ──────────
    $curLine =  nextLine $oneLinerFile;
    $curLine =~ /^BEG_CKSUMS$/   or   die  'could not find CKSUMS section';
    
    while(   $curLine=  nextLine $oneLinerFile   ){
        last   if $curLine =~ /^END_CKSUMS$/;

        # Remove any leading white space and the trailing newline from $curLine.
        chop $curLine;    $curLine  =~  s/^\s*//;

        my ($claimedCksum, $claimedFilesize, $pathname)  =  split / +/, $curLine;

        defined $pathname   or   die  "when parsing cksum lines, could not parse this line: '$curLine'\n";
        $pathname  =  "$inputDir/$pathname";

        die  "pathname: '$pathname' listed twice in one-liner file header, line '$curLine'\n"   if $cksumDone{$pathname};

        -e $pathname   or  die  "file '$pathname' does not exist.  It was listed in one-liner file line: '$curLine'\n";

        my $cksumResult  =  `cksum $pathname`;
        (my $cksum, my $filesize)  =  split / +/, $cksumResult;
        $filesize  eq  $claimedFilesize
            or   die  "file: '$pathname' filesize ($filesize) differs from claimed filesize ($claimedFilesize)\n";
        $cksum  eq  $claimedCksum
            or   die  "file: '$pathname' cksum ($cksum) differs from claimed cksum ($claimedCksum)\n";

        $cksumDone{$pathname}  =  1;
    }
}



while( <$oneLinerFile> ){  #Loop through one-liners file lines
    chop;

    #  ──────────  Process non-core lines  ──────────
    next  if /^$/;
    next  if /^#/;
    if(  /^(\d\d)/  ){
        $oneLinerNumber      =  $1;
        $oneLinerNumberLine  =  $_;
        next;
    }


    #  ━━━━━━━━━━━━━━━━━━━━━━━━━  Process core line  ━━━━━━━━━━━━━━━━━━━━━━━━━

    #  ────────────  Make assignment from core line info  ────────────
    $oneLinerCoreLine  =  $_;
    {
        my $tabCount   =   $oneLinerCoreLine =~ tr/\011/\011/;
        $tabCount > 0
            or   die  "Error; when reading list of one-liners in '$oneLinerFilename': expected at least one tab in line: '$oneLinerCoreLine'\n";
        $tabCount < 3
            or   die  "too many tabs ($tabCount) in line: '$oneLinerCoreLine'\n"
    }

    $oneLinerNumber   or   die  'one-liner number not found';


    my( $commandArgs, $inputFilePathnames, $stdinRedirectPartialPathname )  =  split /\t/, $oneLinerCoreLine;

    $curCommandGrepish  =  $commandArgs =~ /-q|-c|-g/;

    $inputFilePathnames   or   die  "Could not find input filename(s) in line: '$oneLinerCoreLine'";

    if(  $inputFilePathnames  ne  'dummy'  ){
        my @inputFilePathname   =   map {"$inputDir/$_"}   split / +/, $inputFilePathnames;
        $cksumDone{$_}   ||   die  "No cksum declared for input file '$_'\n"    for @inputFilePathname;
        if(  $absolutePathnames_flag  ){
            convertToAbsPath($_)   for @inputFilePathname;
        }
        $inputFilePathnames   =   join  ' ',  @inputFilePathname;
    }else{
        $inputFilePathnames   =   '';
    }

    my $stdinRedirectPathname = '';
    if(  $stdinRedirectPartialPathname  ){
        $stdinRedirectPathname  =  "$inputDir/$stdinRedirectPartialPathname";
        $cksumDone{$stdinRedirectPathname}   ||   die  "No cksum declared for stdin redirect file '$stdinRedirectPathname'\n";
    }

    ++$count;

    #  If command line says to only execute one command, skip over any others.
    if(  @oneLinersToRun  ){
        next   unless  any  {$count == $_}  @oneLinersToRun;    #Next one-liner line   <───  LOOP EXIT
    }

    $oneLinerNumber == $count
        or   die  "$0 Inconsistency: command count = $count, but command number given in one liner file is '$oneLinerNumber'\n";

    $countString  =  sprintf "%02d", $count;

    $expectedOutputPathname  =  "$expectedOutputDir/out${countString}";
    convertToAbsPath( $expectedOutputPathname )   if $absolutePathnames_flag;

    ( -e $expectedOutputPathname )
        or   die  "Error; could not find expected output file '$expectedOutputPathname'\n";


    {   # ─────  Ensure consistency with expected output file header  ─────
        chop(   my $expectedOutputLine1  =  `head -1 $expectedOutputPathname`   );

        if(  $expectedOutputLine1 ne $oneLinerNumberLine ){
            warn  "In one-liners file number one-liner description line was:\n    '$oneLinerNumberLine'\n";
            warn  "but the first line in the expected output file is:\n    '$expectedOutputLine1'\n";
            die  'aborting...'    if $opt_abortOnErrorP;
        }

        chop(   my $expectedOutputLine2  =  `head -2 $expectedOutputPathname | tail -n +2`   );

        if(  $expectedOutputLine2 ne $oneLinerCoreLine  ){
            warn  "In one-liners file core line was:\n    '$oneLinerCoreLine'\n";
            warn  "but in expected output file '$expectedOutputPathname' was:\n    '$expectedOutputLine2'\n";
            die  'aborting...'    if $opt_abortOnErrorP;
        }
    }


    #  ─────  Handle case in which -e, -n, or -v flags are given.
    if(  $printExpected_flag  ){
        say $expectedOutputPathname;
    }

    if(  $verbose_flag  ){
        print  "\n$countString: ";
        print  "$fastaqplPath "    if $justPrint_flag;
        print  "$commandArgs  $inputFilePathnames";
        print  " < $stdinRedirectPathname"   if $stdinRedirectPartialPathname;
        print  "\n";
    }

    next  if $justPrint_flag;  #Early jump to next one-liner  <---------- JUMP


    #  ───────────────  Test fastapql executed directly  ───────────────
    my $command  =  "$fastaqplPath $commandArgs $inputFilePathnames";
    say  "$countString: " . $command;

    $command .= " < $stdinRedirectPathname"   if $stdinRedirectPartialPathname;


    print  'direct test...';
    assertCommandOutputMatchesExpected( $command );

    unless(  $generatePrograms_flag  ){
        say '';  next;    #Go to next one-liner
    }

    
    #  ───────────────  Test program made with '-n' flag  ───────────────
    $command  =  "$fastaqplPath -n $commandArgs";

    my $prog  =  qx( $command );

    print  '   Testing via -n flag...';
    my $progFilename  =  "/tmp/prog${countString}.pl";
    my $progFile;
    open  $progFile, '>', "$progFilename";
    print  $progFile $prog;
    close  $progFile;

    $command  =  "perl $progFilename $inputFilePathnames";

    $command .= " < $stdinRedirectPathname"  if $stdinRedirectPartialPathname;

    assertCommandOutputMatchesExpected( $command );
    say '';

} #  END: while( <$oneLinerFile> )



#  Execute $command and store its results in /tmp directory.
#  Then check that the results are the same as expected by
#  comparing to the file $expectedOutputPathname.
#  Die if the result differs.
sub assertCommandOutputMatchesExpected{
    state  $cwmTime = 0;
    my ($command)  =  @_;

    my $timeBefore  =  Time::HiRes::time;

    my $output  =  qx( $command );
    $output  or   die  'No output generated';

    $cwmTime  +=  Time::HiRes::time - $timeBefore;
    printf  "%5.2fs ", $cwmTime;


    $?   and  die  "Non-zero status '$?' returned"   if !$curCommandGrepish;

    printf  "status=%-3d", $?   if $?;
    
    my $outputPathname  =  "/tmp/out${countString}";
    open   my $outputFile,  '>',  $outputPathname
        or   die  "could not open file '$outputPathname'";

    print  $outputFile  $output;
    close  $outputFile;

    say  "Executing prog: $countString, will compare '$outputPathname' to '$expectedOutputPathname'"  if $verbose_flag;

    #Compare to expected output (skipping first two lines)
    my $diff =
        system  "tail -n +3 $expectedOutputPathname | diff - $outputPathname > /dev/null";

    die  "output generated in '$outputPathname' does not match that expected from '$expectedOutputPathname'"
        if $diff;

    print 'OK.';
}#END assertOutputMatchesExpected()



#  Return next line from $file, skipping empty and comment lines.
sub nextLine($){
    my $file  =  $_[0]   or   die  'nextLine called without args';
    my $line = '';
    $line  =  <$file>    while  $line =~ /^(#|$)/;
    return $line;
}


#  Print $message to STDERR. If stdout is a terminal print to STDOUT as well.
sub printError($){
    my $message  =  'Command line parsing error; '  .  $_[0]  .  "\n";
    print STDERR $message;
    print STDOUT $message   unless -t *STDOUT;
    pod2usage();
}


#  Convert $_[0] to absolute path.
sub convertToAbsPath{    $_[0]  =  Cwd::abs_path(  File::Spec->rel2abs( $_[0] )  );    }


=pod


=head1 NAME

testFastaplOneLiners.pl


=head1 SYNOPSIS

B<fastaqplRegressionTester.pl>  [B<-v>|B<-n>] [numbers-to-run] I<programToTestPathname>  I<one-liners-file>

B<fastaqplRegressionTester.pl>  B<--usage>|B<--help>|B<--man>


=head1 DESCRIPTION

Utility to test B<fastapl> or B<fastqpl> one-liners.

The I<one-liners-file> stipulates fastpl one-liners and
their input files. For each one-liner a description line containing
the script number and a description of what it does, followed by a
line containing the command line to pass to B<fastapl>, a tab, and the
input filename is given.  This is best shown by example:

  01  Truncate sequences to length 39.
  -p  --main '$seq = substr( $seq, 0, 39 )'	fastaFiles/proteins.fasta

This program runs each one-liner script one by one, running fastaqpl directly each time.
When the B<-g|--generate-programs> flag is given, fastaqpl -n is used to also
generate standalone programs implementing the one lines and testing those as well.

The output is compared with the expected output stored in the
expectedOutput directory. The expected output of the first one-liner
should be found in F<expectedOutput/out01>, the second in F<expectedOutput/out02>,
etc. Each of these expected output files should start with two special lines
matching the two lines for that one-liner script found in the I<one-liners-file>.
This program will die with error if the match is not exact.
The contents of the expected output file of our running example might look like:

  01  Truncate sequences to length 39.
  -p  --main '$seq = substr( $seq, 0, 39 )'	fastaFiles/proteins.fasta
  >P52923|AIF1_YEAST
  MTINTKNIVVVGAGVFGVSVANHLYRELGGTYAIKLVTA
  ...

After skipping the first two lines, the remaining lines of the
expected output file are compared to the newly generated output and
this program dies with error if the content does not match exactly.


=head1 ARGUMENTS

=over 8

=item [I<scriptNumberBeg>[..[I<scriptNumberEnd>]]]

Only run given range of one-liners.  '7..9' means run one-liners numbered 7 through 9.  '7' means just run the seventh one-liner.  '7..' means from the seventh one-liner onward.

=item [I<programToTestPathname>]

Pathname to fastapl or fastqpl program to test.

=item I<one-liners-file>

File containing information on one-liners to test.

=back


=head1 OPTIONS

=over 8

=item B<-a|--abortOnErrorP>

Abort immediately when an error is encountered.

=item B<-A|--absolute-pathnames>

Use absolute pathnames for fastapql input files.  Useful if you want to know where the input files are.

=item B<-g|--generate-programs>

Also test programs generated with fastapql -n.

=item B<-n|--just-print>

Echo commands to stdout instead of executing them.

=item B<-v|--verbose>

Echo commands while executing them. Noop when B<-n> option
is already given.

=back


=head1 PREREQUISITES

This program requires the following modules:

  Getopt::Long    for argv parsing.
  Pod::Usage      for --usage|--help|--man message processing.
  File::Basename  to extract directory pathname of one-liners file.

=head1 AUTHOR

Paul Horton <paulh@iscb.org>

=head1 COPYRIGHT

Copyright (C) 2010

=cut
