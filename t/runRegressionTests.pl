#!/usr/bin/perl -w
use Cwd;
use File::Basename;
use File::Spec::Functions;
use feature qw(say);

#  Prepend local modules directory onto Perl lib path
my $modulesPath=   catfile(  dirname(Cwd::abs_path $0)=> '..'=> 'modules/perl5');


$ENV{PERL5LIB}=  exists $ENV{PERL5LIB}?  "$ENV{PERL5LIB}:$modulesPath"  :  $modulesPath;

say  "\n$0 Running fastapl5...";
system  "t/fastaqplRegressionTester.pl -g @ARGV fastapl5 t/regression/oneLinerLists/fastapl5oneLiners.txt"
    and   exit 1;    #Exit if error encountered.


say  "\n$0 Running fastqpl5...";
system  "t/fastaqplRegressionTester.pl -g @ARGV fastqpl5 t/regression/oneLinerLists/fastqpl5oneLiners.txt";


