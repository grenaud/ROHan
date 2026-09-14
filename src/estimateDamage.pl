#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use List::Util qw(min max);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;

my $mock =0;
my $noclean =0;


sub fileExists{
  my ($exeFile) = @_;

  if( -f $exeFile){
    return 1;
  }else{
    return 0;
  }
}

#fileExists() only reports, it does not act. Use this one for the binaries this script
#drives: calling fileExists() and ignoring what it returns prints "success" for a program
#that is not there, and the failure then surfaces much later as an unexplained exit code
sub checkExecutable{
  my ($exeFile,$name,$howtobuild) = @_;

  if( ! -f $exeFile ){
    print STDERR "\nERROR: cannot find ".$name.", it should be at:\n  ".$exeFile."\n".
                 "Run \"".$howtobuild."\" in the ROHan directory to build it.\n";
    exit(1);
  }

  if( ! -x $exeFile ){
    print STDERR "\nERROR: ".$name." is present but not executable:\n  ".$exeFile."\n";
    exit(1);
  }

  return 1;
}


my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);

my $minq=30;
my $minl=35;
my $subsam=1000000;
my $length=10;


my $rohan    = $pathdir."/../bin/rohan";
my $bam2prof = $pathdir."/../bam2prof/bam2prof";

my $all=1;
my $double=0;
my $single=0;
my $map="0";

my $threads=1;


print STDERR "detecting ROHan...";
checkExecutable($rohan,"the rohan binary","make");
print STDERR  "..success\n";
$rohan = abs_path($rohan);

print STDERR  "detecting bam2prof...";
checkExecutable($bam2prof,"bam2prof","make");
print STDERR  "..success\n";
$bam2prof = abs_path($bam2prof);



sub findcommandAlias{#does not work with background processes
  my ($cmdtorun) = @_;
  print STDERR "trying to find executable ". $cmdtorun."..";
  my $command;

  if($mock != 1){
    my $cmdtodetect = "which ".$cmdtorun."";
    my $out = `$cmdtodetect`;

    chomp($out);

    if( $out =~ /^(.+)$/){#found
      #my @ar = split(" ",$out);
      #$command=$ar[$#ar]."\n";
      $command=$1;
      chomp( $command );
      print STDERR ".. found: ".$command."\n";
    }else{
      $cmdtodetect = "bash -i -c \"alias |grep  ".$cmdtorun."=\"";

      $out = `$cmdtodetect`;

      if( $out =~ /^alias (\S+)=(.+)$/){
	$command=$2;
	chomp( $command );
	my @arrcmd = split //, $command;

	#removing beginning chars
	while($arrcmd[0] ne '/'){
	  shift(@arrcmd);
	}
	#removing trailing chars
	while( (ord($arrcmd[ $#arrcmd ]) == 32) || (ord($arrcmd[ $#arrcmd ]) == 39) ){
	  pop(@arrcmd);
	}

	my $command = join("",@arrcmd);


	print STDERR ".. found: ".$command."\n";
      }else{
	print STDERR "cannot find executable ".$cmdtorun."\n";
	die;
      }
    }

  }else{
    print STDERR ".. running as mock: ".$cmdtorun."\n";
    $command = $cmdtorun;
  }
  return $command;

}

sub findcommand{
  my ($cmdtorun) = @_;
  print STDERR "trying to find executable ". $cmdtorun."..";
  my $command;

  if($mock != 1){
    my $cmdtodetect = "which ".$cmdtorun."";
    my $out = `$cmdtodetect`;

    chomp($out);

    if( $out =~ /^(.+)$/){#found
      #my @ar = split(" ",$out);
      #$command=$ar[$#ar]."\n";
      $command=$1;
      chomp( $command );
      print STDERR ".. found: ".$command."\n";
    }else{
      print STDERR "cannot find executable ".$cmdtorun."\n";
      die;
    }

  }else{
    print STDERR ".. running as mock: ".$cmdtorun."\n";
    $command = $cmdtorun;
  }
  return $command;

}

sub runcmd{
  my ($cmdtorun) = @_;
  chomp($cmdtorun);
  print STDERR "running cmd ->". $cmdtorun."<-\n";
  if($mock != 1){
    #my @argstorun = ( "bash", "-i","-c", $cmdtorun );
    my @argstorun = ( $cmdtorun );
    #print STDERR "calling\n";
    if( system(@argstorun) != 0 ){
      my $why;
      if( $? == -1 ){
	$why = "it could not be run at all: ".$!;
      }elsif( $? & 127 ){
	$why = "it was killed by signal ".($? & 127);
      }else{
	my $exitcode = $? >> 8;
	$why = "it exited with code ".$exitcode;
	#127 is what the shell returns both when the program is not there and when the
	#dynamic loader cannot find a library it needs, which are the two usual causes
	if($exitcode == 127){
	  $why .= ", which means the program was not found, or a shared library it needs was not found";
	}
      }
      die "\nERROR: the following command failed, ".$why.":\n".$cmdtorun."\n";
    }
    #}else{
    print STDERR "done\n";
    #}
  }

}

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\n
 This script is a wrapper to estimate the damage patterns with the
 potentially polymorphic positions masked. This script probably will
 not work on samples with coverage less than 5X

\n\n usage:\t".$0." <options> [fasta file] [bam file] \n\n".


  " where\n".
  "  [fasta file]: The reference used for the mapping\n".
  "  [bam file]  : The BAM file for the evaluation of the damage\n".
  "\n".
" Options:\n".
  "\t-o\t\t\t\t\tOutput prefix (default: [bam file])\n".
  "\t--mock\t\t\t\t\tDo nothing, just print the commands that will be run\n".

 "\n".
  " Damage estimate\n".
  " ===================\n".

 		"\t\t--minq\t\t\t\tRequire the base to have at least this quality to be considered (Default: ".$minq.")\n".
		"\t\t--minl\t\t\t\tRequire the read to have at least this mapping quality to be considered (Default: ".$minl.")\n".#TODO
		"\t\t--subsam\t\t\tSubsample this number of reads to consider (Default: ".$subsam.")\n".#TODO

		"\t\t--length\t[length]\tDo not consider bases beyond this length  (Default: ".$length." )\n".
		"\t\t--map\t\t[bed file]\tUse these mappability tracked in bed format  (Default: ".$length." )\n".
		"\t\t--threads\t[# threads]\tNumber of threads to use  (Default: ".$threads." )\n".

		  "\n".
                " \tSubstitutions reported: specify either one of the 3 possible options:\n".

		"\t\t--all\t\t\t\tReport all potential substitutions                      (Default: used)\n".
		"\t\t--single\t\t\tUse the deamination profile of a single strand library  (Default: not on/not used)\n".
		"\t\t--double\t\t\tUse the deamination profile of a double strand library  (Default: not on/not used)\n".

		"\t\t--noclean\t\t\tDo not clean up intermediate files  (Default: not on/not used)\n".

  "\n";
  exit;
}

my $help;
my $reffile = $ARGV[$#ARGV-1];
my $bamfile = $ARGV[$#ARGV];
my $outputprefix;

usage() if ( @ARGV < 1 or
             ! GetOptions('help|?' => \$help,'o=s' => \$outputprefix,'mock' => \$mock,'noclean' => \$noclean,'minq=i' => \$minq,'threads=i' => \$threads,'minl=i' => \$minl,'subsam=i' => \$subsam,'length=i' => \$length,'all' => \$all,'single' => \$single,'double' => \$double,'map=s' => \$map) or defined $help );


if( !(defined $outputprefix) ){
  $outputprefix = substr($bamfile,0,-4)."";
}else{
}
print STDERR  "Using output prefix ... ".$outputprefix."\n";

if($double != 1){
  $all=0;
}

if($single != 1){
  $all=0;
}


#detect commands
my $cmdhead    = "head";
$cmdhead=findcommand($cmdhead);

my $cmdrm    = "rm";
$cmdrm=findcommand($cmdrm);

my $cmdgrep    = "grep";
$cmdgrep=findcommand($cmdgrep);

my $cmdawk    = "awk";
$cmdawk=findcommand($cmdawk);

my $cmdsamtools    = "samtools";
$cmdsamtools=findcommand($cmdsamtools);

my $cmdgzip    = "gzip";
$cmdgzip=findcommand($cmdgzip);

my $cmdcat    = "cat";
$cmdcat=findcommand($cmdcat);



#subsampling
my $commandsubsample = $cmdsamtools."  view -h ".$bamfile." ";
if($map ne "0"){
  $commandsubsample .= " -L ".$map." ";
}
$commandsubsample .= "| ".$cmdhead." -n ".$subsam." | ".$cmdsamtools." view -bS  /dev/stdin > ".$outputprefix."_sub".$subsam.".bam";
runcmd($commandsubsample);

my $commandsubsamplei = $cmdsamtools." index ".$outputprefix."_sub".$subsam.".bam";
runcmd($commandsubsamplei);




my $commanddepth = $cmdsamtools."   depth   ";

if($map ne "0"){
  $commanddepth = $commanddepth." -b ".$map." ";
}

$commanddepth = $commanddepth." ".$outputprefix."_sub".$subsam.".bam ";
$commanddepth = $commanddepth." | ".$cmdawk." 'BEGIN{sum=0}{sum+=\$3} END{print sum/NR}' ";

my $depth;
if($mock != 1){
  $depth = `$commanddepth`;
  chomp($depth);
  print STDERR "Found depth of coverage ".$depth."\n";
}else{
  $depth = "NA";
  print STDERR "Running command $commanddepth\n";
  print STDERR "Found depth of coverage ".$depth."\n";
}



#my $bam2profcmd1 = $cmdsamtools."  view -h ".$bamfile." ";
#if($map ne "0"){
#  $bam2profcmd1 .= " -L ".$map." ";
#}

#$bam2profcmd1 .= "| ".$cmdhead." -n ".$subsam." | ".$cmdsamtools." view -bS  /dev/stdin | ".$bam2prof." -minq ".$minq." -minl ".$minl." -length ".$length." ";
my $bam2profcmd1 = $bam2prof." -minq ".$minq." -minl ".$minl." -length ".$length." ";

if($all){
  #do nothing
}
if($single){
  $bam2profcmd1 = $bam2profcmd1." -single ";
}
if($double){
  $bam2profcmd1 = $bam2profcmd1." -double ";
}

my $deam5p1 = $outputprefix."_1.5p.prof";
my $deam3p1 = $outputprefix."_1.3p.prof";

$bam2profcmd1 = $bam2profcmd1." -5p ".$deam5p1." ";
$bam2profcmd1 = $bam2profcmd1." -3p ".$deam3p1." ";

#no "2> /dev/null" here: bam2prof reports why it gave up on stderr (a missing MD tag,
#for instance) and hiding that leaves the user with a bare exit code and nothing else
$bam2profcmd1 = $bam2profcmd1."  ".$outputprefix."_sub".$subsam.".bam";
runcmd($bam2profcmd1);


print "bam2prof done\n";
# check for nan



if($mock != 1){
  open(FILE5p1, $deam5p1) or die "cannot open ".$deam5p1."\n";
  while(my $line = <FILE5p1>){
    chomp($line);
    my @arrayDeam = split( "\t",$line);

    if($#arrayDeam > 0){
      if( $arrayDeam[ 0 ] eq "A>C"){
	next;
      }
      for(my $i=0;$i<=$#arrayDeam;$i++){
	if($arrayDeam[ $i ] eq "-nan"){
	  die "ERROR: found a not a number (NAN) on line ".$line." in file ".$deam5p1."\n"."rerun with a lower --minq (used: ".$minq.") ";
	}
      }
    }
  }
  close(FILE5p1);
}else{
  print STDERR "opening $deam5p1\n";
}

if($mock != 1){
  open(FILE3p1, $deam3p1) or die "cannot open ".$deam3p1."\n";
  while(my $line = <FILE3p1>){
    chomp($line);
    my @arrayDeam = split( "\t",$line);

    if($#arrayDeam > 0){
      if( $arrayDeam[ 0 ] eq "A>C"){
	next;
      }
      for(my $i=0;$i<=$#arrayDeam;$i++){
	if($arrayDeam[ $i ] eq "-nan"){

	  die "ERROR: found a not a number (NAN) on line ".$line." in file ".$deam3p1."\n"."rerun with a lower --minq (used: ".$minq.") ";
	}
      }
    }
  }
  close(FILE3p1);
}else{
  print STDERR "opening $deam3p1\n";
}




if( fileExists($reffile.".fai") ){
  #fine
}else{
  my $samtoolscmdfai = $cmdsamtools." faidx ".$reffile;
  runcmd($samtoolscmdfai);
}


my $catbed = $cmdcat."  ".$reffile.".fai   | ".$cmdawk." '{print \$1\"\\t\"(1)\"\\t\"\$2}' > ".$outputprefix."_genome.bed";
runcmd($catbed);



my $rohancommand1 = $rohan." --lambda ".$depth." -t ".$threads." -v --gl --first --nohmm --deam5p ".$deam5p1." --deam3p ".$deam3p1."  -o ".$outputprefix." ";

$rohancommand1 .= " --bed ".$outputprefix."_genome.bed ";

if($map ne "0"){
  $rohancommand1 .= " --map ".$map." ";
}

$rohancommand1 = $rohancommand1." ".$reffile." ".$outputprefix."_sub".$subsam.".bam";
runcmd($rohancommand1);

my $extractsegsites=$cmdgzip." -c -d ".$outputprefix.".vcf.gz | ".$cmdgrep." -v \"^#\" | ".$cmdgrep." \"0/1\\|1/1\"  | ".$cmdawk." '{print \$1\"\\t\"(\$2-1)\"\\t\"\$2}' | ".$cmdgzip." -c > ".$outputprefix.".bed.gz";
runcmd($extractsegsites);





my $bam2profcmd2 = $bam2prof." -minq ".$minq." -minl ".$minl." -length ".$length." ";

if($all){
  #do nothing
}
if($single){
  $bam2profcmd2 = $bam2profcmd2." -single ";
}
if($double){
  $bam2profcmd2 = $bam2profcmd2." -double ";
}

my $deam5p2 = $outputprefix.".5p.prof";
my $deam3p2 = $outputprefix.".3p.prof";

$bam2profcmd2 = $bam2profcmd2." -5p ".$deam5p2." ";
$bam2profcmd2 = $bam2profcmd2." -3p ".$deam3p2." ";

$bam2profcmd2 = $bam2profcmd2." -mask ".$outputprefix.".bed.gz  ".$outputprefix."_sub".$subsam.".bam";
runcmd($bam2profcmd2);





#cleaning up
if($noclean != 1){
  my $commandrmsubsample = $cmdrm." -f ".$outputprefix."_sub".$subsam.".bam";
  runcmd($commandrmsubsample);

  my $commandrmsubsamplei = $cmdrm." -f ".$outputprefix."_sub".$subsam.".bam.bai";
  runcmd($commandrmsubsamplei);

  my $commandrm5p = $cmdrm." -f ".$outputprefix."_1.5p.prof";
  runcmd($commandrm5p);

  my $commandrm3p = $cmdrm." -f ".$outputprefix."_1.3p.prof";
  runcmd($commandrm3p);

  my $commandrmest = $cmdrm." -f ".$outputprefix."hEst.gz";
  runcmd($commandrmest);

  my $commandrminfo = $cmdrm." -f ".$outputprefix."rginfo.gz";
  runcmd($commandrminfo);

  my $commandrmvcf = $cmdrm." -f ".$outputprefix."vcf.gz";
  runcmd($commandrmvcf);

  my $commandrmbedgz = $cmdrm." -f ".$outputprefix."bed.gz";
  runcmd($commandrmbedgz);

  my $commandrmbed = $cmdrm." -f ".$outputprefix."_genome.bed";
  runcmd($commandrmbed);
}

if($mock != 1){
  open(FILE5p2, $deam5p2) or die "cannot open ".$deam5p2."\n";
  while(my $line = <FILE5p2>){
    chomp($line);
    my @arrayDeam = split( "\t",$line);

    if($#arrayDeam > 0){
      if( $arrayDeam[ 0 ] eq "A>C"){
	next;
      }
      for(my $i=0;$i<=$#arrayDeam;$i++){
	if($arrayDeam[ $i ] eq "-nan"){
	  die "ERROR: found a not a number (NAN) on line ".$line." in file ".$deam5p2."\n"."rerun with a lower --minq (used: ".$minq.") ";
	}
      }
    }
  }
  close(FILE5p2);
}else{
  print STDERR "opening $deam5p2\n";
}

if($mock != 1){
  open(FILE3p2, $deam3p2) or die "cannot open ".$deam3p2."\n";
  while(my $line = <FILE3p2>){
    chomp($line);
    my @arrayDeam = split( "\t",$line);

    if($#arrayDeam > 0){
      if( $arrayDeam[ 0 ] eq "A>C"){
	next;
      }
      for(my $i=0;$i<=$#arrayDeam;$i++){
	if($arrayDeam[ $i ] eq "-nan"){

	  die "ERROR: found a not a number (NAN) on line ".$line." in file ".$deam3p2."\n"."rerun with a lower --minq (used: ".$minq.") ";
	}
      }
    }
  }
  close(FILE3p2);
}else{
  print STDERR "opening $deam3p2\n";
}


if($mock != 1){
print STDERR "  -- Program done -- \n";
print STDERR " 5' deamination file:  ".$deam5p2."\n";
print STDERR " 3' deamination file:  ".$deam3p2."\n";
}


warn $outputprefix;


