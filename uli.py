#!/usr/local/bin/perl -w
##############################################################################
#
# PROGRAM NAME: uli - it comes after June
# DESCRIPTION: 
# Create the test, validation and training datafiles for the NN analysis
#
# loci = DRB1, A, B (and C and DQB1)
#
# INPUT FILES:
# LOCUS.expert.db - testing DNA/Serology for canonical training alleles
# LOCUS.valid.db  - testing DNA/Serology for canonical validation alleles
#
#
# OUTPUT FILES:
# LOCUS.tng.pat - training pattern file (all alleles in expert.db)
# LOCUS.val.pat - validation pattern file (all alleles in valid.db)
# LOCUS.tst.pat - testing pattern file (all alleles)
#
# WRITTEN BY:   Martin Maiers
#
# REVISION HISTORY:
# REVISION DATE         REVISED BY      DESCRIPTION
# ------- ----------    --------------  -------------------------------------
# ver1    19Jun2000     Maiers          Original
# ver6    28Apr2002     Maiers          Final A/B/DRB1 analysis
# ver7    14Aug2003     Maiers          Add C and DQB1
# ver8                  Biagini         Python Ported
#
##############################################################################
use strict;	# always
my $date = `date +%d%b%Y`; chop $date;

use lib "$ENV{'HOME'}/src/SAP";
use lib "/home/mmaiers/src/SAP";
use File::Basename;
use File::Path "/home/asivasan/NN/DBFiles/";
use SAP;
SAP::loadSAP();

##############################################################################
#
# load DNA/SER data
#
my @loci = ('A', 'B','DRB1', 'C', 'DQB1', 'DPB1');
# @loci = ('C');
my $KIRmode = 0;

foreach my $loc (@loci) {
  my %Vdnatotal; my %Vdnaserdata;
  my %dnatotal; my %dnaserdata;
  my @sertotal;

  # 
  # read training data 
  # 
  loadEXPERT($loc, \%dnatotal, \%dnaserdata, \@sertotal);

  # 
  # read validation data 
  # note: it makes no sense for validation data to contain serology 
  #       outputs that are not in the training data
  # 
  loadVALIDATION($loc, \%Vdnatotal, \%Vdnaserdata, \%dnatotal,\@sertotal);

  #
  # buld SAP tables from polymorphisms in training data
  # because if the polymorphism isn't in the training data
  # then there is no need to represent it at all.
  #

  my %allsap;
  my %allelesap = %{$SAP::locallelesap{$loc}};
  foreach my $allele (keys %{$SAP::locallelesap{$loc}}) {
    foreach my $sap (keys %{$SAP::locallelesap{$loc}{$allele}}) {
      $allsap{$sap}++ if defined $dnatotal{$allele};
    }
  }

  
  foreach my $pos (keys %{$SAP::locpossap{$loc}}) {
    #
    # see if its polymorphic
    #
    my %P=();
    foreach my $allele (keys %dnatotal) {
      my $sap = $SAP::locallelepos{$loc}{$allele}{$pos};
      if (!defined $sap) { 
        print STDERR "undefined sap for $allele at pos $pos\n";
        next;
      }
      $P{$sap}++;
    }
    if (scalar (keys %P) <=1) { 
      # remove all saps at this position
      foreach my $allele (keys %allelesap) {
        my $sap = $SAP::locallelepos{$loc}{$allele}{$pos};
        delete $allelesap{$allele}{$sap};
      }
      foreach my $sap (keys %{$SAP::locpossap{$loc}{$pos}}) {
        delete $allsap{$sap} if defined $allsap{$sap};
      }
    }
  }


  ##############################################################################
  # make output files
  ############################################################################## 

  # placeholder SAPs 
  if ($loc eq "DRB1") {
    $allsap{ZZ1}++;
    $allsap{ZZ2}++;
  }
  my @training = ();
  my @validation = ();
  my $num_training_patterns = 0;
  my $num_validation_patterns = 0;

  push @training, "# ", join(' ', sort byallsap keys %allsap), "\n";
  
  foreach my $allele (sort keys %dnatotal) {
    $num_training_patterns++;
    push @training, "# $allele\n";
    push @training, "# input\n";
    foreach my $sd (sort byallsap keys %allsap) {
      my $in = 0;
      if (defined $allelesap{$allele}{$sd}) {
        $in = 1 if $allelesap{$allele}{$sd}>0;
      }
      push @training, "$in ";
    }
    push @training, "\n";
    push @training, "# output ", join ' ', (@sertotal), "\n";
    foreach my $sd (@sertotal) {
      push @training,sprintf ("%2.2f ", defined $dnaserdata{$allele}{$sd} ? $dnaserdata{$allele}{$sd}:0);
    }
    push @training, "\n";
  }
  foreach my $allele (sort keys %Vdnatotal) {
    $num_validation_patterns++;
    push @validation, "# $allele\n";
    push @validation, "# input\n";
    foreach my $sd (sort byallsap keys %allsap) {
      my $in = 0;
      if (defined $allelesap{$allele}{$sd}) {
        $in = 1 if $allelesap{$allele}{$sd}>0;
      }
      push @validation, "$in ";
    }
    push @validation, "\n";
    push @validation, "# output ", join ' ', @sertotal, "\n";
    foreach my $sd (@sertotal) {
      push @validation,sprintf ("%2.2f ", defined $Vdnaserdata{$allele}{$sd} ? $Vdnaserdata{$allele}{$sd}:0);
    }
    push @validation, "\n";
  }

  #
  # output training pattern file
  #
  my $tngfile = "$loc.tng.pat";
  open(TNGFILE, ">$tngfile") or die "can't open $tngfile for writing";
  select TNGFILE;

 my $num_input    = scalar(keys %allsap);
 my $num_output   = scalar(@sertotal);
  print "SNNS pattern definition file V3.2\n";
  print "generated at $date\n\n\n";
  print "No. of patterns : $num_training_patterns\n";
  print "No. of input units : $num_input\n";
  print "No. of output units : $num_output\n";
 
  foreach (@training) { print; }

  #
  # output validation pattern file
  #
  my $valfile = "$loc.val.pat";
  open(VALFILE, ">$valfile") or die "can't open $valfile for writing";
  select VALFILE;
  print "SNNS pattern definition file V3.2\n";
  print "generated at $date\n\n\n";
  print "No. of patterns : $num_validation_patterns\n";
  print "No. of input units : $num_input\n";
  print "No. of output units : $num_output\n";
 
  foreach (@validation) { print; }


  #
  # output testing pattern file
  #
  my $outfile = "$loc.tst.pat";
  open(OUTFILE, ">$outfile") or die "can't open $outfile for writing";
  select OUTFILE;

  my $num_patterns = scalar (keys %allelesap);
  $num_input    = scalar(keys %allsap);
  $num_output   = scalar(@sertotal);
  print "SNNS pattern definition file V3.2\n";
  print "generated at $date\n\n\n";
  print "No. of patterns : $num_patterns\n";
  print "No. of input units : $num_input\n";
  print "No. of output units : $num_output\n";
  
  foreach my $allele (sort keys %allelesap) { 
    print "# testing $allele\n";
    print "# input\n";
    foreach my $sd (sort byallsap keys %allsap) {
      my $in = 0;
      $in = 1 if defined $allelesap{$allele}{$sd} && $allelesap{$allele}{$sd}>0;
      print "$in ";
    }
    print "\n";
    print "# output ", join ' ', (@sertotal), "\n";
    foreach my $sd (@sertotal) {
      printf ("%2.2f ", 0);
    }
    print "\n";
  }
}

exit 0;

##############################################################################
# Function: byallsap
# Description: sort routine to order by sap (position first)
##############################################################################
sub byallsap {
  $a=~/(\D)(\d+)/;
  my $as = $1;
  my $ap = $2;
  $b=~/(\D)(\d+)/;
  my $bs = $1;
  my $bp = $2;
  return $ap <=> $bp || ($as cmp $bs);
}

##############################################################################
# Function: loadEXPERT
# Description: load expert training data
##############################################################################
sub loadEXPERT {
  my($loc, $rdnatotal, $rdnaserdata, $rsertotal) = @_;

  my %sertotal;
  my $exfile = "$loc.expert.db";
  open(EX, $exfile) or die "can't open $exfile: $!";
  while(<EX>) {
    chomp;
    next if /^#/;
    my ($loc, $allele, $ser) = split/	/;
    
    # comment out for KIR ligand groups 2013-01-30

    if (!$KIRmode) {
      $ser=~tr/[A-Z\*]//d if $ser;
      next unless $ser=~/\d/;
      $ser = int $ser; 
    }

    #next unless $allele=~/\d\d\d\d/;
    $$rdnatotal{$allele} = 1;
    $$rdnaserdata{$allele}{$ser} = 1;
    $sertotal{$ser}++;
  }
  if ($KIRmode) {
    @{$rsertotal} = (sort keys %sertotal);
  } else {
    @{$rsertotal} = (sort {$a<=>$b} keys %sertotal);
  }
  close EX;
}

##############################################################################
# Function: loadVALIDATION
# Description: load validation data
##############################################################################
sub loadVALIDATION {
  my($loc, $rdnatotal, $rdnaserdata, $rtdna, $rtser) = @_;

  my %tser;
  foreach (@{$rtser}) {
    $tser{$_}++;
  }

  my $exfile = "$loc.validation.db";
  open(EX, $exfile) or die "can't open $exfile: $!";
  while(<EX>) {
    chomp;
    next if /^#/;
    my ($loc, $allele, $ser) = split/	/;
    if(!$KIRmode) {
      $ser=~tr/[A-Z\*]//d if $ser;
      next unless $ser=~/\d/;
      $ser = int $ser; 
    }
    #next unless $allele=~/\d\d\d\d/;
    if($$rtdna{$allele}) {
      warn "$loc $allele in training; can't use in validation dataset\n";
      next;
    }
    if(!defined $tser{$ser}) {
      warn "$loc $allele  ser=$ser not in training; can't use in validation\n";
      next;
    }
    $$rdnatotal{$allele} = 1;
    $$rdnaserdata{$allele}{$ser} = 1;
  }
  close EX;
}
__END__
