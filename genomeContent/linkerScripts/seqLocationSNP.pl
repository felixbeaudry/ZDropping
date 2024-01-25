#!/usr/bin/perl

#get genome alignments for a seq
#perl seqLocation.pl genomeSAMfile outfile

use strict;
use warnings;

my ($gfile, $outfile) = @ARGV;

open OUT, ">$outfile" or die;
$,="\t";

my %data = ();

open GFILE, "$gfile" or die;

#print header
print OUT "SNPname\tNewScaff\tProbePos\tProbeOrientation\tSNPpos\tnewSNPname\n";

#read in alignments to genome
while(<GFILE>)
{
	chomp;
	#read in headers
	if ($_ =~ /^\@/)
	{
		next;
	}
	else
	{
		my @line = split "\t";
		my $snp;
		if ($line[1] == 0)
		{
			$snp = $line[3] + 60;
			print OUT "$line[0]\t$line[2]\t$line[3]\t$line[1]\t$snp\t$line[2]p$snp\n";
		}
		elsif ($line[1] == 16)
		{
			$snp = $line[3] - 60;
			print OUT "$line[0]\t$line[2]\t$line[3]\t$line[1]\t$snp\t$line[2]p$snp\n";
		}	
		else
		{
			$snp = ' ';
			print OUT "$line[0]\t$line[2]\t$line[3]\t$line[1]\t$snp\t$line[2]p$snp\n";
		}	
	}
}

close GFILE;