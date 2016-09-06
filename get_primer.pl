#! /usr/bin/perl
use warnings;
use strict;

##########################################################################
#
# Takes a gff3 file and finds the upstream n bases from a named feature
#  
# Outputs fasta with gff field 8 as the name (modified to remove ID=)
#
# USAGE: getprimer.pl genome gff n_bases feature
#
##########################################################################


my $my_FASTA_File = uc(join("",get_file_data(0)));
$my_FASTA_File = trim($my_FASTA_File);
my @scaffold = split(/>[a-zA-Z]*_\d*/,$my_FASTA_File);
my $nbases = $ARGV[2];


my $feature = $ARGV[3];

my @my_searches = get_file_data(1);
 
foreach (@my_searches) {
	if($_=~/$feature/) {
		my @protein = split(/\t|  */);
		$protein[0]=~s/[a-zA-Z].*_//;
		my $primer_sequence="";
		$protein[8]=~s/ID\=//;
		chop($protein[8]);
		if ($protein[6]eq"-") {
			$primer_sequence=reverse substr($scaffold[$protein[0]],$protein[4],$ARGV[2]);
			$primer_sequence=~tr/ACGT/TGCA/;
		} else {
			my $offset=$protein[3]-($ARGV[2]+1);
			my $nbases = $ARGV[2];
			if($offset<0) {
				#print STDERR "$offset\n";
				$nbases = $nbases + $offset + 1;
				$offset=0;
				
			}
			$primer_sequence=substr($scaffold[$protein[0]],$offset,$nbases);
		}
		print ">$protein[8]\n$primer_sequence\n";
	}
}
	

sub get_file_data {
#opens the input file specified from the command line
	my ($inFile) = $ARGV[$_[0]];
    	unless(open(INFILE, $inFile) ) {
		#Instructions("Cannot open input file \"$inFile\"\n\n");
        	exit;
    	}
	my @file = <INFILE>;
	chomp(@file);
	
	close INFILE;
	return @file;
}

sub trim
{

	my $string = shift;
	chomp($string);
	chop($string);
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	$string =~ s/>SCAFFOLD_//g;
	return $string;
}

