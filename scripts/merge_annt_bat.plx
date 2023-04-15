#!/usr/bin/perl
# merge_annt_bat.plx
use strict;
use warnings;

#
my $script_dir = "./scripts";
my $seq_dir = shift;
my $infile  = shift;
if (not $infile){
	die "usage: perl merge_annt_bat.plx <seq_dir> <query_fa_list>\n";
}#if
#my $outfile = $infile . ".log";
my @items;


open IN, "$infile" or die "fail to open $infile:$!\n";
while (<IN>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/;
	
	system("perl $script_dir/merge_annt3.plx $seq_dir $items[0]");
}#
close IN;






sub ReadFasta($\@\@){
	my $filename  = shift;
	my $seqname_r = shift;
	my $seq_r     = shift;
	
	my $currseq = "";
	
	my $lineno   = 0;
	
	open ReadFasta_FH, "$filename" or die "sub ReadFasta(): fail to open file $filename: $!\n";
	while (<ReadFasta_FH>) {
		$lineno ++;
		chomp;
		if(! length $_){
			next;
		}#blank line
		if(m/^>(.+)$/){#this is the name of a new seq
			if ($currseq) {
				push @{$seq_r}, $currseq;
				$currseq = "";
			}#
			push @{$seqname_r},$1;
		}else{#part of the current seq
			$currseq = $currseq . $_;
		}#else
	}#while
	if ($currseq) {
		push @{$seq_r},$currseq;
		$currseq = "";
	}#
	close ReadFasta_FH;
	
	#remove the whitespaces in the seqs
	#remove the asterisk in the seqs
	foreach my $i (0..$#{$seq_r}) {
		$seq_r->[$i] =~ s/\s//g;
		$seq_r->[$i] =~ s/\*//g;
	}#foreach
	
	print "ReadFasta: $filename","\t",scalar @{$seqname_r}, "\t", scalar @{$seq_r},"\n";
}#ReadFasta()


sub WriteFasta($\@\@){#this version does not check the number of seqname/seq
	my $filename   = shift;
	my $seqname_r  = shift;
	my $seq_r      = shift;
	
	print "WriteFasta: $filename","\t",scalar @{$seqname_r}, "\t", scalar @{$seq_r},"\n";

	open WriteFasta_FH, ">$filename" or die "sub WriteFasta(): fail to open file $filename:$!\n";
	foreach my $i (0..$#{$seqname_r}) {
		print WriteFasta_FH ">",$seqname_r->[$i],"\n";
		print WriteFasta_FH $seq_r->[$i],"\n";
	}#foreach
	close WriteFasta_FH;
}#sub WriteFasta()
	
sub rc($){
	my $in =shift;
	my @ch = split //,$in;
	@ch = reverse @ch;
	foreach my $i (0 .. $#ch){
		if ($ch[$i] ne "A" 
		and $ch[$i] ne "T"
		and $ch[$i] ne "G"
		and $ch[$i] ne "C"
		and $ch[$i] ne "a"
		and $ch[$i] ne "t"
		and $ch[$i] ne "g"
		and $ch[$i] ne "c"
		and $ch[$i] ne "-"
		and $ch[$i] ne " " ) {
			$ch[$i] = "N";
		}#if
	}#foreach my $i
	my $out = join "",@ch;
	$out =~ tr/ATGCatgc/TACGtacg/;
	$out;
}#sub rc($)
