#!/usr/bin/perl
#blastn_annt_bat.plx
use strict;
use warnings;
my $script_dir = "./scripts";

my $cpu          = shift;
my $ref_dir      = shift;
my $inputseq_dir = shift;
my $infile       = shift;
if (not $infile){
	die "usage: perl blastn_annt_bat.plx <cpu_num> <ref_dir> <input_seq_dir> <query_fa_list>\n";
}#if
my $outfile = $infile . ".log";
my @items;
my $lbuf;
my $cnt = 0;

open JOB,"$infile" or die "fail to open $infile:$!\n";
#infile format: target_fasta_name	prefix
while(<JOB>){
	chomp;
	next if not $_;
	next if /^#/;
	$cnt ++;
	@items = split /\t/,$_;
	print "\n\nJOB ID $cnt: $items[0]\t$items[1]\n";
	
	system("perl $script_dir/blastn_annt.plx $cpu $ref_dir $inputseq_dir $items[0] $items[1]");

#	system("pigz $inputseq_dir/$items[0].genome.fna");
	#zip blastn, map files
#	system("pigz $inputseq_dir/$items[0].genome.fna.blastn");#zip to save storage, or delete it
#	system("pigz $inputseq_dir/$items[0].genome.fna.blastn.map");#zip to save storage, or delete it
}#
close JOB;

