#!/usr/bin/perl
# maker_ref_blastp_bat2.plx
use strict;
use warnings;

#required program: ncbi_blast+
#required script: blast_out_pick3.plx
my $scipt_dir    = "./scripts";
my $maker_filter_mdl_prot_suffix2= ".maker_filter.prot.clean.faa";

my $cpu_num      = shift;
my $ref_seq_dir  = shift;
my $ref_seq_file = shift;#S288C prot seq
my $maker_seq_dir= shift;
my $infile       = shift;
my $blastout_dir = shift;#default: "tmp_blastout"
if (not $blastout_dir){
	die "usage: perl maker_ref_blastp_bat2.plx <cpu num> <ref_dir> <ref.fa> <maker.faa dir> <infile> <blastout_dir>\n";
}#
my @items;
my $lbuf;

if (not -e $blastout_dir){
	system ("mkdir $blastout_dir");
}#if

system ("makeblastdb -in $ref_seq_dir/$ref_seq_file -dbtype prot");

open IN, $infile or die "fail to open $infile:$!\n";
while (<IN>){
	chomp;
	next if not $_;
	next if /^#/;
	$lbuf = $_;
	@items = split /\t/,$_;
	system ("blastp -db $ref_seq_dir/$ref_seq_file -query $maker_seq_dir/$items[0]$maker_filter_mdl_prot_suffix2  -outfmt '6 qaccver saccver pident qcovs qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 10 -num_threads $cpu_num -out $blastout_dir/$ref_seq_file.$items[0]$maker_filter_mdl_prot_suffix2.blast");
	#filter and convert
	system ("perl $scipt_dir/blast_out_pick3.plx $blastout_dir/$ref_seq_file.$items[0]$maker_filter_mdl_prot_suffix2.blast");
}#while
close IN;
