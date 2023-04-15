#!/usr/bin/perl
#blastn_annt_2_family.plx
use strict;
use warnings;


sub ReadFasta($\@\@);#$filename,@seqname,@seq
sub WriteFasta($\@\@);#$filename,@seqname,@seq
sub rc($);#$inseq

my $blastn_gff_suffix = ".blastn.gff";
my $ref_fa  = "GCF_000146045.2_R64_protein.shortname.faa";
my $ref_tag = "S288C";

my $ref_dir      = shift;
my $inputseq_dir = shift;
my $infile       = shift;
my $outfile      = shift;
if (not $outfile){
	die "usage: perl blastn_annt_2_family.plx <ref_dir> <input_seq_dir> <query_fa_list> <outfile_name>\n";
}#if

my @seqname;
my @seq;
my @items;
my @items2;
my $lbuf;
my $cnt = 0;

my @tag_order = ();

my %ref_2_mdl;#{$ref_prot_id}->{ $spec_tag }->{$proteins in spec $spec_tag} ++

open JOB,"$infile" or die "fail to open $infile:$!\n";
#infile format: prefix	target_fasta_name	
while(<JOB>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/,$_;
	push @tag_order, $items[0];
	print "\n\nJOB ID $cnt: $items[0]\t$items[1]\n";
	open GFF, "$inputseq_dir/$items[0]$blastn_gff_suffix" or die "fail to open $inputseq_dir/$items[0]$blastn_gff_suffix:$!\n";
	while (<GFF>){
		chomp;
		next if not $_;
		next if /^#/;
		@items2 = split /\t/;
		if ($items2[2] eq "CDS"){
			if ($items2[8] =~ /protein_id=([^,;]+).+ref_gene_mdl=[^;]+protein_id:([^;]+)/){
				$ref_2_mdl{ $2 }->{ $items[0] }->{$1} ++;
			}#
		}#if
	}#
	close GFF;
}
close JOB;

open OUT, ">$outfile" or die "fail to open $outfile:$!\n";
print OUT $ref_tag;
foreach my $i (@tag_order){
	print OUT "\t",$i;
}
print OUT "\n";

open IN, "$ref_dir/$ref_fa" or die "fail to open $ref_dir/$ref_fa:$!\n";
while(<IN>){
	chomp;
	next if not $_;
	if ( /^>(.+)$/ ) {
		print OUT $1;
		foreach my $i (0 .. $#tag_order){
			if (exists $ref_2_mdl{ $1 }->{ $tag_order[$i] } ) {
				@items = sort keys %{ $ref_2_mdl{ $1 }->{ $tag_order[$i] } };
				$lbuf = join "|", @items;
				print OUT "\t", $lbuf;
			}else{
				print OUT "\t";
			}#else
		}#
		print OUT "\n";
	}#if
}
close IN;
close OUT;

			
