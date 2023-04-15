#!/usr/bin/perl
# merge_annt3.plx
use strict;
use warnings;
#also convert non-ATGC nt to N
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sub ReadFasta($\@\@);#$filename,@seqname,@seq
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sub WriteFasta($\@\@);#$filename,@seqname,@seq
sub rc($);#$inseq

my $min_aa_len = 15;#the min length for the S288C prot was 16 aa

my $blastn_mdl_gff_suffix = ".blastn.gff";
my $blastn_mdl_gene_suffix= ".blastn.gene.fna";
my $blastn_mdl_mrna_suffix= ".blastn.mrna.fna";
my $blastn_mdl_prot_suffix= ".blastn.prot.faa";
my $blastn_mdl_cds_suffix = ".blastn.cds.fna";

my $maker_mdl_gff_suffix = ".maker.gff";
my $maker_mdl_gene_suffix= ".maker.gene.fna";
my $maker_mdl_mrna_suffix= ".maker.mrna.fna";
my $maker_mdl_prot_suffix= ".maker.prot.faa";
my $maker_mdl_cds_suffix = ".maker.cds.fna";

my $maker_filter_mdl_gff_suffix = ".maker_filter.gff";
my $maker_filter_mdl_gene_suffix= ".maker_filter.gene.fna";
my $maker_filter_mdl_mrna_suffix= ".maker_filter.mrna.fna";
my $maker_filter_mdl_prot_suffix= ".maker_filter.prot.faa";
my $maker_filter_mdl_cds_suffix = ".maker_filter.cds.fna";
my $maker_filter_mdl_prot_suffix2= ".maker_filter.prot.clean.faa";

my $merge_mdl_gff_suffix = ".final.gff";
my $merge_mdl_gene_suffix= ".final.gene.fna";
my $merge_mdl_mrna_suffix= ".final.mrna.fna";
my $merge_mdl_prot_suffix= ".final.prot.faa";
my $merge_mdl_cds_suffix = ".final.cds.fna";
my $merge_mdl_prot_suffix2= ".final.prot.clean.faa";


#seq_id match!!!

my $seq_dir = shift;
my $pref  = shift;
if (not $pref){
	die "usage: perl merge_annt3.plx <seq_dir> <query_fa_prefix>\n";
}#if
my $outfile = $pref . ".log";
my @items;
my @items2;

my @name;
my @seq;
my $lbuf;
my $tmp;
my $ch;
my $tmp_prot_id;
my $tmp_gene_id;
my %protid_2_geneid;

my @protseq;
my @protname;

my %blastn_ctg_pos = ();	#{$ctg}->{$pos} ++
my %maker_pick     = ();	#{$maker_gene_id} ++;

#read the blastn mdl coordinate
open BLN, "$seq_dir/$pref$blastn_mdl_gff_suffix" or die "fail to open $seq_dir/$pref$blastn_mdl_gff_suffix:$!\n";
while (<BLN>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/;
	if ($items[2] eq "gene"){
		foreach my $i ($items[3] .. $items[4]){
			$blastn_ctg_pos{ $items[0] }->{ $i } ++;
		}
	}elsif ($items[2] eq "CDS" ) {
		$tmp_prot_id = "";
		$tmp_gene_id = "";
		if ($items[8] =~ /protein_id=([^,;]+)/ ) {
			$tmp_prot_id = $1;
		}else{
			die "fail to parse protein_id in $items[8]\n";
		}#
		if ($items[8] =~ /locus_tag=([^,;]+)/ ) {
			$tmp_gene_id = $1;
		}else{
			die "fail to parse locus_tag in $items[8]\n";
		}#else
		$protid_2_geneid{ $tmp_prot_id } = $tmp_gene_id;
	}#elsif
}
close BLN;
#read the maker mdl coordinate
open MKR, "$seq_dir/$pref$maker_mdl_gff_suffix" or die "fail to open $seq_dir/$pref$maker_mdl_gff_suffix:$!\n";


while (<MKR>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/;
	
	if ($items[2] eq "CDS"){#
		$tmp_prot_id = "";
		$tmp_gene_id = "";
		if ($items[8] =~ /protein_id=([^,;]+)/ ) {
			$tmp_prot_id = $1;
		}else{
			die "fail to parse protein_id in $items[8]\n";
		}#
		if ($items[8] =~ /locus_tag=([^,;]+)/ ) {
			$tmp_gene_id = $1;
		}else{
			die "fail to parse locus_tag in $items[8]\n";
		}#else
		$protid_2_geneid{ $tmp_prot_id } = $tmp_gene_id;
	}elsif ( $items[2] ne "gene" ) {
		if ($items[8] =~ /locus_tag=([^;]+)/){
			$tmp = 0;
			foreach my $i ($items[3] .. $items[4]){
				$tmp++ if exists $blastn_ctg_pos{ $items[0] }->{$i};
			}#
			$maker_pick{ $1 } ++ if not $tmp;
		}else{
			die "fail to parse gene ID from string $items[8].\n";
		}#
	}#
}#
close MKR;


#read prot seq
ReadFasta ("$seq_dir/$pref$maker_mdl_prot_suffix", @protname, @protseq);
my $cnt_middle_stop = 0;#excluding those containing stop * in the middle of the seq 

my $flag = 0;
my $tmplen = 0;
my $tmpstr;
my $tmpcnt = 0;
foreach my $i (0 .. $#protname){
	$flag = 0;
	#check the prot seq
	$tmp = substr $protseq[$i], 0, (length $protseq[$i]) - 1;
	if ( $tmp =~ /\*/ ) {
		$cnt_middle_stop ++;
		$flag ++;
	}#if
	$tmpstr = $protseq[$i];
	$tmpstr =~ s/\*//g;
	$tmplen = length $tmpstr;
	if ( $flag or $tmplen < $min_aa_len) {
		if (exists $maker_pick{ $protid_2_geneid{ $protname[$i] } } ) {
			delete $maker_pick{ $protid_2_geneid{ $protname[$i] } };
		}#if
	}#if
}#

#produce merged annotation files
## gff
##also get the names of mrna/cds/prot
open FIL, ">$seq_dir/$pref$maker_filter_mdl_gff_suffix" or die "fail to open $seq_dir/$pref$maker_filter_mdl_gff_suffix:$!\n";
open MKR, "$seq_dir/$pref$maker_mdl_gff_suffix" or die "fail to open $seq_dir/$pref$maker_mdl_gff_suffix:$!\n";
while(<MKR>){
	chomp;
	next if not $_;
	next if /^#/;
	$lbuf = $_;
	@items = split /\t/,$_;
	if ($items[8] =~ /locus_tag=([^,;]+)/ ) {
		if (exists $maker_pick{$1}){
			print FIL $lbuf,"\n";
			#also get other id
			if ($items[8] =~ /protein_id=([^,;]+)/ ) {
				$maker_pick{$1} ++;
			}
			if ($items[8] =~ /transcript_id=([^,;]+)/ ) {
				$maker_pick{$1} ++;
			}
			if ($items[8] =~ /^ID=([^,;-]+)/ ) {
				$maker_pick{$1} ++;
			}#if
		}#if
	}#if
}
close MKR;
close FIL;

system ("cat $seq_dir/$pref$blastn_mdl_gff_suffix $seq_dir/$pref$maker_filter_mdl_gff_suffix > $seq_dir/$pref$merge_mdl_gff_suffix");

@name=();
@seq=();
ReadFasta("$seq_dir/$pref$maker_mdl_gene_suffix", @name,@seq);
open OUT, ">$seq_dir/$pref$maker_filter_mdl_gene_suffix" or die "fail to open $seq_dir/$pref$maker_filter_mdl_gene_suffix:$!\n";
foreach my $i (0 .. $#name){
	if (exists $maker_pick{ $name[$i] } ) {
		print OUT ">",$name[$i],"\n",$seq[$i],"\n";
	}#if
}#
close OUT;
system ("cat $seq_dir/$pref$blastn_mdl_gene_suffix $seq_dir/$pref$maker_filter_mdl_gene_suffix > $seq_dir/$pref$merge_mdl_gene_suffix");

	
@name=();
@seq=();
ReadFasta("$seq_dir/$pref$maker_mdl_mrna_suffix", @name,@seq);
open OUT, ">$seq_dir/$pref$maker_filter_mdl_mrna_suffix" or die "fail to open $seq_dir/$pref$maker_filter_mdl_mrna_suffix:$!\n";
foreach my $i (0 .. $#name){
	if (exists $maker_pick{ $name[$i] } ) {
		print OUT ">",$name[$i],"\n",$seq[$i],"\n";
	}#if
}#
close OUT;
system ("cat $seq_dir/$pref$blastn_mdl_mrna_suffix $seq_dir/$pref$maker_filter_mdl_mrna_suffix > $seq_dir/$pref$merge_mdl_mrna_suffix");


	
@name=();
@seq=();
ReadFasta("$seq_dir/$pref$maker_mdl_prot_suffix", @name,@seq);
open OUT, ">$seq_dir/$pref$maker_filter_mdl_prot_suffix" or die "fail to open $seq_dir/$pref$maker_filter_mdl_prot_suffix:$!\n";
open OUT2, ">$seq_dir/$pref$maker_filter_mdl_prot_suffix2" or die "fail to open $seq_dir/$pref$maker_filter_mdl_prot_suffix2:$!\n";
foreach my $i (0 .. $#name){
	if (exists $maker_pick{ $name[$i] } ) {
		print OUT ">",$name[$i],"\n",$seq[$i],"\n";
		$seq[$i] =~ s/\*//g;
		print OUT2 ">",$name[$i],"\n",$seq[$i],"\n";
	}#if
}#
close OUT;
close OUT2;

system ("cat $seq_dir/$pref$blastn_mdl_prot_suffix $seq_dir/$pref$maker_filter_mdl_prot_suffix > $seq_dir/$pref$merge_mdl_prot_suffix");

open IN, "$seq_dir/$pref$merge_mdl_prot_suffix" or die "fail to open $seq_dir/$pref$merge_mdl_prot_suffix:$!\n";
open OUT, ">$seq_dir/$pref$merge_mdl_prot_suffix2" or die "fail to open $seq_dir/$pref$merge_mdl_prot_suffix2:$!\n";
while (<IN>){
	$lbuf = $_;
	$lbuf =~ s/\*//g;	
	print OUT $lbuf;
}
close IN;
close OUT;

	
@name=();
@seq=();
ReadFasta("$seq_dir/$pref$maker_mdl_cds_suffix", @name,@seq);
open OUT, ">$seq_dir/$pref$maker_filter_mdl_cds_suffix" or die "fail to open $seq_dir/$pref$maker_filter_mdl_cds_suffix:$!\n";
foreach my $i (0 .. $#name){
	if (exists $maker_pick{ $name[$i] } ) {
		print OUT ">",$name[$i],"\n",$seq[$i],"\n";
	}#if
}#
close OUT;
system ("cat $seq_dir/$pref$blastn_mdl_cds_suffix $seq_dir/$pref$maker_filter_mdl_cds_suffix > $seq_dir/$pref$merge_mdl_cds_suffix");
	



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
#		$seq_r->[$i] =~ s/\*//g;
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
