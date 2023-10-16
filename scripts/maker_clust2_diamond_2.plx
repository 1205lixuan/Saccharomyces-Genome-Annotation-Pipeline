#!/usr/bin/perl
# maker_clust2_diamond_2.plx
use strict;
use warnings;

sub ReadFasta($\@\@);#$filename,@seqname,@seq

my $script_dir = "./scripts";
#column defination
my $q_id = 0;
my $s_id = 1;
my $iden_id = 2;
my $q_cov_id = 3;
my $q_len_id = 4;
my $s_len_id = 5;
my $len_id = 6;
my $mism_id = 7;
my $gpop_id = 8;
my $q_start_id = 9;
my $q_end_id = 10;
my $s_start_id = 11;
my $s_end_id = 12;
my $evalue_id = 13;
my $score_id = 14;
my $s_cov_id = 15;	
my $id_max   = 15;	

$| = 1;
my $maker_filter_mdl_prot_suffix2= ".maker_filter.prot.clean.faa";
my $cpu_num      = shift;
my $cov_cutoff   = shift;
my $iden_cutoff  = shift;
my $ref_seq_dir  = shift;
my $ref_seq_file = shift;
my $maker_seq_dir= shift;
my $infile       = shift;
my $blastout_dir = "tmp_blastout";
my $job_prefix   = shift;

if (not $job_prefix){
	die "usage: perl maker_clust2_diamond_2.plx <cpu_num> <cov cutoff> <iden cutoff> <ref.faa_dir> <ref.fa> <maker.faa_dir> <infile_list> <job prefix>\n";
}#

my $outfile1 = $job_prefix . ".detail.out";
my $outfile2 = $job_prefix . ".maker_only_cluster.out";
my $unmap_seq = $job_prefix . ".all_unmapped_seq.fa";

my @maker_filename;

my %f_q2s;	#{$filename}->{$query} => $subject
my %seqname2file; #{$seqname} => $file
my @seqname;
my @seq;
my @items;
my @items2;
my $lbuf;
my $lno;

	
open OUT, ">$outfile1" or die "fail to open $outfile1:$!\n";
open UNMAP, ">$unmap_seq" or die "fail to open $unmap_seq:$!\n";
#read file list
open IN, $infile or die "fail to open $infile:$!\n";
while (<IN>){
	chomp;
	next if not $_;
	next if /^#/;
	$lbuf = $_;
	@items = split /\t/,$_;
	push @maker_filename, "$items[0]$maker_filter_mdl_prot_suffix2";
	#pick the top match and check the identity and coverage
#	if (-e "$blastout_dir/$ref_seq_file.$items[0]$maker_filter_mdl_prot_suffix2.blast.pick3.tab" ) {
		open BLAST, "$blastout_dir/$ref_seq_file.$items[0]$maker_filter_mdl_prot_suffix2.blast.pick3.tab" or die "fail to open file $blastout_dir/$ref_seq_file.$items[0]$maker_filter_mdl_prot_suffix2.blast.pick3.tab:$!\n";
#	}else{
#		open BLAST, "$blastout_dir/$ref_seq_file.$items[0]$maker_filter_mdl_prot_suffix2.blast.pick2.tab" or die "fail to open file $blastout_dir/$ref_seq_file.$items[0]$maker_filter_mdl_prot_suffix2.blast.pick2.tab:$!\n";
#	}#else
	while (<BLAST>){
		$lno ++;
		chomp;
		next if not $_;
		next if /^#/;
		@items2 = split /\t/,$_;
		if ($#items2 != $id_max){
			die "wrong format at line $lno in $blastout_dir/$ref_seq_file.$items[0]$maker_filter_mdl_prot_suffix2.blast.pick2/3.tab\n";
		}#if
		if ($items2[$iden_id] > $iden_cutoff and $items2[$q_cov_id] > $cov_cutoff and $items2[$s_cov_id] > $cov_cutoff ) {#
			if (not exists $f_q2s{"$items[0]$maker_filter_mdl_prot_suffix2"}->{ $items2[ $q_id ] } ) {
				$f_q2s{"$items[0]$maker_filter_mdl_prot_suffix2"}->{ $items2[ $q_id ] } = $items2[ $s_id ];
				print OUT $_,"\n";
			}#if
		}#if
	}#
	close BLAST;
	#save unselected seqs in a new file for further analysis (self clust)
	@seqname = ();
	@seq     = ();
	ReadFasta("$maker_seq_dir/$items[0]$maker_filter_mdl_prot_suffix2",@seqname,@seq);
	foreach my $i (0 .. $#seqname){
		if (exists $seqname2file{ $seqname[$i] } ) {#
			die "duplicated seq name $seqname[$i] in files $seqname2file{ $seqname[$i] } and $items[0]$maker_filter_mdl_prot_suffix2.\n";
		}else{
			$seqname2file{$seqname[$i] } = "$items[0]$maker_filter_mdl_prot_suffix2";
		}#else
		if (not exists $f_q2s{"$items[0]$maker_filter_mdl_prot_suffix2"}->{ $seqname[$i] } ) {
#			if ( $mode eq "all"){
				print UNMAP ">$seqname[$i]\n$seq[$i]\n";
#			}#
		}#if
	}#
}#while
close IN;
#if ( $mode eq "all"){
	close UNMAP;
#}
print "diamond_cluster3.2.plx start\n";

system ("perl $script_dir/diamond_cluster3.2.plx $cpu_num $iden_cutoff $cov_cutoff prot $unmap_seq");
print "diamond_cluster3.2.plx done\n";
#read results and parse
my $tmp = "$unmap_seq.clust." . $iden_cutoff . "." . $cov_cutoff . ".info";
open OUT3, ">$outfile2" or die "fail to open $outfile2:$!\n";
print OUT3 "#id";
foreach my $i (0 .. $#maker_filename){
	print OUT3 "\t",$maker_filename[$i];
}#
print OUT3 "\n";
open IN, "$tmp" or die "fail to open $tmp:$!\n";
my %tmph;
my $tmpstr;
while (<IN>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/,$_;
	%tmph = ();
	foreach my $i (2 .. $#items){
		if (exists $seqname2file{ $items[$i] } ){
			push @{ $tmph{ $seqname2file{ $items[$i] } } }, $items[$i];
		}else{
			die "should never happen!\n";
		}#else
	}#
	print OUT3 $items[0];
	foreach my $i (0 .. $#maker_filename){
		print OUT3 "\t";
		if (exists $tmph{ $maker_filename[$i] } ) {#
			print OUT3 $tmph{ $maker_filename[$i] }->[0];
			foreach my $j (1 .. $#{ $tmph{ $maker_filename[$i] } } ) {
				print OUT3 "|",$tmph{$maker_filename[$i]}->[$j];
			}#
		}#if
	}#
	print OUT3 "\n";
}#
close OUT3;

close OUT;


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
		undef $currseq;
	}#
	close ReadFasta_FH;
	
	#remove the whitespaces in the seqs
	#remove the asterisk in the seqs
	foreach my $i (0..$#{$seq_r}) {
		$seq_r->[$i] =~ s/\s//g;
#		$seq_r->[$i] =~ s/\*//g;
	}#foreach
	
	#print "ReadFasta: $filename","\t",scalar @{$seqname_r}, "\t", scalar @{$seq_r},"\n";
}#ReadFasta()
