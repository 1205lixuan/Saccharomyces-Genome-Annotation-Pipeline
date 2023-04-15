#!/usr/bin/perl
# blast_cluster3.1.plx
use strict;
use warnings;


#  This script is used to cluster sequences based on the clean BLAST results.
#  Currently, the single linkage clustering method was implemented to do the analysis.
#  There are three steps to do the analysis.
#  STEP 1. BLAST
#  STEP 2. clean the BLAST results
#  STEP 3. cluster analysis


sub ReadFasta($\@\@);#$filename,@seqname,@seq
sub WriteFasta($\@\@);#$filename,@seqname,@seq
sub ClustSimSingleLink2(\@\%\@);#@list,%qs,@cluster


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



my $cpu_num    = shift;
my $cov_cutoff = shift;
my $iden_cutoff= shift;
my $type = shift;
my $fa = shift;
if (not $fa){
	die "usage: perl blast_cluster3.1.plx <cpu_num> <coverage:0-100> <identity:0-100> <seq_type:prot or nucl> <fasta>\n";
}#if
if ($type ne "nucl" and $type ne "prot"){
	die "wrong type: $type\nusage: perl blast_cluster3.1.plx <cpu_num> <coverage:0-100> <identity:0-100> <seq_type:prot or nucl> <fasta>\n";
}#if
	
my $outfile1 = $fa . ".clust." . $cov_cutoff . "." . $iden_cutoff . ".info";
my $outfile2 = $fa . ".clust." . $cov_cutoff . "." . $iden_cutoff . ".fa";

my @seqname;	#tmp
my @seq;		#tmp
my %name2idx;
my %qs_iden;
my %qs_qcov;
my %qs_scov;
my @clust;

my @items;
my $lno = 0;

ReadFasta($fa, @seqname, @seq);
my $seqnum = @seqname;
foreach my $i (0 .. $#seqname){
	$name2idx{ $seqname[$i] } = $i;
}#
#if (not -e "$fa.selfblast.pick3.tab") {
#	if (not -e "$fa.selfblast" ) {
		if ($type eq "nucl"){
			system ("makeblastdb -in $fa -dbtype nucl");
			system ("blastn -db $fa -query $fa -num_alignments $seqnum -outfmt '6 qaccver saccver pident qcovs qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 10 -num_threads $cpu_num -out $fa.selfblast");
		}elsif($type eq "prot") {
			system ("makeblastdb -in $fa -dbtype prot");
			system ("blastp -db $fa -query $fa -num_alignments $seqnum -outfmt '6 qaccver saccver pident qcovs qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 10 -num_threads $cpu_num -out $fa.selfblast");
		}#
				
#	}#if
	system ("perl $script_dir/blast_out_pick3.plx $fa.selfblast");
#}#

open IN, "$fa.selfblast.pick3.tab" or die "fail to open $fa.selfblast.pick3.tab:$!\n";
while (<IN>){
	$lno ++;
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/,$_;
	if ($#items != $id_max){
		die "wrong format at line $lno in $fa.selfblast.pick3.tab\n";
	}#if
	if ($items[$iden_id] > $iden_cutoff and $items[$q_cov_id] > $cov_cutoff and $items[$s_cov_id] > $cov_cutoff ) {#
		$qs_iden{ $items[ $q_id ] }->{ $items[ $s_id ] } = $items[ $iden_id ];
		$qs_qcov{ $items[ $q_id ] }->{ $items[ $s_id ] } = $items[ $q_cov_id ];
		$qs_scov{ $items[ $q_id ] }->{ $items[ $s_id ] } = $items[ $s_cov_id ];
	}#if
}#
close IN;


#cluster
ClustSimSingleLink2( @seqname, %qs_iden, @clust );


#output clusters and seq of representative of each cluster
open OUT, ">$outfile1" or die "fail to open $outfile1:$!\n";
print OUT "#cov_cutoff:\t",$cov_cutoff,"\n";
print OUT "#iden_cutoff:\t",$iden_cutoff,"\n";
print OUT "#seq_num:\t",$seqnum,"\n";
print OUT "#clust_num:\t",scalar @clust,"\n";
print "cov_cutoff:\t",$cov_cutoff,"\n";
print "iden_cutoff:\t",$iden_cutoff,"\n";
print "seq_num:\t",$seqnum,"\n";
print "clust_num:\t",scalar @clust,"\n";
print OUT "#id\tsize\tlist\n";
foreach my $i ( 0 .. $#clust ) {#
	print OUT $i,"\t",scalar @{$clust[$i]};
	foreach my $j (0 .. $#{ $clust[$i] }){
		print OUT "\t",$clust[$i]->[$j];
	}
	print OUT "\n";
}#foreach my $i

close OUT;

#output the seq
open OUT, ">$outfile2" or die "fail to open $outfile2:$!\n";
foreach my $i (0 .. $#clust){
	print OUT ">", $seqname[ $name2idx{ $clust[$i]->[0] } ],"\n",$seq[ $name2idx{ $clust[$i]->[0] } ],"\n";
}#
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


sub WriteFasta($\@\@){#this version does not check the number of seqname/seq
	my $filename   = shift;
	my $seqname_r  = shift;
	my $seq_r      = shift;
	
	#print "WriteFasta: $filename","\t",scalar @{$seqname_r}, "\t", scalar @{$seq_r},"\n";

	open WriteFasta_FH, ">$filename" or die "sub WriteFasta(): fail to open file $filename:$!\n";
	foreach my $i (0..$#{$seqname_r}) {
		print WriteFasta_FH ">",$seqname_r->[$i],"\n";
		print WriteFasta_FH $seq_r->[$i],"\n";
	}#foreach
	close WriteFasta_FH;
}#sub WriteFasta()
	
	
sub ClustSimSingleLink2(\@\%\@){#@list,%qs,@cluster
	my $inlist_r = shift;
	my $inhash_r = shift;
	my $clust_r = shift;
	
	my $clust_num = 0;
	my %tmpclust;

	my $id1;
	my $id2;
	my %item2clust;
	#init clust, each item represents a different cluster
	foreach my $i (0 .. $#{ $inlist_r } ) {
		$item2clust{ $inlist_r->[$i] } = $i;
		push @{ $tmpclust{ $i } }, $inlist_r->[$i];
	}#
	
	#read each relation
	foreach my $i (keys %{$inhash_r}){
		foreach my $j (keys %{$inhash_r->{$i}}){
			$id1 = $item2clust{ $i };
			$id2 = $item2clust{ $j };
			if ($id1 == $id2){#already in the same cluster, do nothing
				#
			}elsif ($id1 < $id2){#merge cluster id2 to cluster id1
				foreach my $k (0 .. $#{$tmpclust{$id2}} ) {
					$item2clust{ $tmpclust{$id2}->[$k] } = $id1;
				}#
				push @{ $tmpclust{ $id1 } }, @{ $tmpclust{$id2} };
				delete $tmpclust{$id2};
			}else{#id1 > $id2
				foreach my $k (0 .. $#{$tmpclust{$id1}}){
					$item2clust{ $tmpclust{$id1}->[$k] } = $id2;
				}
				push @{ $tmpclust{ $id2 } }, @{ $tmpclust{$id1} };
				delete $tmpclust{$id1};
			}#else
		}#
	}#

	#copy the results
	foreach my $i (sort keys %tmpclust){
		push @{ $clust_r->[ $clust_num++ ] }, @{$tmpclust{$i}};
	}
	
}#sub ClustSimSingleLink2(\@\%\@) #@list,%qs,@cluster
	
		
	
