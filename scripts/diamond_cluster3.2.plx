#!/usr/bin/perl
# diamond_cluster3.2.plx
use strict;
use warnings;

#  This script is used to cluster sequences based on the clean BLAST results.
#  Currently, the single linkage clustering method was implemented to do the analysis.
#  There are three steps to do the analysis.
#  STEP 1. DIAMOND BLAST
#  STEP 2. clean the DIAMOND BLAST results
#  STEP 3. cluster analysis


sub ReadFasta($\@\@);#$filename,@seqname,@seq
sub WriteFasta($\@\@);#$filename,@seqname,@seq

$|=1;
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
	die "usage: perl diamond_cluster3.2.plx <cpu_num> <coverage:0-100> <identity:0-100> <seq_type:prot or nucl> <fasta>\n";
}#if
if ($type ne "nucl" and $type ne "prot"){
	die "wrong type: $type\nusage: perl diamond_cluster3.2.plx <cpu_num> <coverage:0-100> <identity:0-100> <seq_type:prot or nucl> <fasta>\n";
}#if
	
my $outfile1 = $fa . ".clust." . $cov_cutoff . "." . $iden_cutoff . ".info";
my $outfile2 = $fa . ".clust." . $cov_cutoff . "." . $iden_cutoff . ".fa";

my @seqname;	#tmp
my @seq;	#tmp
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
system ("diamond makedb --in $fa -d db.$fa");
system ("diamond blastp -d db.$fa -q $fa -o $fa.selfblast -p $cpu_num --very-sensitive -e 10 -k $seqnum --max-hsps 0 -f 6 qseqid sseqid pident qcovhsp qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore");
system ("perl $script_dir/blast_out_pick3.plx $fa.selfblast");


my %name2clst;#{$seqname} => idx of the corresponding cluster
my @clst;#[idx]->[idx] => name
my $clstmaxidx = -1;
my @coc;#[idx]->[idx] => clust idx; cluster_of_cluster
my $cocmaxidx  = -1;
my @clst2coc;#[idx,i.e. idx of @clst] => idx of @coc


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
	next if $items[$q_id] eq $items[$s_id];####!!!
	if ($items[$iden_id] >= $iden_cutoff and $items[$q_cov_id] >= $cov_cutoff and $items[$s_cov_id] >= $cov_cutoff ) {#
		if (not exists $name2clst{ $items[ $q_id ] } and not exists $name2clst{ $items[ $s_id ] } ){
			#add the two seqs to a same new clst
			$name2clst{$items[ $q_id ]}  = $clstmaxidx + 1;
			$name2clst{$items[ $s_id ]}  = $clstmaxidx + 1;
			@{ $clst[ $clstmaxidx + 1] } = ($items[$q_id], $items[$s_id]);#the first two items
			$coc[ $cocmaxidx + 1]->[0]   = $clstmaxidx + 1;
			$clst2coc[ $clstmaxidx + 1 ] = $cocmaxidx  + 1;
			$clstmaxidx ++;
			$cocmaxidx  ++;
		}elsif (exists $name2clst{ $items[ $q_id ] } and not exists $name2clst{ $items[ $s_id] } ) {
			#add the s_id seq to the clust of the q_id seq
			$name2clst{$items[ $s_id ]}  = $name2clst{ $items[ $q_id ] };
			push @{ $clst[ $name2clst{ $items[ $s_id ] } ] }, $items[ $s_id ];
		}elsif (not exists $name2clst{ $items[ $q_id ] } and exists $name2clst{ $items[ $s_id ] } ) {
			#add the q_id seq to the clust of the s_id seq
			$name2clst{ $items[ $q_id ]} = $name2clst{ $items[ $s_id ] };
			push @{ $clst[ $name2clst{ $items[ $q_id ] } ] }, $items[ $q_id ];
		}elsif (exists $name2clst{ $items[ $q_id ] } and exists $name2clst{ $items[ $s_id ] } ) {
			if ($name2clst{ $items[ $q_id ] } == $name2clst{ $items[ $s_id ] } ){#already in the same clust
				#do nothing;
			}elsif( $name2clst{ $items[ $q_id ] } != $name2clst{ $items[ $s_id ] } 
				and   $clst2coc[ $name2clst{ $items[ $q_id ] } ] == $clst2coc[ $name2clst{ $items[ $s_id ] } ] ) {#already in the same coc
				#do nothing;
			}else{#not in the same coc
				#merge the two coc into one coc
				my $cocid1 = $clst2coc[ $name2clst{ $items[ $q_id ] } ];
				my $cocid2 = $clst2coc[ $name2clst{ $items[ $s_id ] } ];
				if ($cocid1 < $cocid2) {
					foreach my $i ( @{ $coc[ $cocid2] } ){
						$clst2coc[ $i ] = $cocid1;
					}
					push @{ $coc[ $cocid1 ] }, @{ $coc[ $cocid2 ] };
					@{$coc[ $cocid2 ] } = ();
				}else{
					foreach my $i ( @{ $coc[ $cocid1 ] } ){
						$clst2coc[ $i ] = $cocid2;
					}
					push @{ $coc[ $cocid2 ] }, @{ $coc[ $cocid1 ] };
					@{$coc[ $cocid1 ] } = ();
				}
			}#
		}else{
			die "should never be here!\n";
		}

	}#if
}#
close IN;


#cluster
my $num_coc = 0;
foreach my $i (0 .. $#coc){
	if ($#{ $coc[$i] } >= 0){
		$num_coc ++;
	}
}

my $num_singleton_seq = 0;
@items = keys %name2clst;
$num_singleton_seq = $seqnum - (scalar @items);

#output clusters and seq of representative of each cluster
open OUT, ">$outfile1" or die "fail to open $outfile1:$!\n";
open OUT2, ">$outfile2" or die "fail to open $outfile2:$!\n";
print OUT "#cov_cutoff:\t",$cov_cutoff,"\n";
print OUT "#iden_cutoff:\t",$iden_cutoff,"\n";
print OUT "#seq_num:\t",$seqnum,"\n";
#print OUT "#clust_num:\t",scalar @clust,"\n";
print OUT "#clust_num:\t",$num_coc + $num_singleton_seq,"\n";
print OUT "#singleton:\t",$num_singleton_seq,"\n";
print "cov_cutoff:\t",$cov_cutoff,"\n";
print "iden_cutoff:\t",$iden_cutoff,"\n";
print "seq_num:\t",$seqnum,"\n";
#print "clust_num:\t",scalar @clust,"\n";
print "#clust_num:\t",$num_coc + $num_singleton_seq,"\n";
print "#singleton:\t",$num_singleton_seq,"\n";

print OUT "#id\tsize\tlist\n";
my $cnt;
my $newid = -1;
foreach my $i ( 0 .. $#coc ) {#
	next if $#{$coc[$i]} < 0;
	$newid ++;
	$cnt = 0;
	foreach my $j ( 0 .. $#{ $coc[ $i ] } ) {
		$cnt += (scalar @{ $clst[ $coc[ $i ]->[ $j ] ] } );
	}
	print OUT $newid,"\t", $cnt;
	foreach my $j ( 0 .. $#{ $coc[ $i ] } ) {
		foreach my $k ( 0 .. $#{ $clst[ $coc[ $i ]->[ $j ] ] } ) {
			print OUT "\t", $clst[ $coc[$i]->[$j] ]->[$k];
		}
	}
	print OUT "\n";
	print OUT2 ">",$seqname[ $name2idx{ $clst[ $coc[$i]->[0] ]->[ 0 ] } ],"\n",$seq[ $name2idx{ $clst[ $coc[$i]->[0] ]->[ 0 ] } ],"\n";
	
}#

#then, output the singleton
foreach my $i (0 .. $#seqname){
	next if exists $name2clst{ $seqname[$i] };
	print OUT $newid,"\t",1,"\t",$seqname[ $i ],"\n";
	print OUT2 ">",$seqname[$i],"\n",$seq[$i],"\n";
	$newid++;
}

close OUT;
close OUT2;



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

