#!/usr/bin/perl
# blast_out_pick3.plx
use strict;
use warnings;


#  This script is used to process the BLAST output to produce clean results
#  that could be used by blast_cluster3.plx.

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
my $s_cov_id = 15;	#undef
my $id_max   = 14;	#assume there is no s_cov_id in the input file, only process raw blastp output file



my $infile = shift;
if (not $infile){
	die "usage: perl blast_out_pick3.plx <infile>\n";
}#if
my $outfile = $infile . ".pick3.tab";
my $outfile2= $infile . ".discard3.tab";

my $col_max = 0;

my @record_trash;
my %query_subject;	#{query}->{subject}->{q_len/s_len/seg_num/seg_parse_lst/seg_rec_lst}
					#{query}->{subject}->{seg_parse_lst}->[seg_id]->[item_id] => items
my $last_q;

my @tmpparse;
my $tmpiden_len;
my $outstr;

my $lno = 0;
my $lbuf;
my @items;

#if (-e "$infile.pick3.tab"){
#	die "exit - file $infile.pick3.tab already exists. Don't perform the analysis.\n";
#}#if

open OUT, ">$outfile" or die "fail to open $outfile:$!\n";
open OUT2, ">$outfile2" or die "fail to open $outfile2:$!\n";
# read and quick parse
open IN, "$infile" or die "fail to open $infile:$!\n";
$lbuf = <IN>;
chomp $lbuf;
die "wrong format: the first line is blank at line 0 in $infile!" if not $lbuf;
die "wrong format: the first line is remarked at line 0 in $infile!" if $lbuf =~ /^#/;
@items = split /\t/, $lbuf;
$col_max = $#items;
if ( $id_max > $col_max ) {#
	die "wrong format id_max > col_max ($id_max > $col_max) at line 0 in $infile!\n";
}#if
$query_subject{ $items[$q_id] }->{ $items[$s_id] }->{ seg_num } = 1;
$items[ $s_cov_id ] = int (($items[ $s_end_id ] - $items[ $s_start_id ] + 1) / $items[ $s_len_id ] * 100);
$lbuf .= ("\t" . "$items[ $s_cov_id ]");
$query_subject{ $items[$q_id] }->{ $items[$s_id] }->{ seg_rec_lst }->[0] = $lbuf;
@{ $query_subject{ $items[ $q_id ] }->{ $items[$s_id] }->{ seg_parse_lst }->[ 0 ] } = @items;
$last_q = $items[$q_id];

while (<IN>){
	$lno ++;
	chomp;
	next if not $_;
	next if /^#/;
	$lbuf = $_;
	@items = split /\t/;
	if ($col_max != $#items){
		die "wrong format: at line $lno in file $infile!\nPlease check the number of columns\n";
	}#elsif
	if ($items[ $q_id ] ne $last_q){#this is a new query, process and save last query related records
		#save last q related records
		foreach my $i (keys %query_subject){
			foreach my $j (keys %{$query_subject{$i}} ){
				if ($query_subject{ $i }->{ $j }->{seg_num} == 1){
					print OUT $query_subject{ $i }->{ $j }->{seg_rec_lst}->[0],"\n";
					next;
				}#if
				#process if there are more than 1 seg
				# @tmpparse for output format
				$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $q_cov_id ] =
				 ( ($query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $q_end_id ] -
				$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $q_start_id ] + 1)/
				$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $q_len_id ] *100 );

				$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $s_cov_id ] =
			 ( ($query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $s_end_id ] -
				$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $s_start_id ] + 1)/
				$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $s_len_id ] *100 );
		
				@tmpparse = @{ $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ] };
				$tmpiden_len  = $tmpparse[ $iden_id ] * $tmpparse[ $len_id ];
				foreach my $k ( 1.. $query_subject{ $i }->{ $j }->{seg_num} - 1){#
					#print "k:$k\n";
					#print "qend:",$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_end_id ],"\n";
					#print "qstart:",$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_start_id ],"\n";
					#print "qlen:",$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_len_id ],"\n";

					$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_cov_id ] =
					 ( ($query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_end_id ] -
					$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_start_id ] + 1)/
					$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_len_id ] *100);

					$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_cov_id ] =
					 ( ($query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_end_id ] -
					$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_start_id ] + 1)/
					$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_len_id ] *100);

					$tmpparse[ $q_cov_id ] += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_cov_id ];
					$tmpparse[ $s_cov_id ] += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_cov_id ];
					$tmpparse[ $len_id ]   += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $len_id ];
					$tmpparse[ $mism_id ]  += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $mism_id ];
					$tmpparse[ $gpop_id ]  += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $gpop_id ];			
					$tmpparse[ $score_id ] += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $score_id ];	
					$tmpparse[ $q_start_id ] .= ("|" . $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_start_id ] );
					$tmpparse[ $q_end_id ]   .= ("|" . $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_end_id ] );
					$tmpparse[ $s_start_id ] .= ("|" . $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_start_id ] );
					$tmpparse[ $s_end_id ]   .= ("|" . $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_end_id ] );

					$tmpiden_len += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $iden_id ] * $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $len_id ];
				}#foreach my $k
				$tmpparse[ $iden_id ] = $tmpiden_len / $tmpparse[ $len_id ] ;
		
				$outstr = join "\t", @tmpparse;
				print OUT $outstr,"\n";
				foreach my $k ( 0.. $query_subject{ $i }->{ $j }->{seg_num} - 1){#
					print OUT "#",$query_subject{ $i }->{ $j }->{ seg_rec_lst }->[ $k ], "\n";
				}#foreach my $k
			}#foreach my $j
		}#foreach my $i

		foreach my $i (@record_trash){
			print OUT2 $i, "\n";
		}#
		
		#reset var
		%query_subject = ();
		@record_trash  = ();
		
		#save last q
		$last_q = $items[ $q_id ];
	}#
	
	$items[ $s_cov_id ] = int (($items[ $s_end_id ] - $items[ $s_start_id ] + 1) / $items[ $s_len_id ] * 100);
	$lbuf .= ("\t" . "$items[ $s_cov_id ]");
	if ( not exists $query_subject{ $items[$q_id] }->{ $items[ $s_id ] } ) {
		$query_subject{ $items[$q_id] }->{ $items[ $s_id ] }->{ seg_rec_lst }->[0] = $lbuf;
		@{ $query_subject{ $items[ $q_id ] }->{ $items[$s_id] }->{ seg_parse_lst }->[ 0 ] } = @items;
		$query_subject{ $items[$q_id] }->{ $items[$s_id] }->{ seg_num } = 1;
	}else{#already exists
		# check whether this record is compatible with existing records
		my $flag = 0;#0: compatible; 1: contradict
		my $rec_max_id = $query_subject{ $items[$q_id] }->{ $items[$s_id] }->{ seg_num } - 1;
		foreach my $i (0 .. $rec_max_id){
			my $tmp_q_start = $query_subject{ $items[$q_id] }->{ $items[$s_id] }->{ seg_parse_lst }->[ $i ]->[ $q_start_id ];
			my $tmp_q_end   = $query_subject{ $items[$q_id] }->{ $items[$s_id] }->{ seg_parse_lst }->[ $i ]->[ $q_end_id ];
			my $tmp_s_start = $query_subject{ $items[$q_id] }->{ $items[$s_id] }->{ seg_parse_lst }->[ $i ]->[ $s_start_id ];
			my $tmp_s_end   = $query_subject{ $items[$q_id] }->{ $items[$s_id] }->{ seg_parse_lst }->[ $i ]->[ $s_end_id ];
			foreach my $j ( $items[ $q_start_id ] .. $items[ $q_end_id ]) {
				if ( $j >= $tmp_q_start and $j <= $tmp_q_end ) {
					$flag ++;
					last;
				}#if
			}#foreach my $j
			last if $flag;
			foreach my $j ( $items[ $s_start_id ] .. $items[ $s_end_id ]){
				if ( $j >= $tmp_s_start and $j <= $tmp_s_end ) {
					$flag ++;
					last;
				}#if	
			}#foreach my $j
			last if $flag;
		}#foreach my $i

		if ( $flag ) {
			push @record_trash, $lbuf;
		}else{
			push @{ $query_subject{ $items[$q_id] }->{ $items[ $s_id ] }->{ seg_rec_lst }}, $lbuf;
			@{ $query_subject{ $items[ $q_id ] }->{ $items[$s_id] }->{ seg_parse_lst }->[ $rec_max_id + 1 ] } = @items;
			$query_subject{ $items[$q_id] }->{ $items[$s_id] }->{ seg_num } ++;
		}#else
	}#else already exists
	
}#while
close IN;

# save the last part
foreach my $i (keys %query_subject){
	foreach my $j (keys %{$query_subject{$i}} ){
		if ($query_subject{ $i }->{ $j }->{seg_num} == 1){
			print OUT $query_subject{ $i }->{ $j }->{seg_rec_lst}->[0],"\n";
			next;
		}#if
		#process if there are more than 1 seg
		# @tmpparse for output format
		$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $q_cov_id ] =
		 ( ($query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $q_end_id ] -
		$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $q_start_id ] + 1)/
		$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $q_len_id ] *100 );
		$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $s_cov_id ] =
	 ( ($query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $s_end_id ] -
		$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $s_start_id ] + 1)/
		$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ]->[ $s_len_id ] *100 );

		@tmpparse = @{ $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ 0 ] };
		$tmpiden_len  = $tmpparse[ $iden_id ] * $tmpparse[ $len_id ];
		foreach my $k ( 1.. $query_subject{ $i }->{ $j }->{seg_num} - 1){#
			#print "k:$k\n";
			#print "qend:",$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_end_id ],"\n";
			#print "qstart:",$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_start_id ],"\n";
			#print "qlen:",$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_len_id ],"\n";
			$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_cov_id ] =
			 ( ($query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_end_id ] -
			$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_start_id ] + 1)/
			$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_len_id ] *100);
			$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_cov_id ] =
			 ( ($query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_end_id ] -
			$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_start_id ] + 1)/
			$query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_len_id ] *100);
			$tmpparse[ $q_cov_id ] += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_cov_id ];
			$tmpparse[ $s_cov_id ] += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_cov_id ];
			$tmpparse[ $len_id ]   += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $len_id ];
			$tmpparse[ $mism_id ]  += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $mism_id ];
			$tmpparse[ $gpop_id ]  += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $gpop_id ];			
			$tmpparse[ $score_id ] += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $score_id ];	
			$tmpparse[ $q_start_id ] .= ("|" . $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_start_id ] );
			$tmpparse[ $q_end_id ]   .= ("|" . $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $q_end_id ] );
			$tmpparse[ $s_start_id ] .= ("|" . $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_start_id ] );
			$tmpparse[ $s_end_id ]   .= ("|" . $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $s_end_id ] );
			$tmpiden_len += $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $iden_id ] * $query_subject{ $i }->{ $j }->{ seg_parse_lst }->[ $k ]->[ $len_id ];
		}#foreach my $k
		$tmpparse[ $iden_id ] = $tmpiden_len / $tmpparse[ $len_id ] ;

		$outstr = join "\t", @tmpparse;
		print OUT $outstr,"\n";
		foreach my $k ( 0.. $query_subject{ $i }->{ $j }->{seg_num} - 1){#
			print OUT "#",$query_subject{ $i }->{ $j }->{ seg_rec_lst }->[ $k ], "\n";
		}#foreach my $k
	}#foreach my $j
}#foreach my $i

foreach my $i (@record_trash){
	print OUT2 $i, "\n";
}#


close OUT;
close OUT2;

