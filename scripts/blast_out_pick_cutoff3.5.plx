#!/usr/bin/perl
# blast_out_pick_cutoff3.5.plx
use strict;
use warnings;


#  This script is used to process the BLAST output to produce clean results
#  that could be used by blast_wrapper.plx.


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


my $qcov_cut= shift;
my $scov_cut= shift;
my $iden_cut= shift;
my $infile = shift;
if (not $infile){
	die "usage: perl blast_out_pick_cutoff3.5.plx <qcov_cutoff> <scov_cutoff> <iden_cut> <infile>\n";
}#if
my $outfile = $infile . ".$qcov_cut.$scov_cut.$iden_cut.pick3.5.tab";
my $outfile2= $infile . ".$qcov_cut.$scov_cut.$iden_cut.discard3.5.tab";

my $col_max = 0;


my $q_len = -1;
my $s_len = -1;
my $seg_num = -1;
my @seg_rec_lst = ();
my @seg_parse_lst = ();


my $last_q;
my $last_s;

my @tmpparse;
my $tmpiden_len;
my $outstr;


my $flag = 0;
my $rec_max_id;
my $tmp_q_start;
my $tmp_q_end;
my $tmp_s_start;
my $tmp_s_end;

my $lno = 0;
my $lbuf;
my @items;


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
$seg_num = 1;
$items[ $s_cov_id ] = int (($items[ $s_end_id ] - $items[ $s_start_id ] + 1) / $items[ $s_len_id ] * 100);
$lbuf .= ("\t" . "$items[ $s_cov_id ]");
$seg_rec_lst[0] = $lbuf;
@{$seg_parse_lst[0]} = @items;
$q_len  = $items[$q_len_id];
$s_len  = $items[$s_len_id];
$last_q = $items[$q_id];
$last_s = $items[$s_id];

while ($lbuf = <IN>){
##!	$lno ++;
	chomp $lbuf;
##!	next if not $lbuf;
##!	next if $lbuf =~ /^#/;
#	$lbuf = $_;
	@items = split /\t/,$lbuf;
##!	if ($col_max != $#items){
##!		die "wrong format: at line $lno in file $infile!\nPlease check the number of columns\n";
##!	}#elsif
	if ($items[ $q_id ] ne $last_q or $items[$s_id] ne $last_s){#this is a new query, process and save last query related records
		#save last q-s related records
			if ( $seg_num == 1) {
				if ( $seg_parse_lst[0]->[ $q_cov_id ] >= $qcov_cut 
				and $seg_parse_lst[0]->[ $s_cov_id ] >= $scov_cut
				and $seg_parse_lst[0]->[ $iden_id ] >= $iden_cut ){
					print OUT $seg_rec_lst[0],"\n";
				}#
			}else{
				#process if there are more than 1 seg
				# @tmpparse for output format
				$seg_parse_lst[0]->[ $q_cov_id ] = ( $seg_parse_lst[0]->[ $q_end_id ] - $seg_parse_lst[0]->[$q_start_id] + 1 ) / $seg_parse_lst[0]->[ $q_len_id ] * 100;

				$seg_parse_lst[0]->[ $s_cov_id ] = ( $seg_parse_lst[0]->[ $s_end_id ] - $seg_parse_lst[0]->[$s_start_id] + 1 ) / $seg_parse_lst[0]->[ $s_len_id ] * 100;


		
				@tmpparse = @{ $seg_parse_lst[0] };
				$tmpiden_len  = $tmpparse[ $iden_id ] * $tmpparse[ $len_id ];
				foreach my $k ( 1 .. $seg_num - 1){#
					$seg_parse_lst[$k]->[ $q_cov_id ] = ( $seg_parse_lst[$k]->[ $q_end_id ] - $seg_parse_lst[$k]->[ $q_start_id] + 1) / $seg_parse_lst[$k]->[ $q_len_id ] * 100;

					$seg_parse_lst[$k]->[ $s_cov_id ] = ( $seg_parse_lst[$k]->[ $s_end_id ] - $seg_parse_lst[$k]->[ $s_start_id] + 1) / $seg_parse_lst[$k]->[ $s_len_id ] * 100;

					$tmpparse[ $q_cov_id ] += $seg_parse_lst[ $k ]->[ $q_cov_id ];
					$tmpparse[ $s_cov_id ] += $seg_parse_lst[ $k ]->[ $s_cov_id ];
					$tmpparse[ $len_id ]   += $seg_parse_lst[ $k ]->[ $len_id ];
					$tmpparse[ $mism_id ]  += $seg_parse_lst[ $k ]->[ $mism_id ];
					$tmpparse[ $gpop_id ]  += $seg_parse_lst[ $k ]->[ $gpop_id ];			
					$tmpparse[ $score_id ] += $seg_parse_lst[ $k ]->[ $score_id ];	
					$tmpparse[ $q_start_id ] .= ("|" . $seg_parse_lst[ $k ]->[ $q_start_id ] );
					$tmpparse[ $q_end_id ]   .= ("|" . $seg_parse_lst[ $k ]->[ $q_end_id ] );
					$tmpparse[ $s_start_id ] .= ("|" . $seg_parse_lst[ $k ]->[ $s_start_id ] );
					$tmpparse[ $s_end_id ]   .= ("|" . $seg_parse_lst[ $k ]->[ $s_end_id ] );

					$tmpiden_len += $seg_parse_lst[ $k ]->[ $iden_id ] * $seg_parse_lst[ $k ]->[ $len_id ];

				}#foreach my $k
				$tmpparse[ $iden_id ] = $tmpiden_len / $tmpparse[ $len_id ] ;

				if ($tmpparse[ $q_cov_id ] >= $qcov_cut 
				and $tmpparse[ $s_cov_id ] >= $scov_cut
				and $tmpparse[ $iden_id  ] >= $iden_cut ) {		
					$outstr = join "\t", @tmpparse;
					print OUT $outstr,"\n";
					foreach my $k ( 0 .. $seg_num - 1){#
						print OUT "#",$seg_rec_lst[$k],"\n";
					}#foreach my $k
				}
			}#if-else
		
		#save last q,s
		$last_q = $items[ $q_id ];
		$last_s = $items[ $s_id ];

		$items[ $s_cov_id ] = int (($items[ $s_end_id ] - $items[ $s_start_id ] + 1) / $items[ $s_len_id ] * 100);
		$lbuf .= ("\t" . "$items[ $s_cov_id ]");
		$seg_rec_lst[0] = $lbuf;
		@{ $seg_parse_lst[0] } = @items;
		$seg_num = 1;


	}else{
	
	$items[ $s_cov_id ] = int (($items[ $s_end_id ] - $items[ $s_start_id ] + 1) / $items[ $s_len_id ] * 100);
	$lbuf .= ("\t" . "$items[ $s_cov_id ]");
		# check whether this record is compatible with existing records
		$flag = 0;#0: compatible; 1: contradict
		$rec_max_id = $seg_num - 1;
		foreach my $i (0 .. $rec_max_id){
			$tmp_q_start = $seg_parse_lst[ $i ]->[ $q_start_id ];
			$tmp_q_end   = $seg_parse_lst[ $i ]->[ $q_end_id ];
			$tmp_s_start = $seg_parse_lst[ $i ]->[ $s_start_id ];
			$tmp_s_end   = $seg_parse_lst[ $i ]->[ $s_end_id ];

			if ( ($items[ $q_start_id ] >= $tmp_q_start and $items[ $q_start_id ] <= $tmp_q_end)
			or   ($items[ $q_end_id   ] >= $tmp_q_start and $items[ $q_end_id   ] <= $tmp_q_end)
			or   ($tmp_q_start >= $items[ $q_start_id ] and $tmp_q_start <= $items[ $q_end_id ])
			or   ($tmp_q_end   >= $items[ $q_start_id ] and $tmp_q_end   <= $items[ $q_end_id ])
			or   ($items[ $s_start_id ] >= $tmp_s_start and $items[ $s_start_id ] <= $tmp_s_end)
			or   ($items[ $s_end_id   ] >= $tmp_s_start and $items[ $s_end_id   ] <= $tmp_s_end)
			or   ($tmp_s_start >= $items[ $s_start_id ] and $tmp_s_start <= $items[ $s_end_id ])
			or   ($tmp_s_end   >= $items[ $s_start_id ] and $tmp_s_end   <= $items[ $s_end_id ]) ){#overlap
				$flag ++;
				last;
			}#
		}#foreach my $i

		if ( $flag ) {
#			push @record_trash, $lbuf;
			print OUT2 $lbuf, "\n";
		}else{
			$seg_rec_lst[$seg_num] = $lbuf;
			@{ $seg_parse_lst[ $rec_max_id + 1 ] } = @items;
			$seg_num ++;
		}#else
	}#else already exists
	
}#while
close IN;

# save the last part
	if ( $seg_num == 1 ) {
			if ( $seg_parse_lst[0]->[ $q_cov_id ] >= $qcov_cut 
				and $seg_parse_lst[0]->[ $s_cov_id ] >= $scov_cut
				and $seg_parse_lst[0]->[ $iden_id ] >= $iden_cut ){
				print OUT $seg_rec_lst[0],"\n";
			}

	}else{#
		#process if there are more than 1 seg
		# @tmpparse for output format
		$seg_parse_lst[0]->[ $q_cov_id ] = ( $seg_parse_lst[0]->[ $q_end_id ] - $seg_parse_lst[0]->[ $q_start_id ] + 1 ) / $seg_parse_lst[0]->[ $q_len_id ] * 100;

		$seg_parse_lst[0]->[ $s_cov_id ] = ( $seg_parse_lst[0]->[ $s_end_id ] - $seg_parse_lst[0]->[ $s_start_id ] + 1 ) / $seg_parse_lst[0]->[ $s_len_id ] * 100;



		@tmpparse = @{ $seg_parse_lst[0] };
		$tmpiden_len  = $tmpparse[ $iden_id ] * $tmpparse[ $len_id ];

		foreach my $k ( 1 .. $seg_num - 1 ) {#
			$seg_parse_lst[$k]->[ $q_cov_id ] = ( $seg_parse_lst[$k]->[ $q_end_id ] - $seg_parse_lst[$k]->[ $q_start_id] + 1)/ $seg_parse_lst[$k]->[$q_len_id] * 100;
			
			$seg_parse_lst[$k]->[ $s_cov_id ] = ( $seg_parse_lst[$k]->[ $s_end_id ] - $seg_parse_lst[$k]->[ $s_start_id] + 1)/ $seg_parse_lst[$k]->[$s_len_id] * 100;

			$tmpparse[ $q_cov_id ] += $seg_parse_lst[ $k ]->[ $q_cov_id ];
			$tmpparse[ $s_cov_id ] += $seg_parse_lst[ $k ]->[ $s_cov_id ];
			$tmpparse[ $len_id ]   += $seg_parse_lst[ $k ]->[ $len_id ];
			$tmpparse[ $mism_id ]  += $seg_parse_lst[ $k ]->[ $mism_id ];
			$tmpparse[ $gpop_id ]  += $seg_parse_lst[ $k ]->[ $gpop_id ];			
			$tmpparse[ $score_id ] += $seg_parse_lst[ $k ]->[ $score_id ];	
			$tmpparse[ $q_start_id ] .= ("|" . $seg_parse_lst[ $k ]->[ $q_start_id ] );
			$tmpparse[ $q_end_id ]   .= ("|" . $seg_parse_lst[ $k ]->[ $q_end_id ] );
			$tmpparse[ $s_start_id ] .= ("|" . $seg_parse_lst[ $k ]->[ $s_start_id ] );
			$tmpparse[ $s_end_id ]   .= ("|" . $seg_parse_lst[ $k ]->[ $s_end_id ] );
			$tmpiden_len += $seg_parse_lst[ $k ]->[ $iden_id ] * $seg_parse_lst[ $k ]->[ $len_id ];


		}#foreach my $k
		$tmpparse[ $iden_id ] = $tmpiden_len / $tmpparse[ $len_id ] ;

		if ($tmpparse[ $q_cov_id ] >= $qcov_cut 
			and $tmpparse[ $s_cov_id ] >= $scov_cut
			and $tmpparse[ $iden_id  ] >= $iden_cut ) {		

			$outstr = join "\t", @tmpparse;
			print OUT $outstr,"\n";

			foreach my $k ( 0 .. $seg_num - 1){#
				print OUT "#",$seg_rec_lst[$k],"\n";
			}#foreach my $k
		}
	}#if-else


close OUT;
close OUT2;

