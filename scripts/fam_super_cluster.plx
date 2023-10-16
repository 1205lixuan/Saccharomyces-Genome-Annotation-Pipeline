#!/usr/bin/perl
# fam_super_cluster_queue.plx
use strict;
use warnings;

#!!!!!!this program works with the script fam_super_diamond_wrapper.plx and blast_out_pick_cutoff3.5.plx

sub GetSeqFromFasta2($$$);#$file,$offset,$seqname; return $seq
sub GetSeqFromFasta($$);#$file,$seqname; return $seq
sub ReadFasta($\@\@);#$filename,@seqname,@seq

my $script_dir = "./scripts";

$|=1;
#my $max_job = 25;

my $produce_seq_file = 1;
my $produce_repseq_file = 1;
#my $diamond_search_rep = 0;
my $diamond_search_fam = 1;

my $diamond_repseq_cpu = 120;
my $diamond_repseq_evalue = 0.001;
my $diamond_fam_cpu = 8;
my $diamond_fam_evalue = 0.001;

#my $cpu = 64;
#my $min_seq_num = 5;	#used to calculate conserv
#my $max_seq_num = 10000;	#do not calculate the identity for those with too much seqs
my $new_flag = shift;#0: restart previous jobs; 1: new jobs
my $fam_file = shift;
my $str2fa_file  = shift;
my $fa_dir   = shift;
my $msa_dir  = "fam_seq";
my $blast_dir = "tmp_blast";
my $iden = shift;
my $cov  = shift;
my $out_prefix = shift;
if (not $out_prefix){
	die "usage: perl fam_super_cluster_queue.plx <restart/new:0/1> <fam_file> <strainname2seqfilename_file> <fa_dir> <iden:0-100> <cov:0-100> <out_pref>\n";
}#if
my $outfile = $out_prefix . ".tab";
my $logfile = $out_prefix . ".log";
my $clustfile = $out_prefix . ".clust";
my $clustfile2 = $out_prefix . ".clust2";

my $repseqfilename = $out_prefix . ".rep.fa";

my @strain_name;	#[$id] =>str name
my %strname2idx;	#{$str_name} => idx in @strain_name
my %strain_name2file;#{str_name} =>prot file name

my %seqname2file;	#{$seqname} => seq file name
#my %seqname2famname;	#{$seqname} => fam name

my %famsize;		#{$famname} => num of seq in the fam, counted based on seqname
my %famsize2;		#{$famname} => num of seq in the fam, seq saved in the fa file

my @famname;		#[$idx] => famname
my %famname2seqfile;	#{$famname} => seq file (unaligned)
#my %famname2msafile;	#{$famname} => msa file
#my %famname2conservfile;

my @seqfilename;
#my @msafilename;
#my @conservfilename;

my %famname2repseqname;#{$famname} => rep_seq_name
my %repseqname2famname;#{$rep_seq_name} => $fam_name

my $num_fam_pair_to_blast = 0;
my $num_fam_pair_blast_run = 0;
my @tmpname;
my @tmpseq;
my @tmpiden;
my @items;
my @items2;
my $lbuf;
my $tmpstr;

my @tmparr;
my $tmpline;

my $cnt1;
my $cnt2;


my %str_seqname_file_pos;
my $curr;

my $tmplen;
my $maxlen;
my $maxid;

my $repseqnum = 0;

open LOG,">$logfile" or die "fail to open $logfile:$!\n";


if (not -e $fa_dir){
	#system ("mkdir $fa_dir");
	die "seq dir $fa_dir does not exists.\n";
}#if
if (not -e $msa_dir){
	system ("mkdir $msa_dir");
}#if
if (not -e $blast_dir){
	system ("mkdir $blast_dir");
}
#read str2fa_file
open IN, $str2fa_file or die "fail to open $str2fa_file:$!\n";
while(<IN>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/,$_;
	$strain_name2file{ $items[0] } = $items[1];
	
	if ( $produce_seq_file ) {
		#create index for fasta file
		if (not -e "$fa_dir/$items[1].mfaidx"){
			open IN2, "$fa_dir/$items[1]" or die "fail to open $fa_dir/$items[1]:$!\n";
			open OUT, ">$fa_dir/$items[1].mfaidx" or die "fail to open $fa_dir/$items[1].mfaidx:$!\n";
			$curr = 0;
			while (<IN2>){
				if (/^>(.+)$/){#this is a seqname
					$curr = ((tell IN2) - length $_);
					print OUT $1,"\t",$curr,"\n";
					$str_seqname_file_pos{ $items[0] }->{$1} = $curr;
				}#
			}#
			close IN2;
			close OUT;
		}else{#already exists, just read
			open IN2, "$fa_dir/$items[1].mfaidx" or die "fail to open $fa_dir/$items[1].mfaidx:$!\n";
			while (<IN2>){
				chomp;
				next if not $_;
				@items2 = split /\t/,$_;
				$str_seqname_file_pos{$items[0] }->{ $items2[0] } = $items2[1];
			}#
			close IN2;
		}#else
	}#if
}
close IN;

#read fam file, extract seqs, pick one representative seq for each fam
#print OUT "#fam\tseqnum\taln_len\tiden_mean\tiden_std\tiden_median\tiden_min\tiden_max\tiden_1_4\tiden_3_4\n";
if ($produce_repseq_file){
	open REP, ">$repseqfilename" or die "fail to open $repseqfilename:$!\n";
}

open IN, $fam_file or die "fail to open $fam_file:$!\n";
$lbuf = <IN>;
chomp $lbuf;
@strain_name = split /\t/,$lbuf;
foreach my $i (1 ..$#strain_name){##the first item is fam name/id
	$strname2idx{ $strain_name[$i] } = $i;
}#
while(<IN>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/,$_;
	push @famname, $items[0];
	$famname2seqfile{ $items[0] } = $items[0] . ".seq.fa";
	@tmpname = ();
	@tmpseq  = ();

	foreach my $i ( 1 .. $#items){
		if ($items[$i]){
			@items2 = split /\|/, $items[$i];
			foreach my $j (0 .. $#items2){
				$famsize{ $items[0] } ++;
				if ( $produce_seq_file ) {
					$tmpstr = GetSeqFromFasta2( "$fa_dir/$strain_name2file{ $strain_name[$i] }",$str_seqname_file_pos{ $strain_name[$i] }->{$items2[$j]}, $items2[$j]);
					if ($tmpstr){
						push @tmpseq, $tmpstr;
						push @tmpname, $items2[$j];
					}else{
						print "warnings: fail to find seq from $fa_dir/$strain_name2file{ $strain_name[$i] }.\n";
					}#else
				}#if
			}
		}#if
	}#

	if ( $produce_seq_file ) {
		open OUTFA, ">$msa_dir/$famname2seqfile{ $items[0] }" or die "fail to open $msa_dir/$famname2seqfile{ $items[0] }:$!\n";
		foreach my $i (0 .. $#tmpname){
			print OUTFA ">",$tmpname[$i],"\n",$tmpseq[$i],"\n";
		}#
		close OUTFA;
	}#if
	

	#pick representative seq for the fam
	if ( not $produce_seq_file ){#@tmpname,@tmpseq are not ready
		ReadFasta "$msa_dir/$famname2seqfile{ $items[0] }", @tmpname, @tmpseq;
	}#
	$tmplen = length $tmpseq[0];
	$maxid  = 0;
	$maxlen = $tmplen;
	foreach my $i (1 .. $#tmpseq){
		$tmplen = length $tmpseq[$i];
		if ($tmplen > $maxlen){
			$maxlen = $tmplen;
			$maxid  = $i;
		}
	}
	$famname2repseqname{ $items[0] } = $tmpname[ $maxid ];
	$repseqname2famname{ $tmpname[ $maxid ] } = $items[0];
	if ($produce_repseq_file){
		print REP ">",$tmpname[ $maxid ], "\n", $tmpseq[ $maxid ],"\n";
	}#if
	$repseqnum ++;
	
	if ($famsize{ $items[0] } > 0 ) {#more than one seq
		push @seqfilename, "$msa_dir/$famname2seqfile{ $items[0] }";
	}#if
}#
close IN;
if ($produce_repseq_file){
	close REP;
}#

#diamond blast repseqfile

my $tmp_cov_cutoff = $cov - 10;
my $tmp_iden_cutoff = $iden - 10;
if ($new_flag or ( (not $new_flag) and (not -e "$out_prefix.rep_diamondblastp.$tmp_cov_cutoff.$tmp_cov_cutoff.$tmp_iden_cutoff.pick3.5.tab") ) ) {
	system ("diamond makedb --in $repseqfilename -d $repseqfilename > $out_prefix.screenout 2>&1");
	system ("diamond blastp -d $repseqfilename -q $repseqfilename -o $out_prefix.rep_diamondblastp -p $diamond_repseq_cpu --very-sensitive -e $diamond_repseq_evalue -k $repseqnum --max-hsps 0 -f 6 qseqid sseqid pident qcovhsp qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore >> $out_prefix.screenout 2>&1");
	system ("perl $script_dir/blast_out_pick_cutoff3.5.plx $tmp_cov_cutoff $tmp_cov_cutoff $tmp_iden_cutoff $out_prefix.rep_diamondblastp");
}

#read back 
my %fam_fam_stat;#{famname}->{fam_name} => stat: 0, not related; 1, related, to merge
open IN, "$out_prefix.rep_diamondblastp.$tmp_cov_cutoff.$tmp_cov_cutoff.$tmp_iden_cutoff.pick3.5.tab" or die "fail to open $out_prefix.rep_diamondblastp.$tmp_cov_cutoff.$tmp_cov_cutoff.$tmp_iden_cutoff.pick3.5.tab:$!\n";
while (<IN>){
	chomp;
	next if not $_;#should never be here
	next if /^#/;#should never be here
	@items = split /\t/, $_;
	next if $items[0] eq $items[1];
	if (not exists $fam_fam_stat{ $repseqname2famname{ $items[0] } }->{ $repseqname2famname{ $items[1] } }
		and not exists $fam_fam_stat{ $repseqname2famname{ $items[1] } }->{ $repseqname2famname{ $items[0] } }){
			$num_fam_pair_to_blast ++;
	}#
	$fam_fam_stat{ $repseqname2famname{ $items[0] } }->{ $repseqname2famname{ $items[1] } } = 0;
}
close IN;

#diamond blastp
print LOG "num_fam_pair_to_blast:\t",$num_fam_pair_to_blast,"\n";
print  "num_fam_pair_to_blast:\t",$num_fam_pair_to_blast,"\n";

	my $tmpname1;
	my $tmpname2;
	my @queue;	#[] => status, 0/unset, not started, 1, running, 2, completed, 3, problem, 4, skip this job.
	my $queue_len = 0;
my @queue_db_file;
my @queue_query_file;
my @queue_db_fam;
my @queue_query_fam;
my @queue_out_file;
my @queue_out_pref;

	my @signal;	#[] => each job has a UNIQUE singal, used as a file name to flag the job is running
	my %running = ();	#{$idx_of_running_job}
	my $num_running = 0;
	my $num_completed = 0;
	my $num_problem = 0;	#treated as job stopped. do not restart.
	my $num_not_started = 0;
#	print "total fam: ", $num_not_started,"\n";

foreach my $i (0 .. $#famname){
	#to save time, i as db
	if ($new_flag or not -e "$msa_dir/$famname2seqfile{ $famname[$i] }.dmnd"){
			system ("diamond makedb --in $msa_dir/$famname2seqfile{ $famname[$i] } -d $msa_dir/$famname2seqfile{ $famname[$i] } >> $out_prefix.screenout 2>&1");
	}
	foreach my $j ($i+1 .. $#famname){
		if(exists $fam_fam_stat{ $famname[$i] }->{ $famname[$j] }
		or exists $fam_fam_stat{ $famname[$j] }->{ $famname[$i] }){
			#add to queue

			#system("perl $script_dir/fam_super_diamond_wrapper2.plx $signal[$i] $diamond_fam_cpu $diamond_fam_evalue $cov $cov $iden 10000 $msa_dir/$famname2seqfile{ $famname[$i] } $msa_dir/$famname2seqfile{ $famname[$j] } $blast_dir/dq.$famname[$i].$famname[$j] > $out_prefix.tmp_screenout.$num_running 2>&1 &");


			system ("diamond blastp -d $msa_dir/$famname2seqfile{ $famname[$i] } -q $msa_dir/$famname2seqfile{ $famname[$j] } -o $blast_dir/dq.$famname[$i].$famname[$j].diamondblastp -p $diamond_fam_cpu --very-sensitive -e $diamond_fam_evalue -k 10000 --max-hsps 0 -f 6 qseqid sseqid pident qcovhsp qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore >> $out_prefix.screenout 2>&1");
			system ("perl $script_dir/blast_out_pick_cutoff3.5.plx $cov $cov $iden $blast_dir/dq.$famname[$i].$famname[$j].diamondblastp");


			open IN, "$blast_dir/dq.$famname[$i].$famname[$j].diamondblastp.$cov.$cov.$iden.pick3.5.tab" or die "fail to open $blast_dir/dq.$famname[$i].$famname[$j].diamondblastp.$cov.$cov.$iden.pick3.5.tab:$!\n";
			while (<IN>){
				chomp;
				next if not $_;
				next if /^#/;
				$fam_fam_stat{ $famname[$i] }->{ $famname[$j] } ++;
			}
			close IN;
			if( $fam_fam_stat{ $famname[$i] }->{ $famname[$j] } ) {
				print LOG "PAIR_FOUND_HIT:","\t",$famname[$i],"\t",$famname[$j],"\t",$fam_fam_stat{ $famname[$i] }->{ $famname[$j] },"\n";
			}#if


		}
	}#
}



#single-linkage

my %name2clst;#{$famname} => idx of the corresponding cluster
my @clst;#[idx]->[idx] => name
my $clstmaxidx = -1;
my @coc;#[idx]->[idx] => clust idx; cluster_of_cluster
my $cocmaxidx  = -1;
my @clst2coc;#[idx,i.e. idx of @clst] => idx of @coc


foreach my $i (keys %fam_fam_stat){
	foreach my $j (keys %{$fam_fam_stat{$i}}){
		next if $i eq $j;
		next if not $fam_fam_stat{$i}->{$j};
		if (not exists $name2clst{ $i } and not exists $name2clst{ $j } ){
			#add the two fams to a same new clst
			$name2clst{$i}  = $clstmaxidx + 1;
			$name2clst{$j}  = $clstmaxidx + 1;
			@{ $clst[ $clstmaxidx + 1] } = ($i, $j);#the first two items
			$coc[ $cocmaxidx + 1]->[0]   = $clstmaxidx + 1;
			$clst2coc[ $clstmaxidx + 1 ] = $cocmaxidx  + 1;
			$clstmaxidx ++;
			$cocmaxidx  ++;
		}elsif (exists $name2clst{ $i } and not exists $name2clst{ $j } ) {
			#add the j fam to the clust of the i fam
			$name2clst{$j}  = $name2clst{ $i };
			push @{ $clst[ $name2clst{ $j } ] }, $j;
		}elsif (not exists $name2clst{ $i } and exists $name2clst{ $j } ) {
			#add the i fam to the clust of the j fam
			$name2clst{ $i } = $name2clst{ $j };
			push @{ $clst[ $name2clst{ $i } ] }, $i;
		}elsif (exists $name2clst{ $i } and exists $name2clst{ $j } ) {
			if ($name2clst{ $i } == $name2clst{ $j } ){#already in the same clust
				#do nothing;
			}elsif( $name2clst{ $i } != $name2clst{ $j } 
				and   $clst2coc[ $name2clst{ $i } ] == $clst2coc[ $name2clst{ $j } ] ) {#already in the same coc
				#do nothing;
			}else{#not in the same coc
				#merge the two coc into one coc
				my $cocid1 = $clst2coc[ $name2clst{ $i } ];
				my $cocid2 = $clst2coc[ $name2clst{ $j } ];
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
		
	}
}


my $num_coc = 0;
foreach my $i (0 .. $#coc){
	if ($#{ $coc[$i] } >= 0){
		$num_coc ++;
	}
}

my $num_fam_wo_search = 0;
@items = keys %name2clst;
$num_fam_wo_search = (scalar @famname) - (scalar @items);

print "famname:\t",scalar @famname,"\n";

my @clustered_famname;
#output clusters and seq of representative of each cluster
open OUT, ">$clustfile" or die "fail to open $clustfile:$!\n";

print OUT "#cov_cutoff:\t",$cov,"\n";
print OUT "#iden_cutoff:\t",$iden,"\n";
print OUT "#fam_num:\t",$repseqnum,"\n";#@famname
#print OUT "#clust_num:\t",scalar @clust,"\n";
print OUT "#clust_num:\t",$num_coc + $num_fam_wo_search,"\n";
print OUT "#w/o_search:\t",$num_fam_wo_search,"\n";
print "cov_cutoff:\t",$cov,"\n";
print "iden_cutoff:\t",$iden,"\n";
print "fam_num:\t",$repseqnum,"\n";
#print "clust_num:\t",scalar @clust,"\n";
print "#clust_num:\t",$num_coc + $num_fam_wo_search,"\n";
print "#w/o_search:\t",$num_fam_wo_search,"\n";

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
			push @{$clustered_famname[$newid]}, $clst[ $coc[$i]->[$j] ]->[$k];
		}
	}
	print OUT "\n";
	
}#

#then, output the fam w/o_search
foreach my $i (0 .. $#famname){
	next if exists $name2clst{ $famname[$i] };
	$newid++;
	print OUT $newid,"\t",1,"\t",$famname[ $i ],"\n";
	push @{$clustered_famname[$newid]},$famname[ $i ];

}
close OUT;


my %origin_famid_sp_memb;	#{$famid}->[$spec_idx]->[$para_idx] => prot_name
my %origin_famid_sp_memb_cnt;	#{$famid}->[$spec_idx] => prot member number
my $famtitle;
#produce new fam file
open FOUT, ">$outfile" or die "fail to open $outfile:$!\n";
 
##read orginal fam file again

open IN, $fam_file or die "fail to open $fam_file:$!\n";
$famtitle = <IN>;
chomp $famtitle;
while (<IN>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/,$_,-1;
	foreach my $i (1 .. $#items){
		@items2 = split /\|/,$items[$i];
		@{$origin_famid_sp_memb{ $items[0] }->[$i-1]} = @items2;
		$origin_famid_sp_memb_cnt{ $items[0] }->[$i-1] = scalar @items2;
	}
}
close IN;

##merge fam where necessary
###famname order
my %famnameidx;	#{$famname} => idx in @famname
foreach my $i (0 .. $#famname){
	$famnameidx{ $famname[$i] } = $i;
}

###
my @clustered_famname_rep;	#[$idx_as_in_@clustered_famname] => rep fam name with minimal idx in @famname/%famnameidx;
foreach my $i (0 .. $#clustered_famname){
	@items = sort {$famnameidx{$a} <=> $famnameidx{$b}} @{$clustered_famname[$i]};
	$clustered_famname_rep[ $i ] = $items[0];
}


foreach my $i (0 .. $#clustered_famname){
	print LOG $i,"\t",$clustered_famname_rep[$i],"\t","@{$clustered_famname[$i]}","\n";
}

#my @merge;	#same size as the origin_famid_sp_mem
print "output new fam tab\n";
print FOUT $famtitle,"\n";
foreach my $i (sort { $famnameidx{ $clustered_famname_rep[$a] } <=> $famnameidx{$clustered_famname_rep[$b] } } (0 .. $#clustered_famname_rep) ){
	print $i,"\t",$clustered_famname_rep[$i],"\n";
	if ( $#{$clustered_famname[$i]} == 0) {#only one fam in the cluster, output directly
		print FOUT $clustered_famname_rep[$i];
		foreach my $j ( 0 .. $#{$origin_famid_sp_memb{ $clustered_famname[$i]->[0] } } ){#each sp
			if ( $origin_famid_sp_memb_cnt{ $clustered_famname[$i]->[0] }->[$j] ){
				print FOUT "\t",$origin_famid_sp_memb{$clustered_famname[$i]->[0] }->[$j]->[0];
				foreach my $k ( 1 .. $#{ $origin_famid_sp_memb{ $clustered_famname[$i]->[0] }->[$j] } ){
					print FOUT "|",$origin_famid_sp_memb{$clustered_famname[$i]->[0] }->[$j]->[$k];
				}
			}else{
				print FOUT "\t";
			}
		}
		print FOUT "\n";
	}else{#more than one fam mem, merge
		
		print FOUT $clustered_famname_rep[$i];
		foreach my $j ( 0 .. $#{ $origin_famid_sp_memb{ $clustered_famname[$i]->[0] } } ){#each sp
			@items = ();
			foreach my $k ( 0 .. $#{$clustered_famname[$i]} ) {#each fam in the cluster
				push @items, @{ $origin_famid_sp_memb{ $clustered_famname[$i]->[$k] }->[$j] };
			}
			if ($#items == -1){
				print FOUT "\t";
			}elsif ( $#items == 0){
				print FOUT "\t",$items[0];
			}else{
				@items = sort {$a cmp $b} @items;
				print FOUT "\t",$items[0];
				foreach my $k (1 .. $#items){
					print FOUT "|",$items[$k];
				}
			}
		}
	print FOUT "\n";
	}#

}
close FOUT;


close LOG;


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
		$seq_r->[$i] =~ s/\*//g;
	}#foreach
	
	#print "ReadFasta: $filename","\t",scalar @{$seqname_r}, "\t", scalar @{$seq_r},"\n";
}#ReadFasta()


	
sub GetSeqFromFasta($$){#$file,$seqname; return $seq
	my $filename = shift;
	my $seqname  = shift;
	my $currseq = "";
	my $flag = 0;
	my $lbuf;
	open GetSeqFromFasta_FH, "$filename" or die "sub GetSeqFromFasta(): fail to open file $filename: $!\n";
	while (<GetSeqFromFasta_FH>) {
#		$lineno ++;
		chomp;
		next if not $_;
		$lbuf = $_;
#		if(! length $_){
#			next;
#		}#blank line
		if(m/^>(.+)$/){#this is the name of a new seq
			if ($1 eq $seqname) {
				$flag ++;#found the target seq
			}else{#ne seqname
				if ($flag){#already got the seq, means the seq of the target has reached
					last;
				}
			}#else
		}else{
			if ($flag){
				$currseq .= $lbuf;
			}#if
		}#else
	}#while
	close GetSeqFromFasta_FH;
	if (not $flag){
		print "GetSeqFromFasta(): warnings -- fail to find seq $seqname in $filename.\n";
	}#if
	return $currseq;
}#sub GetSeqFromFasta($$){#$file,$seqname; return $seq	



	
sub GetSeqFromFasta2($$$){#$file,$offset,$seqname; return $seq
	my $filename = shift;
	my $offset   = shift;
	my $seqname  = shift;
	my $currseq = "";
	my $flag = 0;
	my $lbuf;
	my $lbuf2;
	open GetSeqFromFasta_FH, "$filename" or die "sub GetSeqFromFasta(): fail to open file $filename: $!\n";
	seek GetSeqFromFasta_FH,$offset,0;
	$lbuf = <GetSeqFromFasta_FH>;
	chomp $lbuf;
	if ($lbuf =~ /^>(.+)$/){
		if ($1 eq $seqname){#correct
			while (<GetSeqFromFasta_FH>) {
				chomp;
				next if not $_;
				$lbuf2 = $_;
				if($lbuf2 =~ /^>/){#this is the name of a new seq
					last;
				}else{
					$currseq .= $lbuf2;
				}#else
			}#while

		}#else
	}#if
	close GetSeqFromFasta_FH;
	if (not $currseq){
		print "GetSeqFromFasta(): warnings -- fail to find seq $seqname in $filename.\n";
	}#if
	return $currseq;
}#sub GetSeqFromFasta2($$$){#$file,$offset,$seqname; return $seq	
