#!/usr/bin/perl
#blastn_annt.plx
use strict;
use warnings;

#################################
#   blastn_annt.plx
#   version 1.1
#   2022.12.26
#   by Bin-Bin Xie
#   email: xbb@sdu.edu.cn
#################################

#  This script is used to call gene models based on blastn matches without indels.
#  It was originally designed to annotate Saccharomyces cerevisiae genomes.
#  The reference genomic sequence file and coordinate file for S. cerevisiae S288c
#  were provided along with this script. Just put the reference files in the dir/folder
#  where you run the annoation job.
#  Theoretically, it can also be used to annotate genomes for other organisms.
#  If you want to do this, your need to prepare the reference genome sequence file and 
#  the associated coordinate file for your reference genome.


sub ReadFasta($\@\@);#$filename,@seqname,@seq
sub WriteFasta($\@\@);#$filename,@seqname,@seq
sub rc($);#$inseq

sub blast2map($$);#$blastout_file,$outmap_file
sub map_annt($$$$);#$tgt_fna,$mapfile,$pref,$dir

my $max_ext = 300;	#max num of codons to extend or shrink; 20210917
my $mit_seq = "NC_001224.1";#for GCF_000146045.2_R64_genomic
my $ref_fa  = "GCF_000146045.2_R64_genomic.shortname.fna";
my $ref_gff = "GCF_000146045.2_R64_genomic.fix2.gff";

my $cpu          = shift;
my $ref_dir      = shift;
my $inputseq_dir = shift;
my $tag          = shift;
my $fafile       = shift;
if (not $fafile){
	die "usage: perl blastn_annt.plx <cpu_num> <ref_dir> <input_seq_dir> <tag> <target_fa>\n";
}#if
my $outfile = $fafile . ".log";
my @seqname = ();
my @seq     = ();
my @items;
my $lbuf;
my $cnt = 0;

if (not -e "$ref_dir/$ref_fa.ndb"){#work with NCBI BLAST+ 2.11.0 (and other versions producing output files of the same format)
	system("makeblastdb -in $ref_dir/$ref_fa -dbtype nucl");
}#if

ReadFasta ("$inputseq_dir/$fafile",@seqname,@seq);
#convert to short name, to work with blast
foreach my $i (0 .. $#seqname){
	if ($seqname[$i] =~ /^(\S+)/){
		$seqname[$i] = $1;
	}#if
}
WriteFasta ("$inputseq_dir/$tag.genome.fna",@seqname,@seq);

#blastn
system("blastn -query $inputseq_dir/$tag.genome.fna -db $ref_dir/$ref_fa -out $inputseq_dir/$tag.genome.fna.blastn -num_threads $cpu");#allow more than one target
#map, generating a coordinate map file
blast2map("$inputseq_dir/$tag.genome.fna.blastn","$inputseq_dir/$tag.genome.fna.blastn.map");
#annt, based on the map file
map_annt("$inputseq_dir/$tag.genome.fna","$inputseq_dir/$tag.genome.fna.blastn.map",$tag,$inputseq_dir);


sub blast2map($$){#$blastout_filename,$outmap_filename
	my $infile = shift;
	my $outfile= shift;
	my $ln1;
	my $ln2;
	my $ln3;
	my $qn;
	my $sn;
	#to save space, do not save the query name
	my %q2s_name = ();	#{$q_pos} => sn
	my %q2s_pos  = ();	#{$q_pos} => pos on s or 0;
	my %q2s_str  = ();	#{$q_pos} => strand

	my $qs;
	my $qe;
	my $ss;
	my $se;
	my $qaln;
	my $saln;

	my @qa;
	my @sa;

	my $s_curr;
	my $q_curr;

	my $s_blk_s;
	my $q_blk_s;

	my $inc;

	my @items;
	my $lno = 0;

	open OUT, ">$outfile" or die "fail to open $outfile:$!\n";
	open IN, $infile or die "fail to open $infile:$!\n";
	while ($ln1 = <IN>){
		$lno ++;
		chomp $ln1;
		if ($ln1 =~ /^Query= (.+)$/){
			#save the info for the last query
			@items = sort { $a <=> $b } keys %q2s_pos;
			#print "|@items|\n";
			if ($#items >= 0 ) {
				#print $i,"\t","|@items|","\n";
				$q_blk_s = $items[0];
				$s_blk_s = $q2s_pos{ $items[0] };
				$sn      = $q2s_name{ $items[0] };
				if ($q2s_str{ $items[0] } eq "+"){
					$inc = 1;
				}else{
					$inc = -1;
				}#else     
				print OUT "SQ","\t",$qn,"\t",$sn,"\t",$q2s_str{$items[0]},"\n";
				foreach my $j (1 .. $#items){
					if ($items[$j] == $items[$j-1] + 1
					and $q2s_name{$items[$j]} eq $q2s_name{$items[$j-1]} 
					and $q2s_pos{$items[$j]} == $q2s_pos{$items[$j-1]} + $inc
					and $q2s_str{$items[$j]} eq $q2s_str{$items[$j-1]} ) {
						#continue
					}else{# save this segment and start a new seg
						print OUT $q_blk_s,"\t",$items[$j-1],"\t",$s_blk_s,"\t",$q2s_pos{$items[$j-1]},"\n";

						$sn = $q2s_name{$items[$j]};
						$q_blk_s = $items[$j];
						$s_blk_s = $q2s_pos{$items[$j]};
						if ($q2s_str{$items[$j]} eq "+"){
							$inc = 1;
						}else{
							$inc = -1;
						}#else

						if ( $q2s_name{$items[$j]} ne $q2s_name{$items[$j-1]} 
						or   $q2s_str{$items[$j]} ne $q2s_str{$items[$j-1]} ) {#a new subject or a new strand
							print OUT "SQ","\t",$qn,"\t",$sn,"\t",$q2s_str{$items[$j]},"\n";
						}#if
					}#else
				}#$j

				#save the last blk
				print OUT $q_blk_s,"\t",$items[-1],"\t",$s_blk_s,"\t",$q2s_pos{$items[-1]},"\n";
			}#

			$qn = $1;
		
			%q2s_name = ();
			%q2s_pos  = ();
			%q2s_str  = ();
		}elsif ($ln1 =~ /^>(.+)$/){
			$sn = $1;
		}elsif ($ln1 =~ /^Query +(\d+) +(\S+) +(\d+)$/){
			$qs = $1;
			$qaln = $2;
			$qe = $3;
			<IN>;
			$lno ++;
			$ln3 = <IN>;
			$lno ++;
			chomp $ln3;
			if ($ln3 =~ /^Sbjct +(\d+) +(\S+) + (\d+)$/){#20210926
				$ss = $1;
				$saln = $2;
				$se = $3;
			}else{
				die "wrong format at line $lno in $infile\n";
			}#else
			#cmp
			$qaln =~ tr/a-z/A-Z/;
			$saln =~ tr/a-z/A-Z/;
			@qa = split //,$qaln;
			@sa = split //,$saln;
			if ($#qa != $#sa){
				die "wrong parsing aln at line $lno in $infile.\n";
			}#
			$s_curr = $ss;
			$q_curr = $qs;

			foreach my $i (0 .. $#qa){
				if ($qa[$i] eq "-"){
					if ($ss < $se) {
						$s_curr ++;	
					}else{
						$s_curr --;
					}#else						
				}elsif ($sa[$i] eq "-"){
					if (not exists $q2s_name{ $q_curr }){
						$q2s_name{ $q_curr } = $sn;
						$q2s_pos{ $q_curr } = 0;
						if ($ss < $se){
							$q2s_str{ $q_curr } = "+";
						}else{
							$q2s_str{ $q_curr } = "-";
						}#else
					}#if
					$q_curr ++;
				}else{#allow mismatch
					if (not exists $q2s_name{$q_curr}){
						$q2s_name{ $q_curr } = $sn;
						$q2s_pos{ $q_curr } = $s_curr;
						if ($ss < $se){
							$q2s_str{ $q_curr } = "+";
						}else{
							$q2s_str{ $q_curr } = "-";
						}#else					
					}#if
					
					$q_curr ++;
					if ($ss < $se) {
						$s_curr ++;	
					}else{
						$s_curr --;
					}#else				
				}#else
			}#$i

		}#elsif
	}#while
	close IN;

	#save the last query
	@items = sort { $a <=> $b } keys %q2s_pos;
	if ($#items >=0){
		$q_blk_s = $items[0];
		$s_blk_s = $q2s_pos{ $items[0] };
		$sn      = $q2s_name{ $items[0] };
		if ($q2s_str{ $items[0] } eq "+"){
			$inc = 1;
		}else{
			$inc = -1;
		}#else     
		print OUT "SQ","\t",$qn,"\t",$sn,"\t",$q2s_str{$items[0]},"\n";
		foreach my $j (1 .. $#items){
			if ($items[$j] == $items[$j-1] + 1
			and $q2s_name{$items[$j]} eq $q2s_name{$items[$j-1]} 
			and $q2s_pos{$items[$j]} == $q2s_pos{$items[$j-1]} + $inc
			and $q2s_str{$items[$j]} eq $q2s_str{$items[$j-1]} ) {
				#continue
			}else{# save this segment and start a new seg
				print OUT $q_blk_s,"\t",$items[$j-1],"\t",$s_blk_s,"\t",$q2s_pos{$items[$j-1]},"\n";

				$sn = $q2s_name{$items[$j]};
				$q_blk_s = $items[$j];
				$s_blk_s = $q2s_pos{$items[$j]};
				if ($q2s_str{$items[$j]} eq "+"){
					$inc = 1;
				}else{
					$inc = -1;
				}#else

				if ( $q2s_name{$items[$j]} ne $q2s_name{$items[$j-1]} 
				or   $q2s_str{$items[$j]} ne $q2s_str{$items[$j-1]} ) {#a new subject or a new strand
					print OUT "SQ","\t",$qn,"\t",$sn,"\t",$q2s_str{$items[$j]},"\n";
				}#if
			}#else
		}#$j

		#save the last blk
		print OUT $q_blk_s,"\t",$items[-1],"\t",$s_blk_s,"\t",$q2s_pos{$items[-1]},"\n";
	}#if
	close OUT;
	print "written to $outfile.\n";
}#sub



sub map_annt($$$$){#$tgt_fna,$mapfile,$pref,$dir
	my $tgt_fna  = shift;
	my $mapfile  = shift;
	my $pref     = shift;
	my $dir      = shift;
	my $outfile  = $dir . "/" . $pref . ".log";
	my $outgff   = $dir . "/" . $pref . ".blastn.gff";
	my $outgene  = $dir . "/" . $pref . ".blastn.gene.fna";
	my $outmrna  = $dir . "/" . $pref . ".blastn.mrna.fna";
	my $outprot  = $dir . "/" . $pref . ".blastn.prot.faa";
	my $outcds   = $dir . "/" . $pref . ".blastn.cds.fna";
	my $outnomdl1= $dir . "/" . $pref . ".blastn.no_mdl_mask.fna";
#	my $outnomdl2= $dir . "/" . $pref . ".blastn.no_mdl_seg.fna";
	my $pseudogff= $dir . "/" . $pref . ".blastn.pseudo.gff";
	my $pseudogene= $dir . "/" . $pref . ".blastn.pseudogene.fna";
	my $pseudomrna= $dir . "/" . $pref . ".blastn.pseudomrna.fna";
	my $pseudoprot= $dir . "/" . $pref . ".blastn.pseudoprot.faa";
	my $pseudocds = $dir . "/" . $pref . ".blastn.pseudocds.fna";

	my $flag;
	my $ch1;
	my $ch2;
	my $method = "";
	my %transl_tab1 = (TTT => "F",
	TTC => "F",
	TTA => "L",
	TTG => "L",
	TCT => "S",
	TCC => "S",
	TCA => "S",
	TCG => "S",
	TAT => "Y",
	TAC => "Y",
	TAA => "*",
	TAG => "*",
	TGT => "C",
	TGC => "C",
	TGA => "*",
	TGG => "W",
	CTT => "L",
	CTC => "L",
	CTA => "L",
	CTG => "L",
	CCT => "P",
	CCC => "P",
	CCA => "P",
	CCG => "P",
	CAT => "H",
	CAC => "H",
	CAA => "Q",
	CAG => "Q",
	CGT => "R",
	CGC => "R",
	CGA => "R",
	CGG => "R",
	ATT => "I",
	ATC => "I",
	ATA => "I",
	ATG => "M",
	ACT => "T",
	ACC => "T",
	ACA => "T",
	ACG => "T",
	AAT => "N",
	AAC => "N",
	AAA => "K",
	AAG => "K",
	AGT => "S",
	AGC => "S",
	AGA => "R",
	AGG => "R",
	GTT => "V",
	GTC => "V",
	GTA => "V",
	GTG => "V",
	GCT => "A",
	GCC => "A",
	GCA => "A",
	GCG => "A",
	GAT => "D",
	GAC => "D",
	GAA => "E",
	GAG => "E",
	GGT => "G",
	GGC => "G",
	GGA => "G",
	GGG => "G"
	);
	my %tt1_start = (TTG => "M",
	CTG => "M",
	ATG => "M"
	);
	my %tt1_stop = (TAA => "*",
	TAG => "*",
	TGA => "*"
	);

	#mito£¬transl_table=3
	my %transl_tab3 = (TTT => "F",
	TTC => "F",
	TTA => "L",
	TTG => "L",
	TCT => "S",
	TCC => "S",
	TCA => "S",
	TCG => "S",
	TAT => "Y",
	TAC => "Y",
	TAA => "*",
	TAG => "*",
	TGT => "C",
	TGC => "C",
	TGA => "W",
	TGG => "W",
	CTT => "T",
	CTC => "T",
	CTA => "T",
	CTG => "T",
	CCT => "P",
	CCC => "P",
	CCA => "P",
	CCG => "P",
	CAT => "H",
	CAC => "H",
	CAA => "Q",
	CAG => "Q",
	CGT => "R",
	CGC => "R",
	CGA => "R",
	CGG => "R",
	ATT => "I",
	ATC => "I",
	ATA => "M",
	ATG => "M",
	ACT => "T",
	ACC => "T",
	ACA => "T",
	ACG => "T",
	AAT => "N",
	AAC => "N",
	AAA => "K",
	AAG => "K",
	AGT => "S",
	AGC => "S",
	AGA => "R",
	AGG => "R",
	GTT => "V",
	GTC => "V",
	GTA => "V",
	GTG => "V",
	GCT => "A",
	GCC => "A",
	GCA => "A",
	GCG => "A",
	GAT => "D",
	GAC => "D",
	GAA => "E",
	GAG => "E",
	GGT => "G",
	GGC => "G",
	GGA => "G",
	GGG => "G"
	);
	my %tt3_start = (ATA => "M",
	ATG => "M",
	GTG => "M"
	);
	my %tt3_stop = (TAA => "*",
	TAG => "*"
	);
	
	
	my @ref_seqname;
	my @ref_seq;
	my %ref_name2idx;
	
	my @tgt_seqname;
	my @tgt_seq;
	my %tgt_name2idx;
	
	my %IDtype;		#{$ID} => type
	my %gene_info;	#{$gene_ID}->{ID/Name/contig/start/end/score/strand/phase/mRNA_num/mRNA_lst} {mRNA_lst}->[$idx] => mRNA_ID 
	my %mRNA_info;	#{$mRNA_ID}->{ID/Name/Parent/exon_lst/CDS_lst} {exon/CDS_lst}->[$idx] => exon/CDS ID
	my %exon_info;	#{$exon_ID}->{ID/Name/Parent/idx_in_exon_lst/contig/start/codon/score/strand/phase}
	my %CDS_info;	#{$CDS_ID}->{ID/Name/Parent/idx_in_CDS_lst/contig/start/codon/score/strand/phase}
	
	my @coor_map;	#[idx_of_record]
	my $map_num = 0;

	#summary of target model
	my @pred_mdl;	#predicted gene models, [idx]->{ref_gene_id,coor_map_rec_idx,gene_id,contig,start,end,strand,phase};
	my $pred_num = 0;
	#target models
	my %tgt_gene;
	my %tgt_mRNA;
	my %tgt_exon;
	my %tgt_CDS;
	my $tgt_gid;
	
	
	my @items;
	my @items2;
	my @items3;
	my $lbuf;
	my %tmph;
	my $lno = 0;
	
	#read ref and target fna
	ReadFasta("$ref_dir/$ref_fa", @ref_seqname, @ref_seq);
	foreach my $i (0 .. $#ref_seqname){
		$ref_name2idx{ $ref_seqname[$i] } = $i;
	}#
	
	ReadFasta($tgt_fna, @tgt_seqname, @tgt_seq);
	foreach my $i (0 .. $#tgt_seqname){
		$tgt_name2idx{ $tgt_seqname[$i] } = $i;
	}#
	#read models from gff
	open IN, "$ref_dir/$ref_gff" or die "fail to open $ref_dir/$ref_gff:$!\n";
	while (<IN>){
		$lno++;
		chomp;
		next if not $_;
		next if /^#/;
		@items = split /\t/,$_;
		@items2 = split /;/,$items[8];
		%tmph = ();
		foreach my $i (0 .. $#items2){
			@items3 = split /=/,$items2[$i];
			$tmph{ $items3[0] } = $items3[1];
		}#
		if (not exists $IDtype{ $tmph{ID} } ) {
			$IDtype{ $tmph{ID} } = $items[2];
		}else{
	#		die "duplicated ID $tmph{ID} at line $lno in $ref_dir/$ref_gff.\n";
		}#else
		if ($items[6] ne "+" and $items[6] ne "-"){
			die "strand not specified $items[6] at line $lno in $ref_dir/$ref_gff.\n";
		}#if
		if ($items[2] eq "gene"){
			if (exists $gene_info{ $tmph{ID} } ) {#
	#			die "duplicated gene ID at line $lno in $ref_dir/$ref_gff.\n";
			}#if
			$gene_info{ $tmph{ ID } }->{ID} = $tmph{ ID };
			$gene_info{ $tmph{ ID } }->{Name}=$tmph{ Name };
			$gene_info{ $tmph{ ID } }->{contig} =$items[0];
			$gene_info{ $tmph{ ID } }->{start}=$items[3];
			$gene_info{ $tmph{ ID } }->{end} = $items[4];
			$gene_info{ $tmph{ ID } }->{strand}= $items[6];
			$gene_info{ $tmph{ ID } }->{phase}= $items[7];
			$gene_info{ $tmph{ ID } }->{gene} = $tmph{ gene };
			$gene_info{ $tmph{ ID } }->{locus_tag} = $tmph{locus_tag};
			$gene_info{ $tmph{ ID } }->{mRNA_num} = 0;
		}elsif ($items[2] eq "mRNA"){
			if (exists $mRNA_info{ $tmph{ID} } ) {#
	#			die "duplicated mRNA ID at line $lno in $ref_dir/$ref_gff.\n";
			}#if
			$mRNA_info{ $tmph{ ID } }->{ID} = $tmph{ ID };
			$mRNA_info{ $tmph{ ID } }->{Name}=$tmph{ Name };
			$mRNA_info{ $tmph{ ID } }->{contig} =$items[0];
			$mRNA_info{ $tmph{ ID } }->{start}=$items[3];
			$mRNA_info{ $tmph{ ID } }->{end} = $items[4];
			$mRNA_info{ $tmph{ ID } }->{strand}= $items[6];
			$mRNA_info{ $tmph{ ID } }->{phase}= $items[7];
			$mRNA_info{ $tmph{ ID } }->{gene} = $tmph{ gene };
			$mRNA_info{ $tmph{ ID } }->{transcript_id} = $tmph{ transcript_id};
			$mRNA_info{ $tmph{ ID } }->{exon_num} = 0;
			$mRNA_info{ $tmph{ ID } }->{CDS_num} = 0;
			$mRNA_info{ $tmph{ ID } }->{locus_tag} = $tmph{ locus_tag };
			$mRNA_info{ $tmph{ ID } }->{product}   = $tmph{ product };
			
			$mRNA_info{ $tmph{ ID } }->{Parent} = $tmph{ Parent };
			if(exists $gene_info{ $tmph{ Parent } } ) {#because other features, e.g. pseudogene also contains mRNA
				push @{ $gene_info{ $tmph{ Parent } }->{mRNA_lst} }, $tmph{ ID };
				$gene_info{ $tmph{ Parent } }->{mRNA_num} ++;
			}else{
	#			print "NOPARENT_mRNA","\t","PARENTtype:",$IDtype{ $tmph{ Parent } },"\t",$lno,"\t",$_,"\n";
			}
		}elsif ($items[2] eq "exon"){
			if (exists $exon_info{ $tmph{ID} } ) {#
	#			die "duplicated exon ID at line $lno in $ref_dir/$ref_gff.\n";
			}#if		
			$exon_info{ $tmph{ ID } }->{ID} = $tmph{ ID };
			$exon_info{ $tmph{ ID } }->{Name}=$tmph{ Name };
			$exon_info{ $tmph{ ID } }->{contig} =$items[0];
			$exon_info{ $tmph{ ID } }->{start}=$items[3];
			$exon_info{ $tmph{ ID } }->{end} = $items[4];
			$exon_info{ $tmph{ ID } }->{strand}= $items[6];
			$exon_info{ $tmph{ ID } }->{phase}= $items[7];
			$exon_info{ $tmph{ ID } }->{gene} = $tmph{ gene };
			$exon_info{ $tmph{ ID } }->{locus_tag} = $tmph{ locus_tag };
			$exon_info{ $tmph{ ID } }->{product}   = $tmph{ product };
			$exon_info{ $tmph{ ID } }->{transcript_id} = $tmph{ transcript_id};
			$exon_info{ $tmph{ ID } }->{Parent} = $tmph{ Parent };
			if(exists $mRNA_info{ $tmph{ Parent } } ) {#other features, e.g. ncRNA, also contains exons
				push @{ $mRNA_info{ $tmph{ Parent } }->{exon_lst} }, $tmph{ ID };
				$mRNA_info{ $tmph{ Parent } }->{exon_num} ++;
			}else{
	#			print "NOPARENT_exon","\t","PARENTtype:",$IDtype{ $tmph{ Parent } },"\t",$lno,"\t",$_,"\n";
			}
	
		}elsif ($items[2] eq "CDS"){
			if (exists $CDS_info{ $tmph{ID} } ) {#
	#			die "duplicated CDS ID at line $lno in $ref_dir/$ref_gff.\n";
			}#if
			$CDS_info{ $tmph{ ID } }->{ID} = $tmph{ ID };
			$CDS_info{ $tmph{ ID } }->{Name}=$tmph{ Name };
			$CDS_info{ $tmph{ ID } }->{contig} =$items[0];
			$CDS_info{ $tmph{ ID } }->{start}=$items[3];
			$CDS_info{ $tmph{ ID } }->{end} = $items[4];
			$CDS_info{ $tmph{ ID } }->{strand}= $items[6];
			$CDS_info{ $tmph{ ID } }->{phase}= $items[7];
			$CDS_info{ $tmph{ ID } }->{gene} = $tmph{ gene };
			$CDS_info{ $tmph{ ID } }->{protein_id} = $tmph{ protein_id};
			$CDS_info{ $tmph{ ID } }->{locus_tag} = $tmph{ locus_tag };
			$CDS_info{ $tmph{ ID } }->{product}   = $tmph{ product };
			$CDS_info{ $tmph{ ID } }->{Note}      = $tmph{ Note };
			$CDS_info{ $tmph{ ID } }->{Parent} = $tmph{ Parent };
			if(exists $tmph{ transl_table } ) {
				$CDS_info{ $tmph{ ID } }->{transl_table} = $tmph{ transl_table };
				if ($tmph{transl_table} != 3){
					print "warning: unrecognized transl_table $tmph{ transl_table } at line $lno in $ref_dir/$ref_gff, please specify.\n";
				}#if
			}#if
			if(exists $mRNA_info{ $tmph{ Parent } } ) {#some pseudogenes also act as Parent of CDS (rather than as Parent of mRNA)
				push @{ $mRNA_info{ $tmph{ Parent } }->{CDS_lst} }, $tmph{ ID };
				$mRNA_info{ $tmph{ Parent } }->{CDS_num} ++;
			}elsif(exists $gene_info{ $tmph{ Parent } } ) {#gene -> CDS; no mRNA record
	#			print "gene-CDS","\t",$lno,"\n";
				push @{ $gene_info{ $tmph{ Parent } }->{CDS_lst} }, $tmph{ ID };
				$gene_info{ $tmph{ Parent } }->{CDS_num} ++;
			}else{
	#			print "NOPARENT_CDS","\t",$lno,"\t","PARENTtype:",$IDtype{ $tmph{ Parent } },"\t",$_,"\n";
			}
	
		}else{
			#do nothing about other records
		}#else
	}#while
	close IN;
	
	
		
	#read map file
	open IN, "$mapfile" or die "fail to open $mapfile:$!\n";
	while (<IN>){
		chomp;
		next if not $_;
		next if /^#/;
		@{ $coor_map[ $map_num ++ ] } = split /\t/,$_;
	}
	close IN;
	
	#check ref gff models
	
	my %mrna_iso;	#{num of iso} => cnt of genes
	my $contain_CDS = 0;
	my $contain_CDS2 = 0;
	print "#gene models contain CDS but no mRNA\n";
	print "#gene\tmRNA_num\tCDS_num\n";
	foreach my $i (sort keys %gene_info){
		if (not $gene_info{$i}->{mRNA_num} ){
			$mrna_iso{ 0 } ++;
			if (exists $gene_info{$i}->{CDS_num} ) {
				$contain_CDS ++;
			}#if
		}else{
			$mrna_iso{ $gene_info{$i}->{mRNA_num} } ++;
		}#
		if ($gene_info{$i}->{CDS_num} ) {
			$contain_CDS2 ++;
			print $i,"\t",$gene_info{$i}->{mRNA_num},"\t",$gene_info{$i}->{CDS_num},"\n";
		}#if
		
	}
	print "contain_CDS","\t",$contain_CDS,"\n";		
	print "contain_CDS2","\t",$contain_CDS2,"\n";		

	print "#mrna_iso_num\tmrna_num\n";
	foreach my $i (sort {$a<=>$b} keys %mrna_iso){
		print $i,"\t",$mrna_iso{$i},"\n";
	}#
	
	#check start codon and stop codon
	my $min;
	my $max;
	my $startcodon;
	my $stopcodon;
	my $tmpseq;
	
	#check seq and get start/stop codon here
	foreach my $i (sort keys %mRNA_info){#
		$tmpseq = "";
		if (exists $mRNA_info{$i}->{CDS_lst}){
			foreach my $j (sort {$CDS_info{ $mRNA_info{$i}->{CDS_lst}->[$a] }->{start} <=> $CDS_info{ $mRNA_info{$i}->{CDS_lst}->[$b] }->{start}
				or $CDS_info{ $mRNA_info{$i}->{CDS_lst}->[$a] }->{end} <=> $CDS_info{ $mRNA_info{$i}->{CDS_lst}->[$b] }->{end} } (0 .. $#{$mRNA_info{$i}->{CDS_lst} } ) ) {#
				$tmpseq .= substr $ref_seq[ $ref_name2idx{ $mRNA_info{ $i }->{contig} } ], $CDS_info{$mRNA_info{$i}->{CDS_lst}->[$j]}->{start} - 1, $CDS_info{$mRNA_info{$i}->{CDS_lst}->[$j]}->{end} - $CDS_info{$mRNA_info{$i}->{CDS_lst}->[$j]}->{start} + 1;
			}
			
			if ( $mRNA_info{$i}->{strand} eq "-"){
				$tmpseq = rc($tmpseq);
			}#if
	#		print ">",$i," ", $mRNA_info{$i}->{strand},"\n",$tmpseq,"\n";
			$startcodon = substr $tmpseq, 0, 3;
			$stopcodon  = substr $tmpseq, (length $tmpseq) - 3, 3;
	
			if ($mRNA_info{ $i }->{contig} eq $mit_seq){#tt3
				if(exists $tt3_start{ uc( $startcodon ) }){
					#ok, do nothing
				}else{
					print "gene $i wrong start codon $startcodon.\n";
				}#else
				if(exists $tt3_stop{ uc( $stopcodon ) } ) {
					#ok, do nothing
				}else{
					print "gene $i wrong stop codon $stopcodon.\n";
				}#else
			}else{#tt1
				if(exists $tt1_start{ uc( $startcodon ) }){
					#ok, do nothing
				}else{
					print "gene $i wrong start codon $startcodon.\n";
				}#else
				if(exists $tt1_stop{ uc( $stopcodon ) } ) {
					#ok, do nothing
				}else{
					print "gene $i wrong stop codon $stopcodon.\n";
				}#else
			}#else tt1
	
		}#if
	}
		
	
	foreach my $i (sort keys %gene_info){#
		$tmpseq = "";
		if (exists $gene_info{$i}->{CDS_lst}){
			foreach my $j (sort {$CDS_info{ $gene_info{$i}->{CDS_lst}->[$a] }->{start} <=> $CDS_info{ $gene_info{$i}->{CDS_lst}->[$b] }->{start}
				or $CDS_info{ $gene_info{$i}->{CDS_lst}->[$a] }->{end} <=> $CDS_info{ $gene_info{$i}->{CDS_lst}->[$b] }->{end} } (0 .. $#{$gene_info{$i}->{CDS_lst} } ) ) {#
				$tmpseq .= substr $ref_seq[ $ref_name2idx{ $gene_info{ $i }->{contig} } ], $CDS_info{$gene_info{$i}->{CDS_lst}->[$j]}->{start} - 1, $CDS_info{$gene_info{$i}->{CDS_lst}->[$j]}->{end} - $CDS_info{$gene_info{$i}->{CDS_lst}->[$j]}->{start} + 1;
			}
			
			if ( $gene_info{$i}->{strand} eq "-"){
				$tmpseq = rc($tmpseq);
			}#if
	#		print ">",$i," ", $gene_info{$i}->{strand},"\n",$tmpseq,"\n";
			$startcodon = substr $tmpseq, 0, 3;
			$stopcodon  = substr $tmpseq, (length $tmpseq) - 3, 3;
	
			if ($gene_info{ $i }->{contig} eq $mit_seq){#tt3
				if(exists $tt3_start{ uc( $startcodon ) }){
					#ok, do nothing
				}else{
					print "gene $i wrong start codon $startcodon.\n";
				}#else
				if(exists $tt3_stop{ uc( $stopcodon ) } ) {
					#ok, do nothing
				}else{
					print "gene $i wrong stop codon $stopcodon.\n";
				}#else
			}else{#tt1
				if(exists $tt1_start{ uc( $startcodon ) }){
					#ok, do nothing
				}else{
					print "gene $i wrong start codon $startcodon.\n";
				}#else
				if(exists $tt1_stop{ uc( $stopcodon ) } ) {
					#ok, do nothing
				}else{
					print "gene $i wrong stop codon $stopcodon.\n";
				}#else
			}#else tt1
	
		}#if
	}
		
	
	
	#map exon/CDS/mRNA/gene
	#first check whether this gff file contains UTR
	my $s2;
	my $e2;
	foreach my $i (keys %mRNA_info){
		if ($mRNA_info{$i}->{exon_num}){
			$s2 = $exon_info{ $mRNA_info{$i}->{exon_lst}->[0] }->{start};
			$e2 = $exon_info{ $mRNA_info{$i}->{exon_lst}->[0] }->{end};
			foreach my $j (1 .. $#{ $mRNA_info{$i}->{exon_lst} } ){
				if ($exon_info{ $mRNA_info{$i}->{exon_lst}->[$j] }->{start} < $s2){
					$s2 = $exon_info{ $mRNA_info{$i}->{exon_lst}->[$j] }->{start};
				}#if
				if ($exon_info{ $mRNA_info{$i}->{exon_lst}->[$j] }->{end} > $e2){
					$e2 = $exon_info{ $mRNA_info{$i}->{exon_lst}->[$j] }->{end};
				}#if
			}#
			if ($s2 != $mRNA_info{$i}->{start} or $e2 != $mRNA_info{$i}->{end} ){#
	#			print "mRNA start/end and exon start/end do not match for mRNA ID $i\n";
			}#if
		}#if
	}
				
	
	foreach my $i (keys %mRNA_info){
		if ($mRNA_info{$i}->{CDS_num}){
			$s2 = $CDS_info{ $mRNA_info{$i}->{CDS_lst}->[0] }->{start};
			$e2 = $CDS_info{ $mRNA_info{$i}->{CDS_lst}->[0] }->{end};
			foreach my $j (1 .. $#{ $mRNA_info{$i}->{CDS_lst} } ){
				if ($CDS_info{ $mRNA_info{$i}->{CDS_lst}->[$j] }->{start} < $s2){
					$s2 = $CDS_info{ $mRNA_info{$i}->{CDS_lst}->[$j] }->{start};
				}#if
				if ($CDS_info{ $mRNA_info{$i}->{CDS_lst}->[$j] }->{end} > $e2){
					$e2 = $CDS_info{ $mRNA_info{$i}->{CDS_lst}->[$j] }->{end};
				}#if
			}#
			if ($s2 != $mRNA_info{$i}->{start} or $e2 != $mRNA_info{$i}->{end} ){#
	#			print "mRNA start/end and CDS start/end do not match for mRNA ID $i\n";
			}#if
		}#if
	}
	
			
	#the gene/mRNA
	foreach my $i (keys %gene_info){
		if ($gene_info{$i}->{mRNA_num}){
			$s2 = $mRNA_info{ $gene_info{$i}->{mRNA_lst}->[0] }->{start};
			$e2 = $mRNA_info{ $gene_info{$i}->{mRNA_lst}->[0] }->{end};
			foreach my $j (1 .. $#{ $gene_info{$i}->{mRNA_lst} } ){
				if ($mRNA_info{ $gene_info{$i}->{mRNA_lst}->[$j] }->{start} < $s2){
					$s2 = $mRNA_info{ $gene_info{$i}->{mRNA_lst}->[$j] }->{start};
				}#if
				if ($mRNA_info{ $gene_info{$i}->{mRNA_lst}->[$j] }->{end} > $e2){
					$e2 = $mRNA_info{ $gene_info{$i}->{mRNA_lst}->[$j] }->{end};
				}#if
			}#
			if ($s2 != $gene_info{$i}->{start} or $e2 != $gene_info{$i}->{end} ){#
	#			print "gene start/end and mRNA start/end do not match for gene ID $i\n";
			}#if
		}#if
	}
	
	#based on above check, start and stop are the same for gene/mRNA/exon/CDS except that some mRNA may have introns
	#therefore, only gene boundary is checked.
	#there may be two types reliable predictions
	#1. the full gene fall with in one line of map2 file
	#2. there are short indels in genes, allow the full gene span adjacent lines, and reading frames will be checked for such models for find possible premature/ frameshift.
	
	#to speed up, create index for each ref nt seq
	my %ref_map_paragraph;	#{ $ref }->[$idx_paragraph]->{ref_str/tgtname/idx_first_rec/idx_last_rec}
	my %ref_map_paragraph_num; #{ $ref } => num of paragraph
	my $tgtname = $coor_map[0]->[1];
	my $refname = $coor_map[0]->[2];
	my $refstr  = $coor_map[0]->[3];
	my $firstrec = 1;
	foreach my $i (1 .. $#coor_map){
		if ($coor_map[$i]->[0] eq "SQ"){#
			#save last paragraph
			if (not exists $ref_map_paragraph{ $refname } ) {#
				$ref_map_paragraph_num{$refname} = 0;
			}#
			$ref_map_paragraph{ $refname }->[ $ref_map_paragraph_num{$refname} ]->{ref_str} = $refstr;
			$ref_map_paragraph{ $refname }->[ $ref_map_paragraph_num{$refname} ]->{tgtname}  = $tgtname;
			$ref_map_paragraph{ $refname }->[ $ref_map_paragraph_num{$refname} ]->{idx_first_rec} = $firstrec;
			$ref_map_paragraph{ $refname }->[ $ref_map_paragraph_num{$refname} ]->{idx_last_rec} = $i - 1;
			$ref_map_paragraph_num{$refname} ++;
			
			$tgtname = $coor_map[$i]->[1];
			$refname = $coor_map[$i]->[2];
			$refstr = $coor_map[$i]->[3];
			$firstrec = $i +1;
		}#if
	}#
	#save last paragraph
	if (not exists $ref_map_paragraph{ $refname } ) {#
		$ref_map_paragraph_num{$refname} = 0;
	}#
	$ref_map_paragraph{ $refname }->[ $ref_map_paragraph_num{$refname} ]->{ref_str} = $refstr;
	$ref_map_paragraph{ $refname }->[ $ref_map_paragraph_num{$refname} ]->{tgtname}  = $tgtname;
	$ref_map_paragraph{ $refname }->[ $ref_map_paragraph_num{$refname} ]->{idx_first_rec} = $firstrec;
	$ref_map_paragraph{ $refname }->[ $ref_map_paragraph_num{$refname} ]->{idx_last_rec} = $#coor_map;
	$ref_map_paragraph_num{$refname} ++;
	
	#output to check
	print "#ref_contig\tr_para_num\tr_para_num2\tidx\tstrand\ttgt_name\tnum_rec\tfirst_rec\tlast_rec\n";
	foreach my $i (sort keys %ref_map_paragraph){
		print $i,"\t",$ref_map_paragraph_num{$i},"\t",scalar @{$ref_map_paragraph{$i}},"\n";
		foreach my $j ( 0 .. $#{ $ref_map_paragraph{$i} } ) {
			print "\t\t\t",$j,"\t",$ref_map_paragraph{$i}->[$j]->{ref_str},
			"\t",$ref_map_paragraph{$i}->[$j]->{tgtname},
			"\t",$ref_map_paragraph{$i}->[$j]->{idx_last_rec} - $ref_map_paragraph{$i}->[$j]->{idx_first_rec} + 1,
			"\t",$ref_map_paragraph{$i}->[$j]->{idx_first_rec},
			"\t",$ref_map_paragraph{$i}->[$j]->{idx_last_rec},"\n";
		}#
	}#
	
	
	
	open OUT, ">$outfile" or die "fail to open $outfile:$!\n";
	#produce homo models
	my $cnt;
	foreach my $i (sort keys %gene_info){
		if (not $gene_info{$i}->{mRNA_num} and not $gene_info{$i}->{CDS_num}) {
			next;
		}#
		#check each record in the coor_map
		#check the appropriate records based on the %ref_map_paragraph
		if (exists $ref_map_paragraph{ $gene_info{$i}->{contig} } ){#
			foreach my $j (0 .. $#{ $ref_map_paragraph{ $gene_info{$i}->{contig} } } ) {#each paragraph
				foreach my $k ($ref_map_paragraph{ $gene_info{$i}->{contig} }->[$j]->{idx_first_rec} .. $ref_map_paragraph{ $gene_info{$i}->{contig} }->[$j]->{idx_last_rec} ) {#each rec
					if ( $gene_info{$i}->{start} >= $coor_map[$k]->[2]
					and  $gene_info{$i}->{start} <= $coor_map[$k]->[3]
					and  $gene_info{$i}->{end}   >= $coor_map[$k]->[2]
					and  $gene_info{$i}->{end}   <= $coor_map[$k]->[3] ) {#yes, type I
						$pred_mdl[ $pred_num ]->{ref_gene_id} = $i;
						$pred_mdl[ $pred_num ]->{ref_locus_tag} = $gene_info{$i}->{locus_tag};
						if (exists $gene_info{$i}->{mRNA_lst}){
							$pred_mdl[ $pred_num ]->{ref_transcript_id} = $mRNA_info{ $gene_info{$i}->{mRNA_lst}->[0] }->{transcript_id};
						}else{
	#						$pred_mdl[ $pred_num ]->{ref_transcript_id} = "?";
						}#else
						
						$pred_mdl[ $pred_num ]->{coor_map_rec_idx} = $k;
						$pred_mdl[ $pred_num ]->{contig}      = $ref_map_paragraph{ $gene_info{$i}->{contig} }->[$j]->{tgtname};
						#strand
		
						if ($gene_info{ $i }->{strand} eq $ref_map_paragraph{ $gene_info{$i}->{contig} }->[$j]->{ref_str} ) {
							$pred_mdl[ $pred_num ]->{strand} = "+";
						}else{#
							$pred_mdl[ $pred_num ]->{strand} = "-";
						}#else
						
						#to speed up the calculation
						#cal the pos directly by calculating the offset between the target contig and ref contig
						if ($ref_map_paragraph{ $gene_info{$i}->{contig} }->[$j]->{ref_str} eq "+"){#
							$pred_mdl[ $pred_num ]->{start } = $gene_info{$i}->{start} + ($coor_map[ $k ]->[ 0 ] - $coor_map[ $k ]->[ 2 ]);
							$pred_mdl[ $pred_num ]->{end }   = $gene_info{$i}->{end}   + ($coor_map[ $k ]->[ 0 ] - $coor_map[ $k ]->[ 2 ]);
							#get CDS and exon, if exists
							if (exists $gene_info{$i}->{CDS_lst} ) {
								$pred_mdl[ $pred_num ]->{CDS_num} = 0;
								$cnt = 0;
								foreach my $m (sort { $CDS_info{$a}->{start} <=> $CDS_info{$b}->{start}
								or $CDS_info{$a}->{end} <=> $CDS_info{$b}->{end} } @{ $gene_info{$i}->{CDS_lst} } ) {#each CDS
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{start} = $CDS_info{$m}->{start}  + ($coor_map[ $k ]->[ 0 ] - $coor_map[ $k ]->[ 2 ]);
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{end}   = $CDS_info{$m}->{end}    + ($coor_map[ $k ]->[ 0 ] - $coor_map[ $k ]->[ 2 ]);
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{phase} = $CDS_info{$m}->{phase};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{Note}  = $CDS_info{$m}->{Note};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{product}=$CDS_info{$m}->{product};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{ref_protein_id} =$CDS_info{$m}->{protein_id};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{ref_locus_tag}  =$CDS_info{$m}->{locus_tag};
									if(exists $CDS_info{$m}->{transl_table}){
										$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{transl_table} = $CDS_info{$m}->{transl_table};
									}#
									$cnt ++;
								}#
								$pred_mdl[ $pred_num ]->{CDS_num} = $cnt;
							}elsif (exists $gene_info{$i}->{mRNA_lst} ){
								#current ref models do not contains alternative splicing models; one-gene-one-mRNA
								$pred_mdl[ $pred_num ]->{CDS_num} = 0;
								$cnt = 0;
								foreach my $m (sort { $CDS_info{$a}->{start} <=> $CDS_info{$b}->{start}
								or $CDS_info{$a}->{end} <=> $CDS_info{$b}->{end} } @{ $mRNA_info{ $gene_info{$i}->{mRNA_lst}->[0] }->{CDS_lst} } ) {#each CDS
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{start} = $CDS_info{$m}->{start}  + ($coor_map[ $k ]->[ 0 ] - $coor_map[ $k ]->[ 2 ]);
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{end}   = $CDS_info{$m}->{end}    + ($coor_map[ $k ]->[ 0 ] - $coor_map[ $k ]->[ 2 ]);
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{phase} = $CDS_info{$m}->{phase};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{Note}  = $CDS_info{$m}->{Note};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{product}=$CDS_info{$m}->{product};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{ref_protein_id} =$CDS_info{$m}->{protein_id};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{ref_locus_tag}  =$CDS_info{$m}->{locus_tag};
									if(exists $CDS_info{$m}->{transl_table}){
										$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{transl_table} = $CDS_info{$m}->{transl_table};
									}#
									$cnt ++;
								}#
								$pred_mdl[ $pred_num ]->{CDS_num} = $cnt;
								#produce exon rec
								$pred_mdl[ $pred_num ]->{exon_num} = 0;
								$cnt = 0;
								foreach my $m (sort { $exon_info{$a}->{start} <=> $exon_info{$b}->{start}
								or $exon_info{$a}->{end} <=> $exon_info{$b}->{end} } @{ $mRNA_info{ $gene_info{$i}->{mRNA_lst}->[0] }->{exon_lst} } ) {#each exon
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{start} = $exon_info{$m}->{start} + ($coor_map[ $k ]->[ 0 ] - $coor_map[ $k ]->[ 2 ]);
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{end}   = $exon_info{$m}->{end}   + ($coor_map[ $k ]->[ 0 ] - $coor_map[ $k ]->[ 2 ]);
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{phase} = $exon_info{$m}->{phase};
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{product}=$exon_info{$m}->{product};
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{ref_transcript_id} = $exon_info{$m}->{transcript_id};
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{ref_locus_tag} = $exon_info{$m}->{locus_tag};
									$cnt ++;
								}#
								$pred_mdl[ $pred_num ]->{exon_num} = $cnt;
								
							}else{
								#do nothing, do not produce detailed models
							}#else
						}else{#-
							$pred_mdl[ $pred_num ]->{start } = $coor_map[ $k ]->[ 2 ] - $gene_info{$i}->{end}   + $coor_map[ $k ]->[ 0 ];
							$pred_mdl[ $pred_num ]->{end }   = $coor_map[ $k ]->[ 1 ] - $gene_info{$i}->{start} + $coor_map[ $k ]->[ 3 ];
							#get CDS and exon, if exists
							if (exists $gene_info{$i}->{CDS_lst} ) {
								$pred_mdl[ $pred_num ]->{CDS_num} = 0;
								$cnt = 0;
								foreach my $m (sort { $CDS_info{$a}->{start} <=> $CDS_info{$b}->{start}
								or $CDS_info{$a}->{end} <=> $CDS_info{$b}->{end} } @{ $gene_info{$i}->{CDS_lst} } ) {#each CDS
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{start} = $coor_map[ $k ]->[ 2 ] - $CDS_info{$m}->{end}   + $coor_map[ $k ]->[ 0 ];
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{end}   = $coor_map[ $k ]->[ 1 ] - $CDS_info{$m}->{start} + $coor_map[ $k ]->[ 3 ];
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{phase} = $CDS_info{$m}->{phase};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{Note}  = $CDS_info{$m}->{Note};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{product}=$CDS_info{$m}->{product};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{ref_protein_id} =$CDS_info{$m}->{protein_id};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{ref_locus_tag}  =$CDS_info{$m}->{locus_tag};
									if(exists $CDS_info{$m}->{transl_table}){
										$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{transl_table} = $CDS_info{$m}->{transl_table};
									}#
									$cnt ++;
								}#
								$pred_mdl[ $pred_num ]->{CDS_num} = $cnt;
							}elsif (exists $gene_info{$i}->{mRNA_lst}){
								#current ref models do not contains alternative splicing models; one-gene-one-mRNA
								$pred_mdl[ $pred_num ]->{CDS_num} = 0;
								$cnt = 0;
								foreach my $m (sort { $CDS_info{$a}->{start} <=> $CDS_info{$b}->{start}
								or $CDS_info{$a}->{end} <=> $CDS_info{$b}->{end} } @{ $mRNA_info{ $gene_info{$i}->{mRNA_lst}->[0] }->{CDS_lst} } ) {#each CDS
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{start} = $coor_map[ $k ]->[ 2 ] - $CDS_info{$m}->{end}   + $coor_map[ $k ]->[ 0 ];
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{end}   = $coor_map[ $k ]->[ 1 ] - $CDS_info{$m}->{start} + $coor_map[ $k ]->[ 3 ];
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{phase} = $CDS_info{$m}->{phase};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{Note}  = $CDS_info{$m}->{Note};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{product}=$CDS_info{$m}->{product};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{ref_protein_id} =$CDS_info{$m}->{protein_id};
									$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{ref_locus_tag}  =$CDS_info{$m}->{locus_tag};
									if(exists $CDS_info{$m}->{transl_table}){
										$pred_mdl[$pred_num]->{CDS_lst}->[$cnt]->{transl_table} = $CDS_info{$m}->{transl_table};
									}#
									$cnt ++;
								}#
								$pred_mdl[ $pred_num ]->{CDS_num} = $cnt;
								#get exon rec
								$pred_mdl[ $pred_num ]->{exon_num} = 0;
								$cnt = 0;
								foreach my $m (sort { $exon_info{$a}->{start} <=> $exon_info{$b}->{start}
								or $exon_info{$a}->{end} <=> $exon_info{$b}->{end} } @{ $mRNA_info{ $gene_info{$i}->{mRNA_lst}->[0] }->{exon_lst} } ) {#each exon
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{start} = $coor_map[ $k ]->[ 2 ] - $exon_info{$m}->{end}   + $coor_map[ $k ]->[ 0 ];
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{end}   = $coor_map[ $k ]->[ 1 ] - $exon_info{$m}->{start} + $coor_map[ $k ]->[ 3 ];
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{phase} = $exon_info{$m}->{phase};
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{product}=$exon_info{$m}->{product};
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{ref_transcript_id} = $exon_info{$m}->{transcript_id};
									$pred_mdl[$pred_num]->{exon_lst}->[$cnt]->{ref_locus_tag}     = $exon_info{$m}->{locus_tag};
									$cnt ++;
								}#
								$pred_mdl[ $pred_num ]->{exon_num} = $cnt;
							}else{
								#do nothing, do not produce detailed models
							}#else
	
						}#else -
						$pred_mdl[ $pred_num ]->{phase}      = $gene_info{$i}->{phase};
						$pred_num ++;
					}#if
				}#foreach my $k
			}#foreach my $j
		}#if
	}#
	
	
	#produce final models
	$tgt_gid = 0;#start from 1
	my $tmpid;
	my $codon;
	my $tmpaa;
	my $tmpgene;
	my $tmpcds;
	my $tmpexon;
	my $tmpaa2;
	my $tmpgene2;
	my $tmpcds2;
	my $tmpexon2;
	my $mdl_stat = 0;#
	my $extaa;
	my $extgene;
	my $extcds;
	my $extexon;
	my $new_start;
	my $new_end;
	my $start_bound;
	my $end_bound;

#	open GFF, "|pigz >$outgff.gz" or die "fail to open $outgff:$!\n";
#	open GENE, "|pigz >$outgene.gz" or die "fail to open $outgene:$!\n";
#	open MRNA, "|pigz >$outmrna.gz" or die "fail to open $outmrna:$!\n";
#	open PROT, "|pigz >$outprot.gz" or die "fail to open $outprot:$!\n";
#	open CDS,  "|pigz >$outcds.gz" or die "fail to open $outcds:$!\n";

	open GFF, ">$outgff" or die "fail to open $outgff:$!\n";
	open GENE, ">$outgene" or die "fail to open $outgene:$!\n";
	open MRNA, ">$outmrna" or die "fail to open $outmrna:$!\n";
	open PROT, ">$outprot" or die "fail to open $outprot:$!\n";
	open CDS,  ">$outcds" or die "fail to open $outcds:$!\n";

	open PGFF, "|pigz >$pseudogff.gz" or die "fail to open $pseudogff.gz:$!\n";
	open PGENE, "|pigz >$pseudogene.gz" or die "fail to open $pseudogene.gz.gz:$!\n";
	open PMRNA, "|pigz >$pseudomrna.gz" or die "fail to open $pseudomrna.gz:$!\n";
	open PPROT, "|pigz >$pseudoprot.gz" or die "fail to open $pseudoprot.gz:$!\n";
	open PCDS,  "|pigz >$pseudocds.gz" or die "fail to open $pseudocds.gz:$!\n";

	foreach my $i (sort { $pred_mdl[$a]->{contig} cmp $pred_mdl[$b]->{contig}
		or $pred_mdl[$a]->{start} <=> $pred_mdl[$b]->{start} 
		or $pred_mdl[$a]->{end}   <=> $pred_mdl[$b]->{end} } (0 .. $#pred_mdl) ){
	
		#print OUT "DEBUG:","\t",$i,"\t",$pred_mdl[$i]->{ref_transcript_id},"\n";
		#get CDS seq
		$tmpcds = "";#for CDS, for some genes, CDS is different from exon
		$tmpexon = "";#for (exon of) mRNA (exon of mRNA and cds are almost identical here, UTRs are ignored)
		$tmpaa  = "";
		$tmpgene = "";
		$extcds = "";
		
		if (exists $pred_mdl[$i]->{CDS_lst}){
			if ($pred_mdl[$i]->{strand} eq "+"){
				if ($pred_mdl[$i]->{CDS_lst}->[0]->{phase} != 0 ) {
					die "phase of the first codon is NOT 0!!!\n";
				}#
			}else{#-
				if ($pred_mdl[$i]->{CDS_lst}->[-1]->{phase} != 0 ) {
					die "phase of the first codon is NOT 0!!!\n";
				}#if
			}#else
			#get gene		
			$tmpgene = substr $tgt_seq[ $tgt_name2idx{ $pred_mdl[$i]->{contig} } ], $pred_mdl[$i]->{start} - 1, $pred_mdl[$i]->{end} - $pred_mdl[$i]->{start} + 1;
			if ($pred_mdl[$i]->{strand} eq "-"){
				$tmpgene = rc ($tmpgene);
			}#
			$tmpgene = uc ($tmpgene);
			#get exon/mRNA
			if (exists $pred_mdl[$i]->{exon_lst} ) {
				foreach my $j (0 .. $#{$pred_mdl[$i]->{exon_lst}}){
					$tmpexon .= substr $tgt_seq[ $tgt_name2idx{ $pred_mdl[$i]->{contig} } ], $pred_mdl[$i]->{exon_lst}->[$j]->{start} - 1, $pred_mdl[$i]->{exon_lst}->[$j]->{end} - $pred_mdl[$i]->{exon_lst}->[$j]->{start} + 1;
				}
				if ($pred_mdl[$i]->{strand} eq "-"){
					$tmpexon = rc($tmpexon);
				}#
				$tmpexon = uc($tmpexon);
			}#if
			#get cds		
			foreach my $j (0 .. $#{$pred_mdl[$i]->{CDS_lst}}){
				$tmpcds .= substr $tgt_seq[ $tgt_name2idx{ $pred_mdl[$i]->{contig} } ], $pred_mdl[$i]->{CDS_lst}->[$j]->{start} - 1, $pred_mdl[$i]->{CDS_lst}->[$j]->{end} - $pred_mdl[$i]->{CDS_lst}->[$j]->{start} + 1;
			}
			if ($pred_mdl[$i]->{strand} eq "-"){
				$tmpcds = rc($tmpcds);
			}#if
			$tmpcds = uc($tmpcds);
			
			#check each codon
			#check start codon
			$startcodon = substr $tmpcds,0,3;
			$stopcodon  = substr $tmpcds, (length $tmpcds) - 3, 3;
			if ($mit_seq eq $gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{contig} ) {#ref mit seq, use tt3
				#start codon
				if (exists $tt3_start{ $startcodon } ) {#start codon ok
					$tmpaa = $tt3_start{ $startcodon };
				}elsif ($startcodon =~ /[^ATGC]/){
					$tmpaa = "X";
					print "warning: start codon contains N.\t";
					print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
					$pred_mdl[$i]->{startcodon_containN} ++;
				}elsif (exists $tt3_stop{ $startcodon } ) {#start codon to stop codon, 20210918
					$tmpaa = "*";
					$pred_mdl[$i]->{premature} ++;
					print "warning: start codon mutated to stop codon, pseudogene\t";
					print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
				}else{
					$tmpaa = $transl_tab3{ $startcodon };
					print "warning: start codon mutated to non-start codon.\t";
					print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
					$pred_mdl[$i]->{startcodon_mutation} ++;
				}#else
				#check other codons, and produce aa seq
				$cnt = (length $tmpcds) / 3;#num of codon
				
				foreach my $j ( 1 .. $cnt - 2){
					$codon = substr $tmpcds, $j*3, 3;
					if (exists $tt3_stop{ $codon } ) {#premature, non-sense mutation
						$pred_mdl[$i]->{premature} ++;
						$tmpaa .= "*";
						print "warning: stop codon in the middle of gene, pseudogene\t";
						print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
					}elsif ($codon =~ /[^ATGC]/){
						$tmpaa .= "X";
						print "warning: N in the middle of gene\t";
						print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
					}else{
						$tmpaa .= $transl_tab3{ $codon };
					}#else
				}#
				#stop codon
				if (exists $tt3_stop{ $stopcodon } ) {#stop codon ok
	#				$tmpaa .= "*";
				}elsif( $stopcodon =~ /[^ATGC]/){
					$pred_mdl[$i]->{stopcodon_containN} ++;
					print "warning: stop codon contains N.\t";
					print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
				}else{
					$pred_mdl[$i]->{stopcodon_mutation} ++;
					print "warning: stop codon mutated to non-stop codon.\t";
					print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
				}#else
				$pred_mdl[$i]->{nt_seq} = $tmpcds;
				$pred_mdl[$i]->{aa_seq} = $tmpaa;
			}else{#transl tab 1
				#start codon
				if (exists $tt1_start{ $startcodon } ) {#start codon ok
					$tmpaa = $tt1_start{ $startcodon };
				}elsif ($startcodon =~ /[^ATGC]/){
					$tmpaa = "X";
					print "warning: start codon contains N.\t";
					print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
					$pred_mdl[$i]->{startcodon_containN} ++;
				}elsif (exists $tt1_stop{ $startcodon } ) {#start codon to stop codon, 20210923
					$tmpaa = "*";
					$pred_mdl[$i]->{premature} ++;
					print "warning: start codon mutated to stop codon, pseudogene\t";
					print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
				}else{
					$tmpaa = $transl_tab1{ $startcodon };
					print "warning: start codon mutated to non-start codon.\t";
					print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
					$pred_mdl[$i]->{startcodon_mutation} ++;
				}#else
				#check other codons, and produce aa seq
				$cnt = (length $tmpcds) / 3;#num of codon
				
				foreach my $j ( 1 .. $cnt - 2){
					$codon = substr $tmpcds, $j*3, 3;
					if (exists $tt1_stop{ $codon } ) {#premature, non-sense mutation
						$pred_mdl[$i]->{premature} ++;
						$tmpaa .= "*";
						print "warning: stop codon in the middle of gene, pseudogene\t";
						print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
					}elsif ($codon =~ /[^ATGC]/){
						$tmpaa .= "X";
						print "warning: N in the middle of gene\t";
						print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
					}else{
						$tmpaa .= $transl_tab1{ $codon };
					}#else
				}#
				#stop codon
				if (exists $tt1_stop{ $stopcodon } ) {#stop codon ok
	#				$tmpaa .= "*";
				}elsif( $stopcodon =~ /[^ATGC]/) {
					$pred_mdl[$i]->{stopcodon_containN} ++;
					print "warning: stop codon contains N.\t";
					print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
				}else{
					$pred_mdl[$i]->{stopcodon_mutation} ++;
					print "warning: stop codon mutated to non-stop codon.\t";
					print "map_coor:\t",$pred_mdl[$i]->{contig},"\t",$pred_mdl[$i]->{start},"\t",$pred_mdl[$i]->{end},"\t",$pred_mdl[$i]->{strand},"\t","map_ref:\t",$pred_mdl[$i]->{ref_gene_id},"\n";
				}#else
				$pred_mdl[$i]->{nt_seq} = $tmpcds;
				$pred_mdl[$i]->{aa_seq} = $tmpaa;
			}#else transl tab 1

			#check above model stat flags based on start/stop codon mutations
			$mdl_stat = 0;
			$extaa = "";
			$extgene = "";
			$extcds = "";
			$extexon= "";
			if(exists $pred_mdl[$i]->{premature}){#pseudogene
				$mdl_stat = 1;
			}elsif (exists $pred_mdl[$i]->{stopcodon_mutation}){#
				#try to extend, extend to the first possible stop codon or to the end of the contig
				if ($pred_mdl[$i]->{strand} eq "+"){#forward
					$flag = 0;
					$end_bound = length $tgt_seq[ $tgt_name2idx{ $pred_mdl[$i]->{contig} } ];
					$cnt = (($end_bound - $pred_mdl[$i]->{end}) - ($end_bound - $pred_mdl[$i]->{end})%3) / 3;#max num of possible condons to extend available at the contig
					if ($cnt > $max_ext){
						$cnt = $max_ext;
					}#if
					foreach my $j (1 .. $cnt){#how to process N???
						$codon     = uc ( substr $tgt_seq[ $tgt_name2idx{ $pred_mdl[$i]->{contig} } ], $pred_mdl[$i]->{end} + ($j - 1) * 3, 3 );#20210917
						$extcds   .= $codon;
						$extgene  .= $codon;
						$extexon  .= $codon;
						
						if ( $mit_seq eq $gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{contig} ) {#ref mit seq, use tt3
							if( exists $tt3_stop{ $codon } ) {#this is a stop codon
								$new_end = $pred_mdl[$i]->{end} + $j * 3;
								$flag = 1;
								last;
							}elsif( $codon =~ /[^ATGC]/){
								$extaa .= "X";
							}else{
								$extaa .= $transl_tab3{ $codon };
							}#else
						}else{#ref not mit seq
							if(exists $tt1_stop{ $codon } ) {#this is a stop codon
								$new_end = $pred_mdl[$i]->{end} + $j * 3;
								$flag = 1;
								last;
							}elsif ( $codon =~ /[^ATGC]/ ) {
								$extaa .= "X";
							}else{
								$extaa .= $transl_tab1{ $codon };
							}#else
						}#else
					}#foreach my $j
					if (not $flag) {#not found the stop codon
						$pred_mdl[$i]->{stopcodon_incomplete}++;
						#$new_end = $pred_mdl[$i]->{end};#20210917, do not change the coord
					}else{
						#revise the end (stop) of the model
						$pred_mdl[$i]->{end} = $new_end;##gene and mRNA, if available
						if (exists $pred_mdl[$i]->{exon_lst}){##exon, if available
							$pred_mdl[$i]->{exon_lst}->[-1]->{end} = $new_end;
						}#
						$pred_mdl[$i]->{CDS_lst}->[-1]->{end}      = $new_end;##CDS
					
						$tmpaa   .= $extaa;
						$tmpcds  .= $extcds;
						$tmpgene .=	$extgene;
						$tmpexon .= $extexon;
					}#else
						
					$mdl_stat = 2;
				}else{#"-"
					$flag = 0;#20210917
					$start_bound = 1;
					$cnt = (($pred_mdl[$i]->{start} - $start_bound) - ($pred_mdl[$i]->{start} - $start_bound)%3) / 3;
					if ($cnt > $max_ext){
						$cnt = $max_ext;
					}#if
					foreach my $j ( 1 .. $cnt){#how to process N???
						$codon     = rc ( uc ( substr $tgt_seq[ $tgt_name2idx{ $pred_mdl[$i]->{contig} } ], $pred_mdl[$i]->{start} - $j * 3 - 1, 3 ) );
						$extcds   .= $codon;
						$extgene  .= $codon;
						$extexon  .= $codon;
						
						if ( $mit_seq eq $gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{contig} ) {#ref mit seq, use tt3
							if( exists $tt3_stop{ $codon } ) {#this is a stop codon
								$new_start = $pred_mdl[$i]->{start} - $j * 3;
								$flag = 1;
								last;
							}elsif ($codon =~ /[^ATGC]/){
								$extaa .= "X";
							}else{
								$extaa .= $transl_tab3{ $codon };
							}#else
						}else{#ref not mit seq
							if(exists $tt1_stop{ $codon } ) {#this is a stop codon
								$new_start = $pred_mdl[$i]->{start} - $j * 3;
								$flag = 1;
								last;
							}elsif ($codon =~ /[^ATGC]/){
								$extaa .= "X";
							}else{
								$extaa .= $transl_tab1{ $codon };
							}#else
						}#else
					}#$j
					if (not $flag){
						$pred_mdl[$i]->{stopcodon_incomplete}++;
						$new_start = $pred_mdl[$i]->{start};#20210917, do not change the coord
					}else{#20210917
						#revise the start (stop) of the model
						$pred_mdl[$i]->{start} = $new_start;
						if(exists $pred_mdl[$i]->{exon_lst}){##
							$pred_mdl[$i]->{exon_lst}->[0]->{start} = $new_start;
						}#if
						$pred_mdl[$i]->{CDS_lst}->[0]->{start}      = $new_start;
					
						$tmpaa   .= $extaa;
						$tmpcds  .= $extcds;
						$tmpgene .=	$extgene;
						$tmpexon .= $extexon;
					}#else
						
					$mdl_stat = 2;
				}#else "-"
			}elsif (exists $pred_mdl[$i]->{startcodon_mutation}){#
				#try to extend, extend to the first possible start codon or to the end of the contig, but there should be no stop codon during the extend, 
				#or else to shrink to the first possible start codon (before meet the stop codon), only in the range of the first CDS
				if ($pred_mdl[$i]->{strand} eq "+"){#forward
					$flag = 0;
					$start_bound = 1;
					$cnt = (($pred_mdl[$i]->{start} - $start_bound) - ($pred_mdl[$i]->{start} - $start_bound)%3) / 3;
					if ($cnt > $max_ext){
						$cnt = $max_ext;
					}#if
					foreach my $j (1 .. $cnt){#how to process N???
						$codon     = uc ( substr $tgt_seq[ $tgt_name2idx{ $pred_mdl[$i]->{contig} } ], $pred_mdl[$i]->{start} - $j * 3 - 1, 3 );
						$extcds    = $codon . $extcds;
						$extgene   = $codon . $extgene;
						$extexon   = $codon . $extexon;
						
						if ( $mit_seq eq $gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{contig} ) {#ref mit seq, use tt3
							if( exists $tt3_stop{ $codon } ) {#this is a stop codon, meet stop codon first
								$flag = 1;#meet stop codon first
								last;
							}elsif ( exists $tt3_start{ $codon } ) {#meet start codon first
								$new_start = $pred_mdl[$i]->{start} - $j * 3;
								$flag = 2;#meet start codon first
								$extaa = $transl_tab3{ $codon } . $extaa;
								last;
							}elsif ( $codon =~ /[^ATGC]/ ) {
								$extaa = "X" . $extaa;
							}else{
								$extaa = $transl_tab3{ $codon } . $extaa;
							}#else
						}else{#ref not mit seq
							if( exists $tt1_stop{ $codon } ) {#this is a stop codon, meet stop codon first
								$flag = 1;#meet stop codon first
								last;
							}elsif ( exists $tt1_start{ $codon } ) {#meet start codon first
								$new_start = $pred_mdl[$i]->{start} - $j * 3;
								$flag = 2;#meet start codon first
								$extaa = $transl_tab1{ $codon } . $extaa;
								last;
							}elsif ( $codon =~ /[^ATGC]/ ) {
								$extaa = "X" . $extaa;
							}else{
								$extaa = $transl_tab1{ $codon } . $extaa;
							}#else
						}#else not mit seq
					}#foreac my $j
					if (not $flag){#not found the start codon, before reach the start of the contig
						$pred_mdl[$i]->{startcodon_incomplete} ++;
						#$new_start = $pred_mdl[$i]->{start};#20210917
					}elsif( $flag == 2){#
						#revise the start of the model
						$pred_mdl[$i]->{start} = $new_start;
						if (exists $pred_mdl[$i]->{exon_lst}){##exon, if available
							$pred_mdl[$i]->{exon_lst}->[0]->{start} = $new_start;
						}#
						$pred_mdl[$i]->{CDS_lst}->[0]->{start}      = $new_start;
						$tmpaa = $extaa . $tmpaa;
						$tmpcds = $extcds . $tmpcds;
						$tmpgene = $extgene . $tmpgene;
						$tmpexon = $extexon . $tmpexon;
						
						$mdl_stat = 3;
					}else{#flag == 1, meet the stop codon first, shrink
						$mdl_stat = 30;
						$pred_mdl[$i]->{startcodon_extend_problem} ++;						
						$tmpcds2  = $tmpcds;
						$tmpgene2 = $tmpgene;
						$tmpexon2 = $tmpexon;
						$tmpaa2   = $tmpaa;

						$flag = 0;
						$start_bound = $pred_mdl[$i]->{CDS_lst}->[0]->{end};	#the right boundary of the first CDS (already sorted based on coordinates)
						if ($start_bound - $pred_mdl[$i]->{start} < 5 ) {#20210917
							$cnt = 0;
						}else{
							$cnt = (($start_bound - 2 - $pred_mdl[$i]->{start}) - ($start_bound - 2 - $pred_mdl[$i]->{start})%3) / 3;
						}#else
						if ($cnt > $max_ext){
							$cnt = $max_ext;
						}#if
						foreach my $j ( 1 .. $cnt){#
							$codon = uc ( substr $tgt_seq[ $tgt_name2idx{ $pred_mdl[$i]->{contig} } ], $pred_mdl[$i]->{start} + $j * 3 - 1, 3 );
							if ( $mit_seq eq $gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{contig} ) {#ref mit seq, use tt3
								if (exists $tt3_start{$codon} ){#this is a start codon, stop here
									substr $tmpcds2, 0, 3, "";#delete last codon
									substr $tmpgene2, 0, 3, "";
									substr $tmpexon2, 0, 3, "";
									substr $tmpaa2, 0 , 1, "";
									$new_start = $pred_mdl[$i]->{start} + $j * 3;
									$flag = 5;#found the new start codon, and keep this codon
									last;
								}else{#delete last codon
									substr $tmpcds2, 0, 3, "";
									substr $tmpgene2, 0, 3, "";
									substr $tmpexon2, 0, 3, "";
									substr $tmpaa2, 0, 1, "";
								}#else
							}else{#ref not mit seq
								if (exists $tt1_start{$codon} ) {#
									substr $tmpcds2, 0, 3, "";#delete last codon
									substr $tmpgene2, 0, 3, "";
									substr $tmpexon2, 0, 3, "";
									substr $tmpaa2, 0 , 1, "";

									$new_start = $pred_mdl[$i]->{start} + $j * 3;
									$flag = 5;
									last;
								}else{
									substr $tmpcds2, 0, 3, "";
									substr $tmpgene2, 0, 3, "";
									substr $tmpexon2, 0, 3, "";
									substr $tmpaa2, 0, 1, "";
								}#else
							}#else
						}#foreach my $j
						if (not $flag){#not found the proper start codon,
							$pred_mdl[$i]->{startcodon_shrink_problem} ++;
							$mdl_stat = 32;
						}else{#shrink ok
							$pred_mdl[$i]->{startcodon_shrink_ok} ++;
							#revise the model
							$pred_mdl[$i]->{start} = $new_start;
							if(exists $pred_mdl[$i]->{exon_lst}){#
								$pred_mdl[$i]->{exon_lst}->[0]->{start} = $new_start;
							}#if
							$pred_mdl[$i]->{CDS_lst}->[0]->{start} = $new_start;
							$tmpaa = $tmpaa2;
							$tmpcds = $tmpcds2;
							$tmpgene = $tmpgene2;
							$tmpexon = $tmpexon2;
							
							$mdl_stat = 31;
						}#else
					}#else
				}else{#"-"
					$flag = 0;
					$end_bound = length $tgt_seq[ $tgt_name2idx{ $pred_mdl[$i]->{contig} } ];
					$cnt = (($end_bound - $pred_mdl[$i]->{end}) - ($end_bound - $pred_mdl[$i]->{end})%3) / 3;#max num of possible condons to extend available at the contig
					if ($cnt > $max_ext){
						$cnt = $max_ext;
					}#if
					foreach my $j (1 .. $cnt){#how to process N???
						$codon     = rc ( uc ( substr $tgt_seq[ $tgt_name2idx{ $pred_mdl[$i]->{contig} } ], $pred_mdl[$i]->{end} + ($j - 1) * 3, 3 ) );
						$extcds    = $codon . $extcds;
						$extgene   = $codon . $extgene;
						$extexon   = $codon . $extexon;
						
						if ($mit_seq eq $gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{contig} ) {#ref mit seq, use tt3
							if (exists $tt3_stop{ $codon } ) {#this is a stop codon, meet stop codon first
								$flag = 1;
								last;
							}elsif (exists $tt3_start{ $codon } ) {#meet start codon first
								$new_end = $pred_mdl[$i]->{end} + $j * 3;
								$flag = 2;
								$extaa = $transl_tab3{ $codon } . $extaa;
								last;
							}elsif ($codon =~ /[^ATGC]/ ) {
								$extaa = "X" . $extaa;
							}else{
								$extaa = $transl_tab3{ $codon } . $extaa;
							}#else
						}else{#ref not mit seq
							if( exists $tt1_stop{ $codon } ){#meet stop codon first
								$flag = 1;
								last;
							}elsif (exists $tt1_start{ $codon } ) {#
								$new_end = $pred_mdl[$i]->{end} + $j * 3;
								$flag =2;
								$extaa = $transl_tab1{ $codon } . $extaa;
								last;
							}elsif ($codon =~ /[^ATGC]/){
								$extaa = "X" . $extaa;
							}else{
								$extaa = $transl_tab1{ $codon } . $extaa;
							}#else
						}#else
					}#foreach my $j
					if (not $flag){#not found the start codon, before meet the end of the contig
						$pred_mdl[$i]->{startcodon_incomplete} ++;
						#$new_end = $pred_mdl[$i]->{end};#20210917
					}elsif($flag == 2){#
						#revise the start/end of the model
						$pred_mdl[$i]->{end}  = $new_end;
						if (exists $pred_mdl[$i]->{exon_lst}){#
							$pred_mdl[$i]->{exon_lst}->[-1]->{end} = $new_end;
						}#if
						$pred_mdl[$i]->{CDS_lst}->[-1]->{end}      = $new_end;
						$tmpaa = $extaa . $tmpaa;
						$tmpcds = $extcds . $tmpcds;
						$tmpgene = $extgene . $tmpgene;
						$tmpexon = $extexon . $tmpexon;
							
						$mdl_stat = 3;
					}else{
						$mdl_stat = 30;
						$pred_mdl[$i]->{startcodon_extend_problem} ++;
						#try to shrink
						$tmpcds2  = $tmpcds;
						$tmpgene2 = $tmpgene;
						$tmpexon2 = $tmpexon;
						$tmpaa2   = $tmpaa;						
						$flag = 0;
						$end_bound = $pred_mdl[$i]->{CDS_lst}->[-1]->{start};
						if ( $pred_mdl[$i]->{end} - $end_bound < 5 ) {#20210917
							$cnt = 0;
						}else{
							$cnt = (($pred_mdl[$i]->{end} - 2 - $end_bound) - ($pred_mdl[$i]->{end} - 2 - $end_bound)%3) / 3;
						}#else
						if ($cnt > $max_ext){
							$cnt = $max_ext;
						}#if
						foreach my $j ( 1 .. $cnt){#
							$codon = rc ( substr $tgt_seq[ $tgt_name2idx{ $pred_mdl[$i]->{contig} } ], $pred_mdl[$i]->{end} - $j * 3 - 2 - 1, 3 ) ;#20210917
							if ( $mit_seq eq $gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{contig} ) {#ref mit seq, use tt3
								if (exists $tt3_start{$codon} ) {#this is a start codon, stop here
									substr $tmpcds2, 0, 3, "";#delete last codon
									substr $tmpgene2, 0, 3, "";
									substr $tmpexon2, 0, 3, "";
									substr $tmpaa2, 0 , 1, "";
									$new_end = $pred_mdl[$i]->{end} - $j * 3;
									$flag = 5;
									last;
								}else{#delete last codon
									substr $tmpcds2, 0, 3, "";
									substr $tmpgene2, 0, 3, "";
									substr $tmpexon2, 0, 3, "";
									substr $tmpaa2, 0, 1, "";
								}#else
							}else{#ref not mit seq
								if (exists $tt1_start{$codon} ) {#
									substr $tmpcds2, 0, 3, "";#delete last codon
									substr $tmpgene2, 0, 3, "";
									substr $tmpexon2, 0, 3, "";
									substr $tmpaa2, 0 , 1, "";
									$new_end = $pred_mdl[$i]->{end} - $j * 3;#20210917
									$flag = 5;
									last;
								}else{
									substr $tmpcds2, 0, 3, "";
									substr $tmpgene2, 0, 3, "";
									substr $tmpexon2, 0, 3, "";
									substr $tmpaa2, 0, 1, "";
								}#else
							}#else
						}#foreach my $j
						if (not $flag){#not found the proper start codon
							$pred_mdl[$i]->{startcodon_shrink_problem} ++;
							$mdl_stat = 32;
						}else{#shrink ok
							$pred_mdl[$i]->{startcodon_shrink_ok} ++;
							#revise the model
							$pred_mdl[$i]->{end} = $new_end;
							if(exists $pred_mdl[$i]->{exon_lst}){#
								$pred_mdl[$i]->{exon_lst}->[-1]->{end} = $new_end;
							}#if
							$pred_mdl[$i]->{CDS_lst}->[-1]->{end} = $new_end;#20210917
							$tmpaa = $tmpaa2;
							$tmpcds = $tmpcds2;
							$tmpgene = $tmpgene2;
							$tmpexon = $tmpexon2;
							
							$mdl_stat = 31;
						}#else shrink ok
					}#else
				}#else "-"
			}else{#model ok, do nothing
			}#else
				
			
			#save this gene model
			$tmpid = sprintf "%06d", ++$tgt_gid;
			
			#output
			if ($pred_mdl[$i]->{premature}){#pseudo gene
				$method = "blastn_chk";
				print PGFF $pred_mdl[$i]->{contig},"\t",
				$method,"\t",
				"pseudogene","\t",
				$pred_mdl[$i]->{start},"\t",
				$pred_mdl[$i]->{end},"\t",
				".","\t",
				$pred_mdl[$i]->{strand},"\t",
				".","\t",
				"ID=",$pref, "_G",$tmpid,";",
				"Parent=",$pref, "_G",$tmpid,";",
				"Name=",$pref, "_G",$tmpid,";",
				"locus_tag=",$pref, "_G",$tmpid,";",
				"ref_gene_mdl=","ID:",$pred_mdl[$i]->{ref_gene_id},",",
				"Name:",$gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{Name},",",
				"locus_tag:",$pred_mdl[$i]->{ref_locus_tag},"\n";
	
				#mRNA and each exon, if exists
				if (exists $pred_mdl[$i]->{exon_lst}){
					print PGFF $pred_mdl[$i]->{contig},"\t",
					$method,"\t",
					"mRNA","\t",
					$pred_mdl[$i]->{start},"\t",
					$pred_mdl[$i]->{end},"\t",
					".","\t",
					$pred_mdl[$i]->{strand},"\t",
					".","\t",
					"ID=",$pref, "_T",$tmpid,".1;",
					"Parent=",$pref, "_G",$tmpid,";",
					"Name=",$pref, "_T",$tmpid,".1;",
					"locus_tag=",$pref, "_G",$tmpid,";",
					"product=",$pred_mdl[$i]->{exon_lst}->[0]->{product},";",
					"ref_gene_mdl=","ID:",$pred_mdl[$i]->{ref_gene_id},",",
					"Name:",$gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{Name},",",
					"locus_tag:",$pred_mdl[$i]->{ref_locus_tag};
					if ($pred_mdl[$i]->{ref_transcript_id}){
						print PGFF ",transcript_id:",$pred_mdl[$i]->{ref_transcript_id},"\n";
					}else{
						print PGFF "\n";
					}#else		
	
					foreach my $j (0 .. $#{$pred_mdl[$i]->{exon_lst}}){#each exon
						print PGFF $pred_mdl[$i]->{contig},"\t",
						$method,"\t",
						"exon","\t",
						$pred_mdl[$i]->{exon_lst}->[$j]->{start},"\t",
						$pred_mdl[$i]->{exon_lst}->[$j]->{end},"\t",
						".","\t",
						$pred_mdl[$i]->{strand},"\t",
						".","\t",
						"ID=",$pref, "_T",$tmpid,".1","-",$j+1,";",
						"Parent=",$pref, "_T",$tmpid,".1",";",
						"locus_tag=",$pref, "_G",$tmpid,";",
						"product=",$pred_mdl[$i]->{exon_lst}->[$j]->{product},";",
						"transcript_id=",$pref, "_T",$tmpid,".1",";",
						"ref_gene_mdl=","ID:",$pred_mdl[$i]->{ref_gene_id},",",
						"Name:",$gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{Name},",",
						"locus_tag:",$pred_mdl[$i]->{ref_locus_tag};
						if ($pred_mdl[$i]->{ref_transcript_id}){
							print PGFF ",transcript_id:",$pred_mdl[$i]->{ref_transcript_id},"\n";
						}else{
							print PGFF "\n";
						}#else
					}#foreach 
				}#if
	
				#each cds
				foreach my $j (0 .. $#{$pred_mdl[$i]->{CDS_lst}}){#each cds
					print PGFF $pred_mdl[$i]->{contig},"\t",
					$method,"\t",
					"CDS","\t",
					$pred_mdl[$i]->{CDS_lst}->[$j]->{start},"\t",
					$pred_mdl[$i]->{CDS_lst}->[$j]->{end},"\t",
					".","\t",
					$pred_mdl[$i]->{strand},"\t",
					$pred_mdl[$i]->{CDS_lst}->[$j]->{phase},"\t",
					"ID=",$pref, "_C",$tmpid,".1","-", $j+1,";";
					if (exists $pred_mdl[$i]->{exon_lst}){
						print PGFF "Parent=",$pref, "_T",$tmpid,".1",";";
					}else{
						print PGFF "Parent=",$pref, "_G",$tmpid,";";
					}#else
					print PGFF "locus_tag=",$pref, "_G",$tmpid,";",
					"product=",$pred_mdl[$i]->{CDS_lst}->[$j]->{product},";",
					"protein_id=",$pref, "_P",$tmpid,".1",";";
					if ($pred_mdl[$i]->{CDS_lst}->[$j]->{Note} ) {
						print PGFF "Note=",$pred_mdl[$i]->{CDS_lst}->[$j]->{Note},";";
					}#if
					print PGFF "ref_gene_mdl=","ID:",$pred_mdl[$i]->{ref_gene_id},",",
					"Name:",$gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{Name},",",
					"locus_tag:",$gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{locus_tag},",",
					"protein_id:",$pred_mdl[$i]->{CDS_lst}->[$j]->{ref_protein_id};
					if (exists $pred_mdl[$i]->{CDS_lst}->[$j]->{transl_table}){
						print PGFF ";transl_table=",$pred_mdl[$i]->{CDS_lst}->[$j]->{transl_table},"\n";
					}else{
						print PGFF "\n";
					}#else
				}# $j
			
				print PGENE ">",$pref, "_G",$tmpid,"\n",$tmpgene,"\n";
				print PMRNA ">",$pref, "_T",$tmpid,".1\n",$tmpexon,"\n" if $tmpexon;
				print PCDS  ">",$pref, "_C",$tmpid,".1\n",$tmpcds,"\n";
				print PPROT ">",$pref, "_P",$tmpid,".1\n",$tmpaa,"\n";
				
			}else{#
				$ch1 = 0;
				$ch2 = 0;
				if (exists $pred_mdl[$i]->{startcodon_shrink_problem} ) {
					$ch1 = 4;
				}elsif (exists $pred_mdl[$i]->{startcodon_shrink_ok} ) {
					$ch1 = 3;
				}elsif (exists $pred_mdl[$i]->{startcodon_incomplete} ) {
					$ch1 = 2;
				}elsif(exists $pred_mdl[$i]->{startcodon_mutation} ){
					$ch1 = 1;
				}#
				if (exists $pred_mdl[$i]->{stopcodon_incomplete} ) {
					$ch2 = 2;
				}elsif (exists $pred_mdl[$i]->{stopcodon_mutation}){
					$ch2 = 1;
				}#
				if ($ch1 == 0 and $ch2 == 0){
					$method = "blastn";
				}else{
					$method = "blastn_ext" . $ch1 . $ch2;
				}#else

				print GFF $pred_mdl[$i]->{contig},"\t",
				$method,"\t",
				"gene","\t",
				$pred_mdl[$i]->{start},"\t",
				$pred_mdl[$i]->{end},"\t",
				".","\t",
				$pred_mdl[$i]->{strand},"\t",
				".","\t",
				"ID=",$pref, "_G",$tmpid,";",
				"Parent=",$pref, "_G",$tmpid,";",
				"Name=",$pref, "_G",$tmpid,";",
				"locus_tag=",$pref, "_G",$tmpid,";",
				"ref_gene_mdl=","ID:",$pred_mdl[$i]->{ref_gene_id},",",
				"Name:",$gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{Name},",",
				"locus_tag:",$pred_mdl[$i]->{ref_locus_tag},"\n";
	
				#mRNA and each exon, if exists
				if (exists $pred_mdl[$i]->{exon_lst}){
					print GFF $pred_mdl[$i]->{contig},"\t",
					$method,"\t",
					"mRNA","\t",
					$pred_mdl[$i]->{start},"\t",
					$pred_mdl[$i]->{end},"\t",
					".","\t",
					$pred_mdl[$i]->{strand},"\t",
					".","\t",
					"ID=",$pref, "_T",$tmpid,".1;",
					"Parent=",$pref, "_G",$tmpid,";",
					"Name=",$pref, "_T",$tmpid,".1;",
					"locus_tag=",$pref, "_G",$tmpid,";",
					"product=",$pred_mdl[$i]->{exon_lst}->[0]->{product},";",
					"ref_gene_mdl=","ID:",$pred_mdl[$i]->{ref_gene_id},",",
					"Name:",$gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{Name},",",
					"locus_tag:",$pred_mdl[$i]->{ref_locus_tag};
					if ($pred_mdl[$i]->{ref_transcript_id}){
						print GFF ",transcript_id:",$pred_mdl[$i]->{ref_transcript_id},"\n";
					}else{
						print GFF "\n";
					}#else		
	
					foreach my $j (0 .. $#{$pred_mdl[$i]->{exon_lst}}){#each exon
						print GFF $pred_mdl[$i]->{contig},"\t",
						$method,"\t",
						"exon","\t",
						$pred_mdl[$i]->{exon_lst}->[$j]->{start},"\t",
						$pred_mdl[$i]->{exon_lst}->[$j]->{end},"\t",
						".","\t",
						$pred_mdl[$i]->{strand},"\t",
						".","\t",
						"ID=",$pref, "_T",$tmpid,".1","-",$j+1,";",
						"Parent=",$pref, "_T",$tmpid,".1",";",
						"locus_tag=",$pref, "_G",$tmpid,";",
						"product=",$pred_mdl[$i]->{exon_lst}->[$j]->{product},";",
						"transcript_id=",$pref, "_T",$tmpid,".1",";",
						"ref_gene_mdl=","ID:",$pred_mdl[$i]->{ref_gene_id},",",
						"Name:",$gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{Name},",",
						"locus_tag:",$pred_mdl[$i]->{ref_locus_tag};
						if ($pred_mdl[$i]->{ref_transcript_id}){
							print GFF ",transcript_id:",$pred_mdl[$i]->{ref_transcript_id},"\n";
						}else{
							print GFF "\n";
						}#else
					}#foreach 
				}#if
	
				#each cds
				foreach my $j (0 .. $#{$pred_mdl[$i]->{CDS_lst}}){#each cds
					print GFF $pred_mdl[$i]->{contig},"\t",
					$method,"\t",
					"CDS","\t",
					$pred_mdl[$i]->{CDS_lst}->[$j]->{start},"\t",
					$pred_mdl[$i]->{CDS_lst}->[$j]->{end},"\t",
					".","\t",
					$pred_mdl[$i]->{strand},"\t",
					$pred_mdl[$i]->{CDS_lst}->[$j]->{phase},"\t",
					"ID=",$pref, "_C",$tmpid,".1","-", $j+1,";";
					if (exists $pred_mdl[$i]->{exon_lst}){
						print GFF "Parent=",$pref, "_T",$tmpid,".1",";";
					}else{
						print GFF "Parent=",$pref, "_G",$tmpid,";";
					}#else
					print GFF "locus_tag=",$pref, "_G",$tmpid,";",
					"product=",$pred_mdl[$i]->{CDS_lst}->[$j]->{product},";",
					"protein_id=",$pref, "_P",$tmpid,".1",";";
					if ($pred_mdl[$i]->{CDS_lst}->[$j]->{Note} ) {
						print GFF "Note=",$pred_mdl[$i]->{CDS_lst}->[$j]->{Note},";";
					}#if
					print GFF "ref_gene_mdl=","ID:",$pred_mdl[$i]->{ref_gene_id},",",
					"Name:",$gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{Name},",",
					"locus_tag:",$gene_info{ $pred_mdl[$i]->{ref_gene_id} }->{locus_tag},",",
					"protein_id:",$pred_mdl[$i]->{CDS_lst}->[$j]->{ref_protein_id};
					if (exists $pred_mdl[$i]->{CDS_lst}->[$j]->{transl_table}){
						print GFF ";transl_table=",$pred_mdl[$i]->{CDS_lst}->[$j]->{transl_table},"\n";
					}else{
						print GFF "\n";
					}#else
				}# $j
			
				print GENE ">",$pref, "_G",$tmpid,"\n",$tmpgene,"\n";
				print MRNA ">",$pref, "_T",$tmpid,".1\n",$tmpexon,"\n" if $tmpexon;
				print CDS  ">",$pref, "_C",$tmpid,".1\n",$tmpcds,"\n";
				print PROT ">",$pref, "_P",$tmpid,".1\n",$tmpaa,"\n";
			}#else not pseudogene
		}#if exists $CDS
	}# 
	close GFF;
	close GENE;
	close MRNA;
	close PROT;
	close CDS;
	close PGFF;
	close PGENE;
	close PMRNA;
	close PPROT;
	close PCDS;


	close OUT;
	
	#not process pseudogene			
	# output report on mutations on start codon /stop codon
	
	#based on the range of gene mdls
	#get the coor of inter-mdl regions, discard too short regions (<300 nt)
	my %inter_mdl;	#{$contig}->[$idx]->[0/1] start end
	my %inter_mdl_num;#{$contig} => num
	my @tgt_seq_mask;
	foreach my $i (0 .. $#tgt_seq){
		$tgt_seq_mask[ $i ] = uc $tgt_seq[$i];
	}#
	
	my $tmpidx;
	foreach my $i (sort { $pred_mdl[$a]->{contig} cmp $pred_mdl[$b]->{contig}
	or $pred_mdl[$a]->{start} <=> $pred_mdl[$b]->{start} 
	or $pred_mdl[$a]->{end}   <=> $pred_mdl[$b]->{end} } (0 .. $#pred_mdl) ){
		next if exists $pred_mdl[$i]->{premature};#exclude pseudogenes
		$s2 = $pred_mdl[$i]->{start};
		$e2 = $pred_mdl[$i]->{end};
		$tmpidx = $tgt_name2idx{ $pred_mdl[$i]->{contig} }; 
		$tmpseq = substr $tgt_seq_mask[ $tmpidx ], $s2 - 1, $e2 - $s2 + 1;
		$tmpseq = lc $tmpseq;
		substr $tgt_seq_mask[ $tmpidx ], $s2 - 1, $e2 - $s2 + 1, $tmpseq;
	}#foreach my $i
	
	# output the sequence for the segments without homologous gene models
	open NOMDL1, ">$outnomdl1" or die "fail to open $outnomdl1:$!\n";
	foreach my $i (0 .. $#tgt_seq_mask){
		print NOMDL1 ">",$tgt_seqname[$i],"_masked","\n",$tgt_seq_mask[$i],"\n";
	}#
	close NOMDL1;
	
#	open NOMDL2, ">$outnomdl2" or die "fail to open $outnomdl2:$!\n";
#	@items2 = sort { $pred_mdl[$a]->{contig} cmp $pred_mdl[$b]->{contig}
#	or $pred_mdl[$a]->{start} <=> $pred_mdl[$b]->{start} 
#	or $pred_mdl[$a]->{end}   <=> $pred_mdl[$b]->{end} } (0 .. $#pred_mdl);
#	@items = ();
#	foreach my $i (@items2){
#		if (not exists $pred_mdl[$i]->{premature}){
#			push @items, $i;
#		}#if
#	}#
#	my $tmpctg = "";
#	$s2 = -1;
#	$e2 = -1;
#	foreach my $i (0 .. $#items){
#		next if exists $pred_mdl[ $items[$i] ]->{premature};#exclude pseudogenes!!!!!!!!!!!!!!!!!
#		if ($tmpctg ne $pred_mdl[$items[$i]]->{contig} ){#this is a new contig
#			if ($tmpctg eq ""){#this is the first contig
#				$tmpctg = $pred_mdl[$items[$i]]->{contig};
#				$s2 = 1;
#				$e2 = $pred_mdl[$items[$i]]->{start} - 1;
#				if ($e2 - $s2 > 300) {
#					#save this seg
#					$tmpseq = substr $tgt_seq_mask[ $tgt_name2idx{ $tmpctg } ], $s2 - 1, $e2 - $s2 + 1;
#					print NOMDL2 ">", $tmpctg, "_", $s2, "_", $e2,"\n",$tmpseq,"\n";
#				}#if
#				
#			}else{#not the first mdl of a contig
#				#save the last seg of the last contig
#				$s2 = $pred_mdl[ $items[$i-1] ]->{end} + 1;
#				$e2 = length $tgt_seq_mask[ $tgt_name2idx{ $tmpctg } ];
#				if ($e2 - $s2 > 300){
#					#save last seg on last ctg
#					$tmpseq = substr $tgt_seq_mask[ $tgt_name2idx{ $tmpctg } ], $s2 - 1, $e2 - $s2 + 1;
#					print NOMDL2 ">",$tmpctg, "_", $s2, "_", $e2,"\n",$tmpseq,"\n";
#				}#if
#				#save the first seg of the this contig
#				$tmpctg = $pred_mdl[$items[$i]]->{contig};
#				$s2 = 1;
#				$e2 = $pred_mdl[$items[$i]]->{start} - 1;
#				if ($e2 - $s2 > 300){
#					$tmpseq = substr $tgt_seq_mask[ $tgt_name2idx{ $tmpctg } ], $s2 - 1, $e2 - $s2 + 1;
#					print NOMDL2 ">",$tmpctg,"_",$s2,"_",$e2,"\n",$tmpseq,"\n";
#				}#if
#			}#else
#		}else{#not a new contig
#			$s2 = $pred_mdl[$items[$i-1]]->{end} + 1;
#			$e2 = $pred_mdl[$items[$i]]->{start} - 1;
#			if ($e2 - $s2 > 300){
#				$tmpseq = substr $tgt_seq_mask[ $tgt_name2idx{ $tmpctg } ], $s2 - 1, $e2 - $s2 + 1;
#				print NOMDL2 ">", $tmpctg, "_",$s2,"_",$e2,"\n",$tmpseq,"\n";
#			}#if
#		}#else
#	}#
#	#check the last mdl of the last ctg
#	$s2 = $pred_mdl[$items[-1]]->{end} + 1;
#	$e2 = length $tgt_seq_mask[ $tgt_name2idx{ $tmpctg } ];
#	if ($e2 - $s2 > 300 ){
#		$tmpseq = substr $tgt_seq_mask[ $tgt_name2idx{ $tmpctg } ], $s2 - 1, $e2 - $s2 + 1;
#		print NOMDL2 ">",$tmpctg,"_",$s2,"_",$e2,"\n",$tmpseq,"\n";
#	}#if
#	close NOMDL2;

}#sub	
	



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
	
	#remove the whitespaces in the seqs? YES
	#remove the asterisk in the seqs? NO
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


	