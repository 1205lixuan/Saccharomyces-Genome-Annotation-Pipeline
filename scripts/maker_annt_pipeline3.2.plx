#!/usr/bin/perl
#maker_annt_pipeline3.2.plx
use strict;
use warnings;
#also convert non-ATGC nt to N
sub ReadFasta($\@\@);#$filename,@seqname,@seq
sub WriteFasta($\@\@);#$filename,@seqname,@seq
sub rc($);#$inseq

sub ReadFixGff($$$$$);#$in_gff_filename,$out_prefix,$id_base,$ref_fa,$seq_dir

my $id_base = 10000;

my $fa_suffix  = ".blastn.no_mdl_mask.fna";
my $new_fa_suffix = ".blastn.no_mdl_mask.hardmask.fna";
my $pred_dir_suffix1 = ".blastn.no_mdl_mask.hardmask.maker.output";
#my $pred_dir_suffix2 = ".blastn.no_mdl_mask.hardmask_datastore";
my $pred_log_suffix  = ".blastn.no_mdl_mask.hardmask_master_datastore_index.log";
my $cpu_num = shift;
my $seq_dir = shift;
my $infile  = shift;
if (not $infile){
	die "usage: perl maker_annt_pipeline3.2.plx <cpu_num> <seq_dir> <query_fa_list>\n";
}#if
my $outfile = $infile . ".log";
my @items;
my @items2;
my @da1;
my @da2;
my @da3;
my @da4;
my $dir;
my $pref;
my $lbuf;
my $tmp;

my @seqname;
my @seq;
my @ch;

my %seq_stat;
my %result_dir;

my $cnt = 0;

open JOB,"$infile" or die "fail to open $infile:$!\n";
#infile format: target_fasta_name	prefix
while(<JOB>){
	chomp;
	next if not $_;
	next if /^#/;
	$cnt ++;
	@items = split /\t/,$_;
	print "\n\nJOB ID $cnt: $items[0]\t$items[1]\n";
	#convert to hard mask
	@seqname = ();
	@seq     = ();
	ReadFasta("$seq_dir/$items[0]$fa_suffix",@seqname,@seq);
	foreach my $i (0 .. $#seq){
		@ch = split //,$seq[$i];
		foreach my $j (0 .. $#ch){
			if ($ch[$j] eq "a" or $ch[$j] eq "t" or $ch[$j] eq "g" or $ch[$j] eq "c"){
				$ch[$j] = "n";
			}#
			if ($ch[$j] ne "A" and $ch[$j] ne "T" and $ch[$j] ne "G" and $ch[$j] ne "C" and $ch[$j] ne "N" 
			and $ch[$j] ne "a" and $ch[$j] ne "t" and $ch[$j] ne "g" and $ch[$j] ne "c" and $ch[$j] ne "n" ) {#unsupported nt
				$ch[$j] = "N";
			}#if
		}#
		$seq[$i] = join "",@ch;
	}#
	WriteFasta("$seq_dir/$items[0]$new_fa_suffix",@seqname,@seq);
				
	#maker
	system("maker -genome $seq_dir/$items[0]$new_fa_suffix -cpus $cpu_num maker_opts_S288C.ctl maker_bopts.ctl maker_exe.ctl > $seq_dir/$items[0].maker.log 2>&1");
	#check whether correctly completed
	%seq_stat = ();
	$tmp = $seq_dir . "/" . $items[0] . $pred_dir_suffix1 . "/" . $items[0] . $pred_log_suffix;
	%result_dir=();
	open ST, $tmp or die "fail to open $tmp:$!\n";
	while (<ST>){
		chomp;
		@items2 = split /\t/, $_;
		$seq_stat{ $items2[0] }->{ $items2[2] } ++;
		$result_dir{ $items2[1] } ++;
	}
	close ST;
	foreach my $i (0 .. $#seqname){
		if(exists $seq_stat{ $seqname[$i] }->{STARTED} and exists $seq_stat{ $seqname[$i] }->{FINISHED}){#ok
			#do nothing
		}else{
			print "warning: MAKER unfinished for $items[0] seq $seqname[$i]\n";
		}#else
	}#$i
		
	#collect results, .gff, prot.faa, transcripts.faa
	if (-e "$seq_dir/$items[0].maker.origin_pick.gff"){
		system ("echo \"\" > $seq_dir/$items[0].maker.origin_pick.gff");
	}#
	foreach my $i (sort keys %result_dir){
		system ("grep -e \"	maker	\" $seq_dir/$items[0]$pred_dir_suffix1/$i*.gff >> $seq_dir/$items[0].maker.origin_pick.gff");
	}#

	#and rename seq id ?
	#based on gff file, and produce other files
	ReadFixGff( "$seq_dir/$items[0].maker.origin_pick.gff", $items[0], $id_base, "$items[0]$new_fa_suffix",$seq_dir);
	
	
	#and merge with blastn annot
	
	
}#
close JOB;









sub ReadFixGff($$$$$){#$in_gff_filename,$out_prefix,$id_base,$ref_fa,$seq_dir
	my $in_gff = shift;
	my $pref = shift;
	my $id_base= shift;
	my $ref_fa  = shift;
	my $seq_dir = shift;

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

	
	my $lno = 0;
	my @items;
	my @items2;
	my @items3;
	my @order;
	my %tmph;
	my $cnt;

	my %IDtype;		#{$ID} => type
	my %gene_info;	#{$gene_ID}->{ID/Name/contig/start/end/score/strand/phase/mRNA_num/mRNA_lst} {mRNA_lst}->[$idx] => mRNA_ID 
	my %mRNA_info;	#{$mRNA_ID}->{ID/Name/Parent/exon_lst/CDS_lst} {exon/CDS_lst}->[$idx] => exon/CDS ID
	my %exon_info;	#{$exon_ID}->{ID/Name/Parent/idx_in_exon_lst/contig/start/codon/score/strand/phase}
	my %CDS_info;	#{$CDS_ID}->{ID/Name/Parent/idx_in_CDS_lst/contig/start/codon/score/strand/phase}

	my @ref_seqname;
	my @ref_seq;
	my %name2seq;

	my $tmpid;
	my $tmpcds;
	my $tmpexon;
	my $tmpaa;
	my $tmpgene;
	my $tmpseq;
	
	my $codon;
	
	ReadFasta( "$seq_dir/$ref_fa", @ref_seqname, @ref_seq );
	foreach my $i ( 0 .. $#ref_seqname){
		$name2seq{ $ref_seqname[$i] } = $ref_seq[$i];
	}#
	
	#read models from gff
	open ReadFixGff_IN, "$in_gff" or die "fail to open $in_gff:$!\n";
	while (<ReadFixGff_IN>){
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
	#		die "duplicated ID $tmph{ID} at line $lno in $in_gff.\n";
		}#else
		if ($items[6] ne "+" and $items[6] ne "-"){
			die "strand not specified $items[6] at line $lno in $in_gff.\n";
		}#if
		if ($items[2] eq "gene"){
			if (exists $gene_info{ $tmph{ID} } ) {#
	#			die "duplicated gene ID at line $lno in $in_gff.\n";
			}#if
			$gene_info{ $tmph{ ID } }->{ID} = $tmph{ ID };
			$gene_info{ $tmph{ ID } }->{Name}=$tmph{ Name };
			$gene_info{ $tmph{ ID } }->{contig} =$items[0];
			$gene_info{ $tmph{ ID } }->{start}=$items[3];
			$gene_info{ $tmph{ ID } }->{end} = $items[4];
			$gene_info{ $tmph{ ID } }->{score}= $items[5];
			$gene_info{ $tmph{ ID } }->{strand}= $items[6];
			$gene_info{ $tmph{ ID } }->{phase}= $items[7];
#			$gene_info{ $tmph{ ID } }->{gene} = $tmph{ gene };
#			$gene_info{ $tmph{ ID } }->{locus_tag} = $tmph{locus_tag};
			$gene_info{ $tmph{ ID } }->{mRNA_num} = 0;
		}elsif ($items[2] eq "mRNA"){
			if (exists $mRNA_info{ $tmph{ID} } ) {#
	#			die "duplicated mRNA ID at line $lno in $in_gff.\n";
			}#if
			$mRNA_info{ $tmph{ ID } }->{ID} = $tmph{ ID };
			$mRNA_info{ $tmph{ ID } }->{Name}=$tmph{ Name };
			$mRNA_info{ $tmph{ ID } }->{contig} =$items[0];
			$mRNA_info{ $tmph{ ID } }->{start}=$items[3];
			$mRNA_info{ $tmph{ ID } }->{end} = $items[4];
			$mRNA_info{ $tmph{ ID } }->{score}= $items[5];
			$mRNA_info{ $tmph{ ID } }->{strand}= $items[6];
			$mRNA_info{ $tmph{ ID } }->{phase}= $items[7];
#			$mRNA_info{ $tmph{ ID } }->{gene} = $tmph{ gene };
#			$mRNA_info{ $tmph{ ID } }->{transcript_id} = $tmph{ transcript_id};
			$mRNA_info{ $tmph{ ID } }->{exon_num} = 0;
			$mRNA_info{ $tmph{ ID } }->{CDS_num} = 0;
#			$mRNA_info{ $tmph{ ID } }->{locus_tag} = $tmph{ locus_tag };
#			$mRNA_info{ $tmph{ ID } }->{product}   = $tmph{ product };
			
			$mRNA_info{ $tmph{ ID } }->{Parent} = $tmph{ Parent };
			if(exists $gene_info{ $tmph{ Parent } } ) {#because other features, e.g. pseudogene also contains mRNA
				push @{ $gene_info{ $tmph{ Parent } }->{mRNA_lst} }, $tmph{ ID };
				$gene_info{ $tmph{ Parent } }->{mRNA_num} ++;
			}else{
	#			print "NOPARENT_mRNA","\t","PARENTtype:",$IDtype{ $tmph{ Parent } },"\t",$lno,"\t",$_,"\n";
			}
		}elsif ($items[2] eq "exon"){
			if (exists $exon_info{ $tmph{ID} } ) {#
	#			die "duplicated exon ID at line $lno in $in_gff.\n";
			}#if		
			$exon_info{ $tmph{ ID } }->{ID} = $tmph{ ID };
			$exon_info{ $tmph{ ID } }->{Name}=$tmph{ Name };
			$exon_info{ $tmph{ ID } }->{contig} =$items[0];
			$exon_info{ $tmph{ ID } }->{start}=$items[3];
			$exon_info{ $tmph{ ID } }->{end} = $items[4];
			$exon_info{ $tmph{ ID } }->{score}= $items[5];
			$exon_info{ $tmph{ ID } }->{strand}= $items[6];
			$exon_info{ $tmph{ ID } }->{phase}= $items[7];
#			$exon_info{ $tmph{ ID } }->{gene} = $tmph{ gene };
#			$exon_info{ $tmph{ ID } }->{locus_tag} = $tmph{ locus_tag };
#			$exon_info{ $tmph{ ID } }->{product}   = $tmph{ product };
#			$exon_info{ $tmph{ ID } }->{transcript_id} = $tmph{ transcript_id};
			$exon_info{ $tmph{ ID } }->{Parent} = $tmph{ Parent };
			if(exists $mRNA_info{ $tmph{ Parent } } ) {#other features, e.g. ncRNA, also contains exons
				push @{ $mRNA_info{ $tmph{ Parent } }->{exon_lst} }, $tmph{ ID };
				$mRNA_info{ $tmph{ Parent } }->{exon_num} ++;
			}else{
	#			print "NOPARENT_exon","\t","PARENTtype:",$IDtype{ $tmph{ Parent } },"\t",$lno,"\t",$_,"\n";
			}
	
		}elsif ($items[2] eq "CDS"){
			$tmph{ID} .= "__$lno";	#ID are identical for different CDS of the same gene in maker output gff file, maybe a bug
			if (exists $CDS_info{ $tmph{ID} } ) {#
	#			die "duplicated CDS ID at line $lno in $in_gff.\n";
			}#if
			$CDS_info{ $tmph{ ID } }->{ID} = $tmph{ ID };
			$CDS_info{ $tmph{ ID } }->{Name}=$tmph{ Name };
			$CDS_info{ $tmph{ ID } }->{contig} =$items[0];
			$CDS_info{ $tmph{ ID } }->{start}=$items[3];
			$CDS_info{ $tmph{ ID } }->{end} = $items[4];
			$CDS_info{ $tmph{ ID } }->{score}= $items[5];
			$CDS_info{ $tmph{ ID } }->{strand}= $items[6];
			$CDS_info{ $tmph{ ID } }->{phase}= $items[7];
			$CDS_info{ $tmph{ ID } }->{gene} = $tmph{ gene };
#			$CDS_info{ $tmph{ ID } }->{protein_id} = $tmph{ protein_id};
#			$CDS_info{ $tmph{ ID } }->{locus_tag} = $tmph{ locus_tag };
#			$CDS_info{ $tmph{ ID } }->{product}   = $tmph{ product };
#			$CDS_info{ $tmph{ ID } }->{Note}      = $tmph{ Note };
			$CDS_info{ $tmph{ ID } }->{Parent} = $tmph{ Parent };
			if(exists $tmph{ transl_table } ) {
				$CDS_info{ $tmph{ ID } }->{transl_table} = $tmph{ transl_table };
				if ($tmph{transl_table} != 3){
					print "warning: unrecognized transl_table $tmph{ transl_table } at line $lno in $in_gff, please specify.\n";
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
	close ReadFixGff_IN;
	
	#fix id and output seqs
	my $outgff   = $seq_dir . "/" . $pref . ".maker.gff";
	my $outgene  = $seq_dir . "/" . $pref . ".maker.gene.fna";
	my $outmrna  = $seq_dir . "/" . $pref . ".maker.mrna.fna";
	my $outprot  = $seq_dir . "/" . $pref . ".maker.prot.faa";
	my $outcds   = $seq_dir . "/" . $pref . ".maker.cds.fna";
	open ReadFixGff_GFF, "|pigz >$outgff.gz" or die "fail to open $outgff:$!\n";
	open ReadFixGff_GENE, "|pigz >$outgene.gz" or die "fail to open $outgene:$!\n";
	open ReadFixGff_MRNA, "|pigz >$outmrna.gz" or die "fail to open $outmrna:$!\n";
	open ReadFixGff_PROT, "|pigz >$outprot.gz" or die "fail to open $outprot:$!\n";
	open ReadFixGff_CDS,  "|pigz >$outcds.gz" or die "fail to open $outcds:$!\n";

	$cnt = $id_base;
	foreach my $i (sort { $gene_info{ $a }->{contig} cmp $gene_info{$b}->{contig}
	or $gene_info{ $a }->{start} <=> $gene_info{ $b }->{start}
	or $gene_info{ $a }->{end}   <=> $gene_info{ $b }->{end} } keys %gene_info ) {#each gene mdl
		$tmpid = sprintf "%06d", ++$cnt;
		$tmpcds = "";#for CDS, for some genes, CDS is different from exon
		$tmpexon = "";#for (exon of) mRNA (exon of mRNA and cds are almost identical here, UTRs are ignored)
		$tmpaa  = "";
		$tmpgene = "";
		
		#get the gene mdl
		$tmpgene = substr $name2seq{ $gene_info{ $i }->{contig} }, $gene_info{ $i }->{start} - 1, $gene_info{ $i }->{end} - $gene_info{ $i }->{start} + 1;
		if ($gene_info{$i}->{strand} eq "-"){
			$tmpgene = rc($tmpgene);
		}elsif ($gene_info{$i}->{strand} eq "+"){
			#correct, do nothing
		}else{
			die "wrong strand.\n";
		}#else
		print ReadFixGff_GFF $gene_info{ $i }->{contig},"\t",
		"maker","\t",
		"gene","\t",
		$gene_info{ $i }->{start},"\t",
		$gene_info{ $i }->{end},"\t",
		$gene_info{ $i }->{score},"\t",
		$gene_info{ $i }->{strand},"\t",
		".","\t",
		"ID=",$pref, "_G",$tmpid,";",
		"Parent=",$pref, "_G",$tmpid,";",
		"Name=",$pref, "_G",$tmpid,";",
		"locus_tag=",$pref, "_G",$tmpid,"\n";
		
		if ($gene_info{ $i }->{ mRNA_num } ) {#
			#suppose contain only one mRNA, no AS
			$tmp =  $gene_info{ $i }->{mRNA_lst}->[0];
			print ReadFixGff_GFF $gene_info{ $i }->{contig},"\t",
			"maker","\t",
			"mRNA","\t",
			$mRNA_info{$tmp}->{start},"\t",
			$mRNA_info{$tmp}->{end},"\t",
			$mRNA_info{$tmp}->{score},"\t",
			$mRNA_info{$tmp}->{strand},"\t",
			".","\t",
			"ID=",$pref,"_T",$tmpid,".1;",
			"Parent=",$pref,"_G",$tmpid,";",
			"Name=",$pref,"_T",$tmpid,".1;",
			"locus_tag=",$pref,"_G",$tmpid,"\n";
					
			if ( $mRNA_info{ $tmp }->{exon_num} ) {
				@order = sort { $exon_info{ $a }->{ start } <=> $exon_info{$b}->{start}
				or $exon_info{ $a }->{end} <=> $exon_info{$b}->{end} } @{ $mRNA_info{ $tmp }->{exon_lst} };
				foreach my $j ( 0 .. $#order ) {
					$tmpseq = substr $name2seq{ $exon_info{ $order[$j] }->{contig} }, $exon_info{ $order[$j] }->{start} - 1, $exon_info{ $order[$j] }->{end} - $exon_info{ $order[$j] }->{start} + 1;
					$tmpexon .= $tmpseq;
					print ReadFixGff_GFF $gene_info{$i}->{contig},"\t",
					"maker","\t",
					"exon","\t",
					$exon_info{ $order[$j] }->{start},"\t",
					$exon_info{ $order[$j] }->{end},"\t",
					$exon_info{ $order[$j] }->{score},"\t",
					$exon_info{ $order[$j] }->{strand},"\t",
					".","\t",
					"ID=",$pref, "_T",$tmpid, ".1","-",$j+1,";",
					"Parent=",$pref, "_T",$tmpid,".1",";",
					"locus_tag=",$pref, "_G",$tmpid,";",
					"transcript_id=",$pref, "_T",$tmpid,".1","\n";
				}#
				if ($mRNA_info{$tmp}->{strand} eq "+"){
					#ok, do nothing
				}elsif ($mRNA_info{$tmp}->{strand} eq "-"){
					$tmpexon = rc ($tmpexon);
				}else{
					die "wrong strand.\n";
				}#else
				#print seq
				
			}else{
				print "warning: no exon in gene/mRNA.\n";
			}#else						
			if ( $mRNA_info{ $tmp }->{CDS_num} ) {
				@order = sort { $CDS_info{ $a }->{ start } <=> $CDS_info{$b}->{start}
				or $CDS_info{ $a }->{end} <=> $CDS_info{$b}->{end} } @{ $mRNA_info{ $tmp }->{CDS_lst} };
				foreach my $j ( 0 .. $#order ) {
					$tmpseq = substr $name2seq{ $CDS_info{ $order[$j] }->{contig} }, $CDS_info{ $order[$j] }->{start} - 1, $CDS_info{ $order[$j] }->{end} - $CDS_info{ $order[$j] }->{start} + 1;
					$tmpcds .= $tmpseq;
					print ReadFixGff_GFF $gene_info{$i}->{contig},"\t",
					"maker","\t",
					"CDS","\t",
					$CDS_info{ $order[$j] }->{start},"\t",
					$CDS_info{ $order[$j] }->{end},"\t",
					$CDS_info{ $order[$j] }->{score},"\t",
					$CDS_info{ $order[$j] }->{strand},"\t",
					$CDS_info{ $order[$j] }->{phase},"\t",
					"ID=",$pref,"_C",$tmpid,".1","-",$j+1,";",
					"Parent=",$pref,"_T",$tmpid,".1",";",
					"locus_tag=",$pref, "_G",$tmpid,";",
					"protein_id=",$pref,"_P",$tmpid,".1";
					if(exists $CDS_info{ $order[$j] }->{transl_table} ) {
						print ReadFixGff_GFF ";transl_table=",$CDS_info{ $order[$j] }->{transl_table},"\n";
					}else{
						print ReadFixGff_GFF "\n";
					}#else
				}#

				if ($mRNA_info{$tmp}->{strand} eq "+"){
					if ($CDS_info{ $order[0] }->{phase} != 0 ) {#
						die "non-zero phase for mdl $i.\n";
					}#if
				}elsif ($mRNA_info{$tmp}->{strand} eq "-"){
					if ($CDS_info{ $order[-1] }->{phase} != 0 ) {#
						if ($#order > 0){
							foreach my $j (0 .. $#order){
								print $j,"\t",$CDS_info{ $order[$j] }->{contig},"\t",$CDS_info{$order[$j] }->{start},"\t",$CDS_info{$order[$j]}->{end},
								"\t",$CDS_info{$order[$j]}->{strand},"\t",$CDS_info{$order[$j]}->{phase},"\t",$CDS_info{$order[$j]}->{ID},"\t",$CDS_info{$order[$j]}->{Parent},"\n";
							}#
						}#if
						
						die "non-zero phase for mdl $i.\n";
					}#if
					$tmpcds = rc ($tmpcds);
				}else{
					die "wrong strand.\n";
				}#else
				$tmpcds = uc($tmpcds);
				#generate aa seq
				$codon = substr $tmpcds, 0, 3;
				if ( exists $CDS_info{ $order[0] }->{transl_table} ) {
					if ( $CDS_info{ $order[0] }->{transl_table} == 3 ) {
						if (exists $tt3_start{$codon}){
							$tmpaa = $tt3_start{$codon};
						}elsif (exists $transl_tab3{$codon}){#incomplete
							$tmpaa = $transl_tab3{$codon};
						}else{
							$tmpaa = "X";
							print "warning: unrecognized start codon $codon.\n";
						}#else
					}elsif( $CDS_info{ $order[0] }->{transl_table} == 1 ) {
						if (exists $tt1_start{$codon}){
							$tmpaa = $tt1_start{$codon};
						}elsif(exists $transl_tab1{$codon}){
							$tmpaa = $transl_tab1{$codon};
						}else{
							$tmpaa = "X";
							print "warning: unrecognized start codon $codon.\n";
						}#else
					}else{
						print "warning: unsupport transl_table.\n";
					}#else
				}else{#suppose to be 1
					if (exists $tt1_start{$codon}){
						$tmpaa = $tt1_start{$codon};
					}elsif (exists $transl_tab1{ $codon } ) {
						$tmpaa = $transl_tab1{ $codon };
					}else{
						$tmpaa = "X";
						print "warning: unrecognized start codon $codon.\n";
					}#else
				}#else
				#tranl the middle and stop codon
				foreach my $j (1 .. ( (length $tmpcds) - (length $tmpcds) % 3 ) / 3 - 1){
					$codon = substr $tmpcds, $j*3, 3;
					
					if ( exists $CDS_info{ $order[0] }->{transl_table} ) {
						if ( $CDS_info{ $order[0] }->{transl_table} == 3 ) {
							if (exists $transl_tab3{$codon}){
								$tmpaa .= $transl_tab3{$codon};
							}else{
								$tmpaa .= "X";
								print "warning: unrecognized codon $codon.\n";
							}#else
						}elsif( $CDS_info{ $order[0] }->{transl_table} == 1 ) {
							if (exists $transl_tab1{$codon}){
								$tmpaa .= $transl_tab1{$codon};
							}else{
								$tmpaa .= "X";
								print "warning: unrecognized codon $codon.\n";
							}#else
						}else{
							print "warning: unsupport transl_table.\n";
						}#else
					}else{#suppose to be 1
						if (exists $transl_tab1{$codon}){
							$tmpaa .= $transl_tab1{$codon};
						}else{
							$tmpaa .= "X";
							print "warning: unrecognized codon $codon.\n";
						}#else
					}#else
				}#foreach my $j
			}else{
				print "warning: no CDS in gene/mRNA.\n";
			}#else
		}else{
			print "warning: no mRNA in gene.\n";
		}#else	
		#output seq
		print ReadFixGff_GENE ">", $pref, "_G",$tmpid,"\n",$tmpgene,"\n";
		print ReadFixGff_MRNA ">", $pref, "_T",$tmpid,".1\n",$tmpexon,"\n";
		print ReadFixGff_CDS  ">", $pref, "_C",$tmpid,".1\n",$tmpcds,"\n";
		print ReadFixGff_PROT ">", $pref, "_P",$tmpid,".1\n",$tmpaa,"\n";
	}#foreach my $i				
	close ReadFixGff_GENE;
	close ReadFixGff_MRNA;
	close ReadFixGff_CDS;
	close ReadFixGff_PROT;
	close ReadFixGff_GFF;
}#sub ReadFixGff 


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


