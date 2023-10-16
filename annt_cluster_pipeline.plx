#!/usr/bin/perl
#annt_cluster_pipeline.plx
use strict;
use warnings;

# set the number of cpu threads to use by blastn/p and diamond blastp:
my $cpu_num  = 10;

# predefined directories and reference file
# Please do NOT change these parameters:
my $script_dir = "./scripts";
my $ref_dir = "ref_genome";
my $genome_dir = "genome_seq";
my $blastout_dir = "tmp_blastout";

my $ref_prot_faa = "GCF_000146045.2_R64_protein.shortname.faa";



# get parameter from command line,
my $infile  = shift;#fna_list file in the working directory. including the genome sequence to annotate. all fna files should be placed in the genome_dir
my $infile2 = shift;#faa_list file in the working directory. including both the annotated genomes and reference S288C. all faa files shoudl be placed in the genome_dir
if (not $infile2){
	die "usage: perl annt_cluster_pipeline.plx <query_fa_list> <query_faa_list>\n";
}#if
my $outfile = $infile . ".log";


##############################
#  1. Batch genome annotation
##############################

# 1.1. run the blastn annt pipeline
system ("perl $script_dir/blastn_annt_bat.plx $cpu_num $ref_dir $genome_dir $infile");

# 1.2. run the maker annt pipeline
#system ("perl $script_dir/maker_annt_pipeline3.2.plx $cpu_num $genome_dir $infile");

# 1.3. merge annt, based on overlap on the chr (also exclude suspicious models, e.g. too short)
system ("perl $script_dir/merge_annt_bat.plx $genome_dir $infile");

##############################
#  2. Gene family clustering
##############################

# set the parameters for family clustering and output file name
my $cov = 50;
my $iden = 50;
my $runpref = "rundiamondc50i50test5";
my $fampref = "DSCF50TEST5";
my $fampref2= "DSCF50TEST5S";

# 2.1. get blastn homolog info
system ("perl $script_dir/blastn_annt_2_family.plx $ref_dir $genome_dir $infile $infile.blastn_annt_homolog");
# 2.2. cluster the ref proteome with blastp
system ("perl $script_dir/blast_cluster3.1.plx $cpu_num $cov $iden prot $ref_dir/$ref_prot_faa");
# 2.3. match the maker models to the reference models
system ("perl $script_dir/maker_ref_blastp_bat2.plx $cpu_num $ref_dir $ref_prot_faa $genome_dir $infile $blastout_dir");
# 2.4. cluster the maker models not matched to the reference models, use diamond blastp to speed up
system ("perl $script_dir/maker_clust2_diamond_2.plx $cpu_num $cov $iden $ref_dir $ref_prot_faa $genome_dir $infile $runpref");
# 2.5. create candidate families based on the information obtained above
system ("perl $script_dir/create_family.plx $ref_dir/$ref_prot_faa.clust.$cov.$iden.info $infile.blastn_annt_homolog $runpref.detail.out $runpref.maker_only_cluster.out $fampref");
# 2.6. create final families
system ("perl $script_dir/fam_super_cluster.plx 1 $fampref.tab $infile2 $genome_dir $iden $cov $fampref2");
