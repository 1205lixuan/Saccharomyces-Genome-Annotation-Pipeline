A Saccharomyes genome annotation pipeline


Introduction:

This pipeline is used to perform large-scale genome annotation 
for Saccharomyes cerevisiae and its relatives. It performs two 
tasks: (1) generating gene models with a two-step annotation 
procedure, and (2) clustering of predicted genes into protein
families.

The Perl scripts used to perform each task are placed in the 
subdir ./scripts. Users can call these scripts by running the 
master Perl script in the current dir ./annt_cluster_pipeline.plx.

Besides, please keep other subdirs ./ref_genome and ./tmp_blastout.
The former includes the reference files (default for strain S288c).
The latter is for the blast result files.


How to run the pipeline:

Inputs:

Before running the pipeline, the genomes sequences to annotate
should be placed in the subdir ./genome_seq. The predicted gene 
models will also be placed in this subdir by the pipeline. 
Meanwhile, two files containing the sequence file information
should be provided in the current dir (./). (1). A file containing 
the file names of these genome sequences and a tag for each genome set, 
which is used for the genome annotation task (a sample file is 
./genome3_test.list). (2). A file containing the file names of protein
sequences and a tag for each set of the protein sequences, which is
used for the gene (protein) family cluster task (a sample file is
./faa3_test.list). 


Command:

$perl ./annt_cluster_pipeline.plx genome3_test.list faa3_test.list

It will takes very long time for thousands of genomes and the 
time consuming steps are MAKER annotation and clustering of 
tens of millions of genes.

Parameters:
	can be specified in the script file ./annt_cluster_pipeline.plx.

	$cpu_num:
		integer (> 0), number of threads to use performing sequence search
		with blast and diamond.
	$cov:
		integer (1 - 100), coverage cutoff to clustering gene families.
	$iden:
		integer (1 - 100), dentity cutoff to clustering gene families.
	$runpref:
		a prefix to tag the current run.
	$fampref:
		a prefix for the output file of the candidate gene families and for
		naming all obtained candidate gene families.
	$fampref2:
		a prefix for the output file of the gene families and for naming all 
		obtained gene families.

Outputs:
	Gene models will be placed in the subdir ./genome_seq.
	Gene/protein families will be placed in the current dir ($fampref2.tab).
	
Other things:
	The pipeline employs the MAKER program to annotate the genomes. Therefore,
	users should install MAKER correctly before running this pipeline, or 
	have to comment out the command (1.2) in the file ./annt_cluster_pipeline.plx
	by placing a '#' before the command (as shown in the current file).
	Meanwhile, users should provide a series of gene model files (with .maker.
	in the file names as shown in the subdir ./genome_seq). That is to tell
	the pipeline to get gene models produced by MAKER from these files.
	A few protein reference files were placed in the subdir ./files_for_maker_annotation.
	Users need to placed the references file and configure files for MAKER program
	in the correct dir based on their own system.
	
