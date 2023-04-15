#!/usr/bin/perl
# create_family.plx
use strict;
use warnings;

$|=1;
my $ref_cluster_file= shift;
my $blastn_annt_homolog = shift;
my $ref_map_file    = shift;
my $maker_only_cluster_file = shift;
my $out_prefix      = shift;
if (not $out_prefix){
	die "usage: perl create_family.plx <ref_cluster_file> <blastn_annt_homolog tab> <ref_maker_map_file.blastp> <maker_only_cluster_file> <out_pref>\n";
}#if
my $outfile = $out_prefix . ".tab";

my $fam_num = 0;

my @items;
my @items2;
my $lbuf;
my @strain_name;
my %strname2idx;
my %ref_annt;	#{$ref_id}->[$str_idx]->[idx] => seq in strain $str_idx	
my %ref_maker;	#{ref}->[ $str_idx ]->[$idx] =>	seq in strain $str_idx	
my @tmparr;
my $tmpline;
open IN, $blastn_annt_homolog or die "fail to open $blastn_annt_homolog:$!\n";
$lbuf = <IN>;

chomp $lbuf;

@strain_name = split /\t/,$lbuf;

foreach my $i (0 .. $#strain_name){
	$strname2idx{ $strain_name[ $i ] } = $i;
}#
print "strain_name size: ",scalar @strain_name,"\n";
foreach my $i (sort keys %strname2idx){
	print $i,"\t",$strname2idx{$i},"\n";
}#
while (<IN>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/, $_;
	foreach my $i (0 .. $#items){
		next if not $items[$i];
		@items2 = split /\|/,$items[$i];
		foreach my $j (0 .. $#items2){
			push @{ $ref_annt{ $items[0] }->[ $i ] }, $items2[$j];
		}#
	}#
}#
close IN;

open IN, $ref_map_file or die "fail to open $ref_map_file:$!\n";
while (<IN>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/, $_;
	@items2 = split /_/, $items[0];#get the strain name
	#print "@items2\n";
	push @{ $ref_maker{ $items[1] }->[ $strname2idx{ $items2[0] } ] }, $items[0];
}#
close IN;

open OUT, ">$outfile" or die "fail to open $outfile:$!\n";
print OUT "id","\t",$strain_name[0];
foreach my $i (1 .. $#strain_name){
	print OUT "\t",$strain_name[$i];
}
print OUT "\n";
open IN, $ref_cluster_file or die "fail to open $ref_cluster_file:$!\n";
while (<IN>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/, $_;
	printf OUT "%s%s%05d\t",$out_prefix,"_",++$fam_num;
	print OUT $items[2];
	foreach my $i (3 .. $#items){
		print OUT "|",$items[$i];
	}
	foreach my $i (1 .. $#strain_name){#each strain
		@tmparr = ();
		foreach my $j (2 .. $#items){#each ref prot id
			if (exists $ref_annt{ $items[$j] }->[ $i ] ) {
				push @tmparr, @{ $ref_annt{ $items[$j] }->[ $i ] };
			}
			if (exists $ref_maker{ $items[$j] }->[ $i ] ) {
				push @tmparr, @{ $ref_maker{$items[$j] }->[ $i ] };
			}#
		}#
		if (scalar @tmparr == 0 ){
			print OUT "\t";
		}else{
			$tmpline = join "|",@tmparr;
			print OUT "\t",$tmpline;
		}#else
	}#foreach my $i
	print OUT "\n";
}#
close IN;


open IN, $maker_only_cluster_file or die "fail to open $maker_only_cluster_file:$!\n";
while (<IN>){
	chomp;
	next if not $_;
	next if /^#/;
	@items = split /\t/, $_;
	@tmparr = ();
	foreach my $i (0 .. $#strain_name){
		$tmparr[$i] = "";
	}#
	foreach my $i (1 .. $#items){#
		if ($items[$i]){
			@items2 = split /_/,$items[$i];#get the strain name
			$tmparr[ $strname2idx{ $items2[0] } ] = $items[$i];
		}#if
	}#
	$tmpline = join "\t",@tmparr;
	printf OUT "%s%s%05d\t",$out_prefix,"_",++$fam_num;
	print OUT $tmpline,"\n";
}#while
close IN;
			
close OUT;			

