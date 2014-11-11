#!/usr/local/bin/perl

#########
# Developed by Mattia Belli 2013
#########
use warnings;
use strict;
use Statistics::RankCorrelation;
use Statistics::Distributions;
use List::Util qw(shuffle);
#use diagnostics;

print "File codonw: ";
my $file = <STDIN>;
chomp $file;
open (CODONW_FILE, $file) or die "Impossibile aprire $file.";


print "File tRNA codon stats: ";
my $file2 = <STDIN>;
chomp $file2;
open (TRNA_STATS, $file2) or die "Impossibile aprire $file2 .";

print "P-value calculation method (1 T-student (Fast, default), 2 Permutatons (slow), 1 or 2 ?: ";
my $answer = <STDIN>;
chomp $answer  ;
$answer = 1 unless $answer;
$answer = 2 if $answer == 2;
	
	


open (OUT, ">", "$file - CORRELATION.txt") or die "Impossibile aprire.";
open (OUT2, ">", "$file - CORRELATION SUMMARY.txt") or die "Impossibile aprire.";
open (OUT3, ">", "$file - CORRELATION RBCL ATPB.txt") or die "Impossibile aprire.";

if ($answer == 1){
	print OUT2 "NAME	DIVISION	CLASSIFICATION	RHO	T-Value	PVALUE	RHO WOBBLING	T-Value WOBBLING	P-VALUE WOBBLING	RHO superwobbling_AtoI	T-Value superwobbling_AtoI	P-VALUE superwobbling_AtoI\n";
	print OUT3 "NAME	GENE	DIVISION	CLASSIFICATION	RHO	T-Value	PVALUE	RHO WOBBLING	T-Value WOBBLING	P-VALUE WOBBLING	RHO superwobbling_AtoI	T-Value superwobbling_AtoI	P-VALUE superwobbling_AtoI\n";
} elsif ($answer == 2) {
	print OUT2 "NAME	DIVISION	CLASSIFICATION	RHO	P-Perm-VALUE	RHO WOBBLING	P-Perm-VALUE WOBBLING	RHO superwobbling_AtoI	P-Perm-VALUE superwobbling_AtoI\n";
	print OUT3 "NAME	GENE	DIVISION	CLASSIFICATION	RHO	P-Perm-VALUE	RHO WOBBLING	P-Perm-VALUE WOBBLING	RHO superwobbling_AtoI	P-Perm-VALUE superwobbling_AtoI\n";
}


my %codon2AA = (
	'GCA' => 'Ala',
	'GCC' => 'Ala',
	'GCG' => 'Ala',
	'GCT' => 'Ala',
	'AGA' => 'Arg',
	'AGG' => 'Arg',
	'CGA' => 'Arg',
	'CGC' => 'Arg',
	'CGG' => 'Arg',
	'CGT' => 'Arg',
	'AAC' => 'Asn',
	'AAT' => 'Asn',
	'GAC' => 'Asp',
	'GAT' => 'Asp',
	'TGC' => 'Cys',
	'TGT' => 'Cys',
	'CAA' => 'Gln',
	'CAG' => 'Gln',
	'GAA' => 'Glu',
	'GAG' => 'Glu',
	'GGA' => 'Gly',
	'GGC' => 'Gly',
	'GGG' => 'Gly',
	'GGT' => 'Gly',
	'CAC' => 'His',
	'CAT' => 'His',
	'ATA' => 'Ile',
	'ATC' => 'Ile',
	'ATT' => 'Ile',
	'CTA' => 'Leu',
	'CTC' => 'Leu',
	'CTG' => 'Leu',
	'CTT' => 'Leu',
	'TTA' => 'Leu',
	'TTG' => 'Leu',
	'AAA' => 'Lys',
	'AAG' => 'Lys',
	'ATG' => 'Met',
	'TTC' => 'Phe',
	'TTT' => 'Phe',
	'CCA' => 'Pro',
	'CCC' => 'Pro',
	'CCG' => 'Pro',
	'CCT' => 'Pro',
	'AGC' => 'Ser',
	'AGT' => 'Ser',
	'TCA' => 'Ser',
	'TCC' => 'Ser',
	'TCG' => 'Ser',
	'TCT' => 'Ser',
	'TAA' => 'STOP',
	'TAG' => 'STOP',
	'TGA' => 'STOP',
	'ACA' => 'Thr',
	'ACC' => 'Thr',
	'ACG' => 'Thr',
	'ACT' => 'Thr',
	'TGG' => 'Trp',
	'TAC' => 'Tyr',
	'TAT' => 'Tyr',
	'GTA' => 'Val',
	'GTC' => 'Val',
	'GTG' => 'Val',
	'GTT' => 'Val',
);
my %codon_tot;
my %codon_spec;
my %codon_gene;
my %codon_spec_list;
my $codon_n;



while (<CODONW_FILE>) {
	next if $. == 1;
	++$codon_n;
	my @raw = split ("\t", $_);
	my $name = shift(@raw);
	$raw[-1] =~ s/\n//g;
	$codon_tot{$name} = join ("\t",@raw);
	my $species = $name;
	my $gene;
	if ($name =~ m!(.+)\s\[(.+)\]!){
		$species = $1;
		$gene = $2;
		$codon_gene{$name} = $gene;
	}
	$codon_spec{$name} = $species;
	$codon_spec_list{$species} = 1;
}

my %trna_stats;
my %classification;
my %tax_data;

while (<TRNA_STATS>) {
	next if $. == 1;
	my @raw = split ("\t", $_);
	my $species = shift(@raw);
	# next unless exists $codon_spec_list{$species};
	$tax_data{$species} = shift(@raw);
	shift(@raw);
	$classification{$species} = shift(@raw);
	shift(@raw);
	shift(@raw);
	$raw[-1] =~ s/\n//g;
	$trna_stats{$species} = join ("\t",@raw);
	# print "\n",join ("\t",@raw),"\n";
	# <>;
}

close TRNA_STATS;
close CODONW_FILE;

my $proc_codon;
foreach (sort keys %codon_tot) {
	++$proc_codon;
	my $name = $_;
	my $species = $codon_spec{$name};
	
	next unless exists 	$trna_stats{$species};	
	
	my %trna_stats_data;
	my %codon_tot_data;
	my %codon_tot_data_wobbling;
	my %codon_tot_data_superwobbling_AtoI;
	
	print "Processing $proc_codon..$name\n";
	
	while ($trna_stats{$codon_spec{$name}} =~ m!([ATCG]{3})\s([\d\.]+)!g){
		# print "$_ $1 $2\n";
		# <>;
		$trna_stats_data{$1} = $2;
	}
	while ($codon_tot{$name} =~ m!([ATCG]{3})\s([\d\.]+)!g){
		next unless exists $trna_stats_data{$1};
		my $triplet = $1;
		my $count = $2;
		$codon_tot_data{$triplet} += $count;
		$codon_tot_data_wobbling{$triplet} += $count;
		$codon_tot_data_superwobbling_AtoI{$triplet} += $count;
		my @triplet = split ("", $triplet);
		my $triplet_wobbled;


		
		##Crick wobbling rules implementation
		#Anticodons with G in first position (GNN) can pair with NNC (NNG<->NNC) and NNT (NNG<->NNT) codons
		#Anticodons with U in first position (UNN) can pair with NNA (NNU<->NCA) and NNG (NNU<->NNG) but stop codon TGA codon is not recognized by tRNA-Trp
		#Exceptions: Cytidine to Lysidine conversion tRNA-Ile with anticodon CAT (CAU) can pair with AUA codons but not with AUG codons (for Methionine)
		#Exceptions: Adenosine to Inosine conversion tRNA-Arg (ACG) can pair with CGT, CGC, CGA.

		if ($triplet[2] eq 'T') {
			$triplet_wobbled = $triplet[0].$triplet[1].'C';	#GXX anticodons read XXU and XXC codons
			$codon_tot_data_wobbling{$triplet_wobbled} += $count;
			$codon_tot_data_superwobbling_AtoI{$triplet_wobbled} += $count;
		}
		if (($triplet[2] eq 'G') && ($triplet ne 'ATG') && ($triplet ne 'TGG')) { #
			$triplet_wobbled = $triplet[0].$triplet[1].'A';	#UXX anticodons (where U is modified) read XXA and XXG codons except for stop codon UGA
			$codon_tot_data_wobbling{$triplet_wobbled} += $count;
			$codon_tot_data_superwobbling_AtoI{$triplet_wobbled} += $count;			
		}
		if (($triplet eq 'CGC')|($triplet eq 'CGA')) { #Arginine anticodons tRNA-Arg (ACG) can pair with CGT, CGC, CGA because Adenosine is modified in Inosine
			$triplet_wobbled = 'CGT';
			$codon_tot_data_superwobbling_AtoI{$triplet_wobbled} += $count;			
		}
		# print "$triplet	$codon_tot_data{$triplet}\n";
		# print "$triplet_wobbled	$codon_tot_data_wobbling{$triplet}\n";
		# <>;
	}		

	my @trna_data;
	my @codon_data;
	my @trna_count;
	my @codon_count;
	my @codon_data_wobbling;
	my @codon_data_superwobbling_AtoI;
	foreach (sort keys %trna_stats_data) {
		push (@trna_count, "$_ $codon2AA{$_} $trna_stats_data{$_}");
		push (@codon_count, "$_ $codon2AA{$_} $codon_tot_data{$_}");
		push (@trna_data, $trna_stats_data{$_});
		push (@codon_data, $codon_tot_data{$_});
		push (@codon_data_wobbling, $codon_tot_data_wobbling{$_});
		push (@codon_data_superwobbling_AtoI, $codon_tot_data_superwobbling_AtoI{$_});
	}
	# print "\n$name\n@trna_data\n";
	# print "\n@codon_data\n";
	# print "\n@codon_data_wobbling\n";
	# <>;
	my $n1 = scalar @trna_data;
	my $n2 = scalar @codon_data;
	die "Number of trna elements is not equal to the number of codon elements!\n" unless $n1 == $n2;
	
	my $spearman;
	my $t_value;
	my $pvalue;
	my $spearman_wobbling;
	my $t_value_wobbling;
	my $pvalue_wobbling;
	my $spearman_superwobbling_AtoI;
	my $t_value_superwobbling_AtoI;
	my $pvalue_superwobbling_AtoI;
	
	if ($answer == 1){
		my $c = Statistics::RankCorrelation->new( \@trna_data, \@codon_data, sorted => 1 );
		$spearman = $c->spearman;
		# print "s $spearman\n";
		# <>;
		$t_value = $spearman*sqrt(($n1-2)/(1-$spearman**2));
		$pvalue = Statistics::Distributions::tprob(($n1-2), $t_value); #Upper probability, one-tail"
		#$pvalue = Statistics::Distributions::tprob(($n1-2), abs $t_value); #Upper probability, one-tail"
		#$pvalue *= 2; #Two-tail pvalue
		
		my $c_wobbling = Statistics::RankCorrelation->new( \@trna_data, \@codon_data_wobbling, sorted => 1 );
		$spearman_wobbling = $c_wobbling->spearman;
		$t_value_wobbling = $spearman_wobbling*sqrt(($n1-2)/(1-$spearman_wobbling**2));
		$pvalue_wobbling = Statistics::Distributions::tprob(($n1-2), $t_value_wobbling); #Upper probability, one-tail"
		#$pvalue_wobbling = Statistics::Distributions::tprob(($n1-2), abs $t_value_wobbling); #Upper probability, one-tail"
		#$pvalue_wobbling *= 2; #Two-tail pvalue
		
		my $c_superwobbling_AtoI = Statistics::RankCorrelation->new( \@trna_data, \@codon_data_superwobbling_AtoI, sorted => 1 );
		$spearman_superwobbling_AtoI = $c_superwobbling_AtoI->spearman;
		$t_value_superwobbling_AtoI = $spearman_superwobbling_AtoI*sqrt(($n1-2)/(1-$spearman_superwobbling_AtoI**2));
		$pvalue_superwobbling_AtoI = Statistics::Distributions::tprob(($n1-2), $t_value_superwobbling_AtoI); #Upper probability, one-tail"
		#$pvalue_superwobbling_AtoI = Statistics::Distributions::tprob(($n1-2), abs $t_value_superwobbling_AtoI); #Upper probability, one-tail"
		#$pvalue_superwobbling_AtoI *= 2; #Two-tail pvalue
		
	} elsif ($answer == 2) {
		($spearman, $pvalue) = &spearman_test(\@trna_data, \@codon_data);
		($spearman_wobbling, $pvalue_wobbling) = &spearman_test(\@trna_data, \@codon_data_wobbling);	
		($spearman_superwobbling_AtoI, $pvalue_superwobbling_AtoI) = &spearman_test(\@trna_data, \@codon_data_superwobbling_AtoI);
	}
	
	my $division = 'NO DIVISION DATA';
	$division = $tax_data{$species} if exists $tax_data{$species};
	
	# print "NAME |$name|\n";
	# print "SPECIES |$species|\n";
	# print "CLASSIFICATION $classification{$species}\n";
	# print "Spearman $spearman\n";
	# print "T-value $t_value\n";
	# print "Division $division\n";
	# <>;
	if ($answer == 1){
		print OUT "$name\t$division\t$classification{$species}\n";
		print OUT "\ttRNA count\t".join ("\t",@trna_count)."\n\tCodon count\t".join ("\t",@codon_count)."\n\ttRNA data\t".join ("\t",@trna_data)."\n";
		print OUT "\tCodon data\t".join ("\t",@codon_data)."\trho:$spearman\tt:$t_value\tpv:$pvalue\n";
		print OUT "\tCodon data Wobbling\t".join ("\t",@codon_data_wobbling)."\trhow:$spearman_wobbling\ttw:$t_value_wobbling\tpvw:$pvalue_wobbling\n";
		print OUT "\tCodon data superwobbling_AtoI\t".join ("\t",@codon_data_superwobbling_AtoI)."\trhow:$spearman_superwobbling_AtoI\ttw:$t_value_superwobbling_AtoI\tpvw:$pvalue_superwobbling_AtoI\n";
		print OUT2 "$name\t$division\t$classification{$species}\t$spearman\t$t_value\t$pvalue\t$spearman_wobbling\t$t_value_wobbling\t$pvalue_wobbling\t$spearman_superwobbling_AtoI\t$t_value_superwobbling_AtoI\t$pvalue_superwobbling_AtoI\n";
		print OUT3 "$species\t$codon_gene{$name}\t$division\t$classification{$species}\t$spearman\t$t_value\t$pvalue\t$spearman_wobbling\t$t_value_wobbling\t$pvalue_wobbling\t$spearman_superwobbling_AtoI\t$t_value_superwobbling_AtoI\t$pvalue_superwobbling_AtoI\n" if (exists $codon_gene{$name} && $codon_gene{$name} =~ m/(rbcl|atpb)/i);
	} elsif ($answer == 2) {
		print OUT "$name\t$division\t$classification{$species}\n";
		print OUT "\ttRNA count\t".join ("\t",@trna_count)."\n\tCodon count\t".join ("\t",@codon_count)."\n\ttRNA data\t".join ("\t",@trna_data)."\n";
		print OUT "\tCodon data\t".join ("\t",@codon_data)."\trho:$spearman\t\tpv:$pvalue\n";
		print OUT "\tCodon data Wobbling\t".join ("\t",@codon_data_wobbling)."\trhow:$spearman_wobbling\tpvw:$pvalue_wobbling\n";
		print OUT "\tCodon data superwobbling_AtoI\t".join ("\t",@codon_data_superwobbling_AtoI)."\trhow:$spearman_superwobbling_AtoI\tpvw:$pvalue_superwobbling_AtoI\n";		
		print OUT2 "$name\t$division\t$classification{$species}\t$spearman\t$pvalue\t$spearman_wobbling\t$pvalue_wobbling\t$spearman_superwobbling_AtoI\t$pvalue_superwobbling_AtoI\n";
		print OUT3 "$species\t$codon_gene{$name}\t$division\t$classification{$species}\t$spearman\t$pvalue\t$spearman_wobbling\t$pvalue_wobbling\t$spearman_superwobbling_AtoI\t$pvalue_superwobbling_AtoI\n" if (exists $codon_gene{$name} && $codon_gene{$name} =~ m/(rbcl|atpb)/i);
	}

}

close OUT;
close OUT2;
close OUT3;


print "\n\n****FATTO****\n\n";

sub spearman_test {

	my ($trna_ref, $cod_ref) = @_;

	my $c = Statistics::RankCorrelation->new( $trna_ref, $cod_ref, sorted => 1 );
	my $observed_r = $c->spearman;

	my $maxnumber=999; #9999 Random Rho + the true rho = 1000

	my $i=0;
	my $sum = 0;

	while ($i<$maxnumber) {
		my @perm = shuffle @$trna_ref;
		# print "$i\.	@perm\n";
		# <>;
		my $random_c = Statistics::RankCorrelation->new( \@perm, $cod_ref, sorted => 1 );
		my $rc = $random_c->spearman;
		$sum += 1 if ($rc>= $observed_r); #One tail, upper
		#$sum += 1 if (abs $rc>= abs $observed_r); #Two tail
		# push(@set, "@perm",);
		++$i;
	}	
	my $pvalue = $sum/($maxnumber+1);
	# print "$observed_r, $pvalue\n";
	return $observed_r, $pvalue;
}
  



	# my $R = Statistics::R->new();
	# $R->startR;
	# my $x = join(',',@trna_data);
	# my $y = join(',',@codon_data);
	
	# $R->send(qq`x<-c($x)`);
	# $R->send(qq`y<-c($y)`);
	
	# $R->send(q`library(coin)`);
	# $R->send(q`pvalue(spearman_test(x~y, distribution=approximate(B=999)))`);
	# my $pvalue = $R->read;
	
	# $R->stopR();

print "\nFATTO";
<>;