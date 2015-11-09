#!/bin/perl

#########
# Developed by Mattia Belli 2013
#########
use warnings;
use strict;
use Statistics::RankCorrelation;
use Statistics::Distributions;
use List::Util qw(shuffle max);
use Statistics::R 0.33;
use Benchmark;
#use Array::Shuffle qw(shuffle_array);
#use diagnostics;

print "File codonw: ";
my $file = <STDIN>;
chomp $file;
open (CODONW_FILE, $file) or die "Impossibile aprire $file.";


print "File tRNA codon stats: ";
my $file2 = <STDIN>;
chomp $file2;
open (TRNA_STATS, $file2) or die "Impossibile aprire $file2 .";

print "p-value calculation method:

		(1) T-student approximation (Fast, default)
		(2) R Monte Carlo resampling (slower but more reliable)
		(3) Perl Permutation (very slow)
		
	   Choose one method (1,2,3): ";
my $answer = <STDIN>;
chomp $answer  ;
$answer = 1 unless ($answer == 2 or $answer == 3);

my $start_time = new Benchmark;

if ($answer == 1){
	open (OUT, ">", "$file - CORRELATION.txt") or die "Impossibile aprire.";
	open (OUT2, ">", "$file - CORRELATION SUMMARY.txt") or die "Impossibile aprire.";
	open (OUT3, ">", "$file - tRNA coverage.txt") or die "Impossibile aprire.";
	print OUT2 "NAME	GENE	DIVISION	CLASSIFICATION	RHO	T-Value	PVALUE	RHO WOBBLING	T-Value WOBBLING	P-VALUE WOBBLING	RHO editing_AtoI_CtoL	T-Value editing_AtoI_CtoL	P-VALUE editing_AtoI_CtoL	RHO SUPERWOBBLING	T-Value SUPERWOBBLING	P-VALUE SUPERWOBBLING\n";
} elsif ($answer == 2) {
	open (OUT, ">", "$file - CORRELATION R PERM.txt") or die "Impossibile aprire.";
	open (OUT2, ">", "$file - CORRELATION R PERM SUMMARY.txt") or die "Impossibile aprire.";
	open (OUT3, ">", "$file - tRNA coverage.txt") or die "Impossibile aprire.";
	print OUT2 "NAME	GENE	DIVISION	CLASSIFICATION	RHO	P-Perm-VALUE	RHO WOBBLING	P-Perm-VALUE WOBBLING	RHO editing_AtoI_CtoL	P-Perm-VALUE editing_AtoI_CtoL	RHO SUPERWOBBLING	P-Perm-VALUE SUPERWOBBLING\n";
} elsif ($answer == 3) {
	open (OUT, ">", "$file - CORRELATION Perl PERM.txt") or die "Impossibile aprire.";
	open (OUT2, ">", "$file - CORRELATION Perl PERM SUMMARY.txt") or die "Impossibile aprire.";
	open (OUT3, ">", "$file - tRNA coverage.txt") or die "Impossibile aprire.";
	print OUT2 "NAME	GENE	DIVISION	CLASSIFICATION	RHO	P-Perm-VALUE	RHO WOBBLING	P-Perm-VALUE WOBBLING	RHO editing_AtoI_CtoL	P-Perm-VALUE editing_AtoI_CtoL	RHO SUPERWOBBLING	P-Perm-VALUE SUPERWOBBLING\n";
}

my $R = Statistics::R->new() if ($answer == 2);
$R->startR if ($answer == 2);
$R->run(q`library(coin)`) if ($answer == 2);

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


my %AA2codons = (			# Codons from NCBI The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
   'Ala' => [qw/GCA GCC GCG GCT/],
   'Gly' => [qw/GGA GGC GGG GGT/],
   'Pro' => [qw/CCA CCC CCG CCT/],
   'Thr' => [qw/ACA ACC ACG ACT/],
   'Val' => [qw/GTA GTC GTG GTT/],
   
   'Ser' => [qw/AGC AGT TCA TCC TCG TCT/],
   'Arg' => [qw/AGA AGG CGA CGC CGG CGT/],
   'Leu' => [qw/CTA CTC CTG CTT TTA TTG/],
   
   'Phe' => [qw/TTC TTT/],
   
   'Asn' => [qw/AAC AAT/],
   'Lys' => [qw/AAA AAG/],
   
   'Asp' => [qw/GAC GAT/],
   'Glu' => [qw/GAA GAG/],
   
   'His' => [qw/CAC CAT/],
   'Gln' => [qw/CAA CAG/],
   
   'Tyr' => [qw/TAC TAT/],
   # 'STOP' => [qw/TAA TAG TGA/],
   
   'Ile' => [qw/ATA ATC ATT/],
   'Met' => [qw/ATG/],
   
   'Cys' => [qw/TGC TGT/],
   'Trp' => [qw/TGG/],
   # 'SelCys' => [qw/TGA/]
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
	$codon_gene{$name} = '';
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

	my $name = $_;
	my $species = $codon_spec{$name};
	
	next unless exists 	$trna_stats{$species};	
	++$proc_codon;
	
	my %trna_stats_data;
	my %trna_stats_data_wob;
	my %trna_stats_data_superwob;
	my %codon_tot_data;
	my %codon_tot_data_wobbling;
	my %codon_tot_data_editing_AtoI_CtoL;
	my %codon_tot_data_superwobbling;
	
	foreach (keys %codon2AA){
		$trna_stats_data{$_} = 0;
		$trna_stats_data_wob{$_} = 0;
		$trna_stats_data_superwob{$_} = 0;
	}
	
	
	print "Processing $proc_codon..$name\n";
	
	while ($trna_stats{$codon_spec{$name}} =~ m!([ATCG]{3})\s([\d\.]+)!g){
		my $triplet = $1;
		$trna_stats_data{$triplet} = $2;
		$trna_stats_data_wob{$triplet} = 0;
		$trna_stats_data_superwob{$triplet} = 0;
		# print "$triplet = $trna_stats_data_superwob{$triplet}\n";
	}
	while ($codon_tot{$name} =~ m!([ATCG]{3})\s([\d\.]+)!g){
		next unless exists $trna_stats_data{$1};
		my $triplet = $1;
		my $count = $2;
		$codon_tot_data{$triplet} += $count;
		$codon_tot_data_wobbling{$triplet} += $count;
		$codon_tot_data_editing_AtoI_CtoL{$triplet} += $count;
		$codon_tot_data_superwobbling{$triplet} += $count;
		my @triplet = split ("", $triplet);
		my $triplet_wobbled;
		
		# print "triplet = $triplet\n";
		# <>;

		
		##Crick's wobbling rules implementation
		#Anticodons with G in first position (GNN) can pair with NNC (NNG<->NNC) and NNT (NNG<->NNT) codons
		#Anticodons with U in first position (UNN) can pair with NNA (NNU<->NCA) and NNG (NNU<->NNG) but stop codon TGA codon is not recognized by tRNA-Trp
		#Exception: Cytidine to Lysidine conversion tRNA-Ile with anticodon CAT->LAT can pair with AUA codons but not with AUG codons (for Methionine)
		#Exceptions: Adenosine to Inosine conversion tRNA-Arg (ACG) can pair with CGT, CGC, CGA.

		if ($triplet[2] eq 'T') {
			$triplet_wobbled = $triplet[0].$triplet[1].'C';	#GXX anticodons read XXU and XXC codons
			$codon_tot_data_wobbling{$triplet_wobbled} += $count;
			$codon_tot_data_editing_AtoI_CtoL{$triplet_wobbled} += $count;
			$codon_tot_data_superwobbling{$triplet_wobbled} += $count;
			$trna_stats_data_wob{$triplet} = $trna_stats_data{$triplet_wobbled};
			#$trna_stats_data_superwob{$triplet} = $trna_stats_data{$triplet_wobbled};
		}
		if (($triplet[2] eq 'G') && ($triplet ne 'ATG') && ($triplet ne 'TGG')) { #
			$triplet_wobbled = $triplet[0].$triplet[1].'A';	#UXX anticodons (where U is modified) read XXA and XXG codons except for stop codon UGA
			$codon_tot_data_wobbling{$triplet_wobbled} += $count;
			$codon_tot_data_editing_AtoI_CtoL{$triplet_wobbled} += $count;	
			$codon_tot_data_superwobbling{$triplet_wobbled} += $count if $triplet =~ m!^(TT|TA|TG|CA|AT|AA|AG|GA)!;
			$trna_stats_data_wob{$triplet} = $trna_stats_data{$triplet_wobbled};	
		}
		if (($triplet eq 'CGC')|($triplet eq 'CGA')) { #Arginine anticodons tRNA-Arg (ACG) can pair with CGT, CGC, CGA because Adenosine is modified in Inosine
			$triplet_wobbled = 'CGT';
			$codon_tot_data_editing_AtoI_CtoL{$triplet_wobbled} += $count;	
			$codon_tot_data_superwobbling{$triplet_wobbled} += $count;
			$trna_stats_data_wob{$triplet} = $trna_stats_data{$triplet_wobbled};
		#	$trna_stats_data_superwob{$triplet} = $trna_stats_data{$triplet_wobbled};
		}
		if ($triplet eq 'ATA') { #tRNA-Ile with anticodon CAT->LAT can pair with AUA codons but not with AUG codons (for Methionine)
			$triplet_wobbled = 'ATG';
			$codon_tot_data_editing_AtoI_CtoL{$triplet_wobbled} += $count;	
			$codon_tot_data_superwobbling{$triplet_wobbled} += $count;
			$trna_stats_data_wob{$triplet} = $trna_stats_data{$triplet_wobbled};
		#	$trna_stats_data_superwob{$triplet} = $trna_stats_data{$triplet_wobbled};
		}
		# print "$triplet	$codon_tot_data{$triplet}\n";
		# print "$triplet_wobbled	$codon_tot_data_wobbling{$triplet}\n";
		# <>;
		
		
		#Superwobble rules, unmodified Uracil recognize XXA but also XXC, XXG, XXU
		
		if (($triplet[2] ne 'A') and ($triplet !~ m!^(TT|TA|TG|CA|AT|AA|AG|GA)\w!)){ #only four codons AA
			$triplet_wobbled = $triplet[0].$triplet[1].'A';
			$codon_tot_data_superwobbling{$triplet_wobbled} += $count;
			$trna_stats_data_superwob{$triplet} += $trna_stats_data{$triplet_wobbled};
		}

		
	}		

	my @trna_data;
	my @trna_data_wob;
	my @trna_data_superwob;
	my @codon_data;
	my @codon_data_scaled;
	my @trna_count;
	my @codon_count;
	my @codon_data_wobbling;
	my @codon_data_editing_AtoI_CtoL;
	my @codon_data_superwobbling;
	my @AAlist;
	my @AAlist2;	
	# foreach (sort keys %trna_stats_data) {
		# push (@trna_count, "$_ $codon2AA{$_} $trna_stats_data{$_}");
		# push (@codon_count, "$_ $codon2AA{$_} $codon_tot_data{$_}");
		# push (@trna_data, $trna_stats_data{$_});
		# push (@trna_data_wob, $trna_stats_data_wob{$_});
		# push (@codon_data, $codon_tot_data{$_});
		# push (@codon_data_wobbling, $codon_tot_data_wobbling{$_});
		# push (@codon_data_editing_AtoI_CtoL, $codon_tot_data_editing_AtoI_CtoL{$_});
	# }
	
	foreach (sort keys %AA2codons){
		my $AA = $_;
		my $ref = $AA2codons{$_};
		foreach (sort @$ref){		
			my $anticodon = &reverse_complement($_);
			push (@AAlist, "tRNA($anticodon)-$AA $_");
			push (@AAlist2, "tRNA($anticodon)-$AA($_)	tRNA($anticodon)-$AA($_) wob	tRNA($anticodon)-$AA($_) superwob	tRNA($anticodon)-$AA($_) Codons	");
			push (@trna_count, "$_ $AA $trna_stats_data{$_}");
			push (@codon_count, "$_ $AA $codon_tot_data{$_}");
			push (@trna_data, $trna_stats_data{$_});
			push (@trna_data_wob, $trna_stats_data_wob{$_});
			push (@trna_data_superwob, $trna_stats_data_superwob{$_});
			push (@codon_data, $codon_tot_data{$_});
			push (@codon_data_wobbling, $codon_tot_data_wobbling{$_});
			push (@codon_data_editing_AtoI_CtoL, $codon_tot_data_editing_AtoI_CtoL{$_});
			push (@codon_data_superwobbling, $codon_tot_data_superwobbling{$_});
		}
	}

	# print "\n$name\n@trna_data\n";
	# print "\n@codon_data\n";
	# print "\n@codon_data_wobbling\n";
	# <>;
	my $n1 = scalar @trna_data;
	my $n2 = scalar @codon_data;
	die "Number of tRNAs is not equal to the number of codons!\n" unless $n1 == $n2;
	
	my $spearman;
	my $t_value;
	my $pvalue;
	my $spearman_wobbling;
	my $t_value_wobbling;
	my $pvalue_wobbling;
	my $spearman_editing_AtoI_CtoL;
	my $t_value_editing_AtoI_CtoL;
	my $pvalue_editing_AtoI_CtoL;
	my $spearman_superwobbling;
	my $t_value_superwobbling;
	my $pvalue_superwobbling;
	
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
		
		my $c_editing_AtoI_CtoL = Statistics::RankCorrelation->new( \@trna_data, \@codon_data_editing_AtoI_CtoL, sorted => 1 );
		$spearman_editing_AtoI_CtoL = $c_editing_AtoI_CtoL->spearman;
		$t_value_editing_AtoI_CtoL = $spearman_editing_AtoI_CtoL*sqrt(($n1-2)/(1-$spearman_editing_AtoI_CtoL**2));
		$pvalue_editing_AtoI_CtoL = Statistics::Distributions::tprob(($n1-2), $t_value_editing_AtoI_CtoL); #Upper probability, one-tail"
		#$pvalue_editing_AtoI_CtoL = Statistics::Distributions::tprob(($n1-2), abs $t_value_editing_AtoI_CtoL); #Upper probability, one-tail"
		#$pvalue_editing_AtoI_CtoL *= 2; #Two-tail pvalue		
		
		my $c_superwobbling = Statistics::RankCorrelation->new( \@trna_data, \@codon_data_superwobbling, sorted => 1 );
		$spearman_superwobbling = $c_superwobbling->spearman;
		$t_value_superwobbling = $spearman_superwobbling*sqrt(($n1-2)/(1-$spearman_superwobbling**2));
		$pvalue_superwobbling = Statistics::Distributions::tprob(($n1-2), $t_value_superwobbling); #Upper probability, one-tail"
		#$pvalue_superwobbling = Statistics::Distributions::tprob(($n1-2), abs $t_value_superwobbling); #Upper probability, one-tail"
		#$pvalue_superwobbling *= 2; #Two-tail pvalue
		
	} elsif ($answer == 2) {
		($spearman, $pvalue) = &spearman_test_R(\@trna_data, \@codon_data, $R);
		($spearman_wobbling, $pvalue_wobbling) = &spearman_test_R(\@trna_data, \@codon_data_wobbling, $R);	
		($spearman_editing_AtoI_CtoL, $pvalue_editing_AtoI_CtoL) = &spearman_test_R(\@trna_data, \@codon_data_editing_AtoI_CtoL, $R);
		($spearman_superwobbling, $pvalue_superwobbling) = &spearman_test_R(\@trna_data, \@codon_data_superwobbling, $R);	
	} elsif ($answer == 3) {
		my $start_time = new Benchmark;
		($spearman, $pvalue) = &spearman_test(\@trna_data, \@codon_data);
		($spearman_wobbling, $pvalue_wobbling) = &spearman_test(\@trna_data, \@codon_data_wobbling);	
		($spearman_editing_AtoI_CtoL, $pvalue_editing_AtoI_CtoL) = &spearman_test(\@trna_data, \@codon_data_editing_AtoI_CtoL);
		($spearman_superwobbling, $pvalue_superwobbling) = &spearman_test(\@trna_data, \@codon_data_superwobbling);	
		my $end_time   = new Benchmark;
		my $difference = timediff($end_time, $start_time);
		print "Benchmark: it took ", timestr($difference), "\n";
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
	
	my $max = max @codon_data;
	foreach (@codon_data){
		my $scaled_value = $_/$max;
		push(@codon_data_scaled, $scaled_value);
	}
	
	if ($answer == 1){
		print OUT "$name\t$codon_gene{$name}\t$division\t$classification{$species}\n";
		print OUT "\tAA codon\t".join ("\t",@AAlist)."\n";
		print OUT "\ttRNA count\t".join ("\t",@trna_count)."\n\tCodon count\t".join ("\t",@codon_count)."\n\ttRNA data\t".join ("\t",@trna_data)."\n\ttRNA data_wob\t".join ("\t",@trna_data_wob)."\n";
		print OUT "\tCodon data SCALED\t".join ("\t",@codon_data_scaled)."\n";		
		print OUT "\tCodon data\t".join ("\t",@codon_data)."\trho:$spearman\tt:$t_value\tpv:$pvalue\n";
		print OUT "\tCodon data Wobbling\t".join ("\t",@codon_data_wobbling)."\trhow:$spearman_wobbling\ttw:$t_value_wobbling\tpvw:$pvalue_wobbling\n";
		print OUT "\tCodon data editing_AtoI_CtoL\t".join ("\t",@codon_data_editing_AtoI_CtoL)."\trhow:$spearman_editing_AtoI_CtoL\ttw:$t_value_editing_AtoI_CtoL\tpvw:$pvalue_editing_AtoI_CtoL\n";
		print OUT "\tCodon data superwobbling\t".join ("\t",@codon_data_superwobbling)."\trhow:$spearman_superwobbling\ttw:$t_value_superwobbling\tpvw:$pvalue_superwobbling\n";
		print OUT2 "$name\t$codon_gene{$name}\t$division\t$classification{$species}\t$spearman\t$t_value\t$pvalue\t$spearman_wobbling\t$t_value_wobbling\t$pvalue_wobbling\t$spearman_editing_AtoI_CtoL\t$t_value_editing_AtoI_CtoL\t$pvalue_editing_AtoI_CtoL\t$spearman_superwobbling\t$t_value_superwobbling\t$pvalue_superwobbling\n";
		print OUT3 "NAME	GENE	DIVISION	CLASSIFICATION	@AAlist2\n" unless $proc_codon > 1 ;
		print OUT3 "$name\t$codon_gene{$name}\t$division\t$classification{$species}\t";
		foreach (0..$#trna_data) {
			print OUT3 "$trna_data[$_]\t$trna_data_wob[$_]\t$trna_data_superwob[$_]\t$codon_data_scaled[$_]\t";
		}
		print OUT3 "\n";
	} elsif ($answer == 2 or $answer == 3) {
		print OUT "$name\t$codon_gene{$name}\t$division\t$classification{$species}\n";
		print OUT "\tAA codon\t".join ("\t",@AAlist)."\n";
		print OUT "\ttRNA count\t".join ("\t",@trna_count)."\n\tCodon count\t".join ("\t",@codon_count)."\n\ttRNA data\t".join ("\t",@trna_data)."\n\ttRNA data_wob\t".join ("\t",@trna_data_wob).."\n\ttRNA data_superwob\t".join ("\t",@trna_data_superwob)."\n";
		print OUT "\tCodon data SCALED\t".join ("\t",@codon_data_scaled)."\n";	
		print OUT "\tCodon data\t".join ("\t",@codon_data)."\trho:$spearman\tpv:$pvalue\n";
		print OUT "\tCodon data Wobbling\t".join ("\t",@codon_data_wobbling)."\trhow:$spearman_wobbling\tpvw:$pvalue_wobbling\n";
		print OUT "\tCodon data editing_AtoI_CtoL\t".join ("\t",@codon_data_editing_AtoI_CtoL)."\trhow:$spearman_editing_AtoI_CtoL\tpvw:$pvalue_editing_AtoI_CtoL\n";		
		print OUT "\tCodon data superwobbling\t".join ("\t",@codon_data_superwobbling)."\trhow:$spearman_superwobbling\tpvw:$pvalue_superwobbling\n";
		print OUT2 "$name\t$codon_gene{$name}\t$division\t$classification{$species}\t$spearman\t$pvalue\t$spearman_wobbling\t$pvalue_wobbling\t$spearman_editing_AtoI_CtoL\t$pvalue_editing_AtoI_CtoL\t$spearman_superwobbling\t$pvalue_superwobbling\n";
		print OUT3 "NAME	GENE	CLASSIFICATION	@AAlist2\n" unless $proc_codon > 1 ;
		print OUT3 "$name\t$codon_gene{$name}\t$division\t$classification{$species}\t";		
		foreach (0..$#trna_data) {
			print OUT3 "$trna_data[$_]\t$trna_data_wob[$_]\t$trna_data_superwob[$_]\t$codon_data_scaled[$_]\t";
		}	
		print OUT3 "\n";		
	} 

}
	
$R->stop() if ($answer == 2);

close OUT;
close OUT2;
close OUT3;

my $end_time   = new Benchmark;
my $difference = timediff($end_time, $start_time);
print "It took ", timestr($difference), "\n";

print "\n\n****DONE****\n\n";

sub spearman_test_R {
	

	my ($trna_ref, $cod_ref, $R) = @_;
	
	# my $start_time = new Benchmark;	
	my $x = join(',',@$trna_ref);
	my $y = join(',',@$cod_ref);
	
	################################## IMPORTANT ####################################
	# The R instance is initialized in the main program to speed up the procedure ###
	
	# my $R = Statistics::R->new();
	# $R->startR;
	# $R->run(q`library(coin)`);
	#################################################################################	
	
	# $R->set('x', $x);
	# $R->set('y', $y);
	
	$R->run(
		q`rm(list=ls())`,
		qq`x<-c($x)`,
		qq`y<-c($y)`,
		q`coeff<-cor.test(x,y, method="spearman", alternative="greater")`, #rho value
		q`stat<-spearman_test(x~y,alternative="greater", distribution=approximate(B=9999))` #9999 Random Rho + 1 observed rho = distribution of 10000 rho values
		);

	my $rho = $R->get('coeff$estimate[[1]]');		
	$rho = @{$rho} if (ref($rho) eq "ARRAY");
	my $pvalue = $R->get('pvalue(stat)[1]');
	$pvalue = @{$pvalue} if (ref($pvalue) eq "ARRAY");
	################################## IMPORTANT ####################################
	# R is stopped in the main program to speed up the procedure ###
	# $R->stopR() if ($answer == 2);
	#################################################################################
	
	# my $end_time   = new Benchmark;
	# my $difference = timediff($end_time, $start_time);
	# print "Benchmark: it took ", timestr($difference), "\n";
	# print "rho $rho pvalue $pvalue\n";
	return ($rho, $pvalue);
}

sub spearman_test {

	my ($trna_ref, $cod_ref) = @_;

	my $c = Statistics::RankCorrelation->new( $trna_ref, $cod_ref, sorted => 1 );
	my $observed_r = $c->spearman;

	my $maxnumber=999; #999 Random Rho + the true rho = 1000

	my $i=0;
	my $sum = 0;

	while ($i<$maxnumber) {
		my @perm = shuffle @$trna_ref;
		# my @perm = shuffle_array @$trna_ref;
		# print "$i\.	|@perm|\n";
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
	return ($observed_r, $pvalue);
}
  
sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTUacgtu/TGCAAtgcaa/;
        return $revcomp;
}

<>;