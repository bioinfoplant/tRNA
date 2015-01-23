#!/usr/local/bin/perl

#########
# Developed by Mattia Belli 2013
#########
use warnings;
use strict;
#use diagnostics;


#3.1 Added RSCU calculation

print "File GenBank FLAT da analizzare: ";
my $file = <STDIN>;
chomp $file;

print "Specify the Gene to analyze, if any: ";
my $gene_to_analyze = <STDIN>;
chomp $gene_to_analyze;
$gene_to_analyze = undef unless $gene_to_analyze =~m/\w/;

my $RSCU_switch = 0;
# my $ALL_switch = 1;

open (DATA, $file) or die "Impossibile aprire $file .";
open (OUT, ">", "$file - CODONW TOT.txt") or die "Impossibile scrivere il file con le sequenze.";
open (OUT1, ">", "$file - CODONW TOT R Ready.txt") or die "Impossibile scrivere il file con le sequenze.";
open (OUT2, ">", "$file - CODONW TOT RSCU.txt") or die "Impossibile scrivere il file con le sequenze." if $RSCU_switch;
open (OUT3, ">", "$file - CODONW $gene_to_analyze.txt") or die "Impossibile scrivere il file con le sequenze." if $gene_to_analyze;
open (OUT4, ">", "$file - CODONW $gene_to_analyze R Ready.txt") or die "Impossibile scrivere il file con le sequenze." if $gene_to_analyze;
print OUT  "NAME	";
print OUT1  "NAME	";
print OUT2  "NAME	" if $RSCU_switch;
print OUT3  "NAME	" if $gene_to_analyze;
print OUT4  "NAME	" if $gene_to_analyze;
my $n;
my $accession_data;

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
   'STOP' => [qw/TAA TAG TGA/],
   
   'Ile' => [qw/ATA ATC ATT/],
   'Met' => [qw/ATG/],
   
   'Cys' => [qw/TGC TGT/],
   'Trp' => [qw/TGG/],
   # 'SelCys' => [qw/TCA/]
);
# my %AA_isoacceptor_types = (	
   # 'Ala' ,4,
   # 'Gly' ,4,
   # 'Pro' ,4,
   # 'Thr' ,4,
   # 'Val' ,4,
   
   # 'Ser' ,6,
   # 'Arg' ,6,
   # 'Leu' ,6,
   
   # 'Phe' ,2,
   
   # 'Asn' ,2,
   # 'Lys' ,2,
   
   # 'Asp' ,2,
   # 'Glu' ,2,
   
   # 'His' ,2,
   # 'Gln' ,2,
   
   # 'Tyr' ,2,
   # 'STOP' ,3,
   
   # 'Ile' ,3,
   # 'Met' ,1,
   
   # 'Cys' ,2,
   # 'Trp' ,1,
   # 'SelCys' ,1,
# );

foreach (sort keys %AA2codons){
	my $AA = $_;
	my $ref = $AA2codons{$_};
	foreach (sort @$ref){	
		print OUT "$_ $AA	";
		print OUT1 "$_ $AA	";
		print OUT2 "$_ $AA	" if $RSCU_switch;
		print OUT3 "$_ $AA	" if $gene_to_analyze;
		print OUT4 "$_ $AA	" if $gene_to_analyze;
	}
}

print OUT "\n";
print OUT1 "\n";
print OUT2 "\n" if $RSCU_switch;
print OUT3 "\n" if $gene_to_analyze;
print OUT4 "\n" if $gene_to_analyze;


chdir "CodonW" or die "Impossibile aprire cartella CodonW";

while (<DATA>) {	
	$accession_data .= $_; 
	next unless ($_ =~ m|^\/\/|);	#Carica i dati di input e li elabora 1 alla volta
	++$n;
	

	
	my ($name, $organism, $date, $id, $sequence) = '';
	if ($accession_data =~ m|VERSION\s+(.+?)\s|g){
		$id = $1;
	} 	
	
	# if ($accession_data =~ m|DEFINITION\s+(.+),|){
		# $name = $1;
	# }
	if ($accession_data =~ m|ORGANISM\s+(.+?)\n(.+?)\.\n|sg){
		$name = $1;
		$organism = $2;
	}
	
	$name .= " $id"; 
	
	#Skips if there are not CDS annotations
	unless ($accession_data =~ m!\n\s{5,}(CDS\s{5,}.+?)\/translation!sg) {
		print "$name -> no CDSs found!\n";
		undef $accession_data;
		next;
	}
	
	$organism =~ s/^[\s]+//g;
	$organism =~ s/\s{2,}/ /g;
	
	$sequence = $1 if ($accession_data =~ m/ORIGIN([\W\w]+)\n\/\//);
	$sequence =~ s/[\W\d]+//g;
	
	
	print "Processing ($n)..$name ";
	open (all_CDS_seq, ">", "all_CDS_SEQ.txt") or die "Impossibile scrivere il file con le sequenze.";
	if ($gene_to_analyze) {
		open (gene_CDS_seq, ">", "gene_CDS_SEQ.txt") or die "Impossibile scrivere il file con le sequenze.";	
	}
	
	my @name_list;
	my $n_cds;
	my $found_gene;
	my $gene_found;
	while ($accession_data =~ m!\n\s{5,}(CDS\s{5,}.+?)\/translation!sg){
		++$n_cds;

		my $cds = $1;
		my $seq_position;
		my $annotation;
		my $CDS_seq;
		my $gene_name = 'Unknown';
		$gene_name = $1 if ($cds =~ m|\/db_xref="(.+)"|);
		$gene_name = $1 if ($cds =~ m|\/locus_tag="(.+)"|);
		$gene_name = $1 if ($cds =~ m|\/protein_id="(.+)"|);
		$gene_name = $1 if ($cds =~ m|\/gene="(.+)"|);
	
		if ($cds =~ m|(CDS\s+(?:complement\()?join\((.+?)\))\n|s){
			$annotation = $1;
			$seq_position = $2;
			$annotation =~ s/\n\s+//g;
			my @parts = split (",", $seq_position);
			foreach (@parts){
				my ($start, $end) = '';
				if ($_ =~ m|complement\((\d+)\.\.(\d+)\)|){
					$start = $1;
					$end = $2;
					my $part_seq = substr ($sequence, $start-1, $end-$start+1);
					$CDS_seq .= &reverse_complement($part_seq);
				} elsif ($_ =~ m|^(\d+)\$|) {
					$CDS_seq .= substr ($sequence, $1, 1);
				} elsif ($_ =~ m|(\d+)\.\.(\d+)|){
					$start = $1;
					$end = $2;
					my $part_seq = substr ($sequence, $start-1, $end-$start+1);
					$CDS_seq .= $part_seq;
				}		
			}	
		} elsif ($cds =~ m|(CDS.+?>?(\d+)\.\.>?(\d+))|){
			$annotation = $1;
			my $start = $2;
			my $end = $3;
			$seq_position = "$start..$end";
			$CDS_seq = substr ($sequence, $start-1, $end-$start+1);
		}
		if ($cds =~ m|CDS\s+complement|) {
			$CDS_seq = &reverse_complement($CDS_seq);		
		}
		$annotation =~ s/\s{2,}/ /g;
		
		print all_CDS_seq ">$name $gene_name $annotation\n$CDS_seq\n\n";

		my $regex = $gene_to_analyze if $gene_to_analyze;
		$regex =~ s/\s+?/\|/g;
		if ($gene_to_analyze and $gene_name =~ m|$regex|i){
			++$found_gene;
			push (@name_list, "$name [$gene_name]");
			print gene_CDS_seq ">$name $gene_name $annotation\n$CDS_seq\n\n";
		}
		# print ">$gene_name $annotation\n$CDS_seq\n\n";
		# <>;
	}
	print "$n_cds CDSs FOUND\n";
	close all_CDS_seq;
	close gene_CDS_seq if ($gene_to_analyze);	

	
	system ('CodonW.exe all_CDS_SEQ.txt -nomenu -silent -total >NUL 2>&1');
	open (CODONW, 'all_CDS_SEQ.blk') or die "Impossibile aprire il file output di CODONW.";
	my $codonW = join ("",<CODONW>);
	close CODONW;
	
	# my %tot_codon_per_AA;
	my %codon_table;
	my %codon_table_rscu;
	# print $codonW;
	while ($codonW =~ m!([AUCG]+)(?:\s+)?(\d+)\s*([\d\.]+)!g){
		my $triplet = $1;
		my $codon_number = $2;
		my $rscu = $3;
		# print "\n$triplet";
		$triplet =~ s/U/T/g;
		# print "\n$triplet";
		# print "\n$codon2AA{$triplet}";
		$codon_table{$triplet} = $codon_number;
		$codon_table_rscu{$triplet} = $rscu;
		# print "\n$codon_table{$triplet}";
		# $tot_codon_per_AA{$codon2AA{$triplet}} += $codon_number;
		# print "\n$tot_codon_per_AA{$codon2AA{$triplet}}";
	}
	print OUT "$name	";
	print OUT1 "$name	";
	print OUT2 "$name	" if $RSCU_switch;
	
	foreach my $AA (sort keys %AA2codons){
		my $ref = $AA2codons{$AA};
		foreach (sort @$ref){		#Da aggiustare perchè così mi fa il sort degli anticodoni e poi li trasforma in anticodoni
			my $triplet = $_;
			# my $rscu_score;
			# if ($tot_codon_per_AA{$AA} == 0) {
				# $rscu_score = 0;
			# } else {
				# $rscu_score = ($codon_table{$triplet}/$tot_codon_per_AA{$AA})*(scalar @$ref);
			# }
			print OUT "$triplet $codon_table{$triplet}	";
			print OUT1 "$codon_table{$triplet}	";
			print OUT2 "$triplet $codon_table_rscu{$triplet}	" if $RSCU_switch;
			# print "$triplet $codon_table_rscu{$triplet}	";
			# print "\n\n$AA";
			# print "\n$triplet";
			# print "\n$codon_table{$triplet}";
			# print "\n$tot_codon_per_AA{$AA}";		
			# print "\n ", scalar @$ref;
			# print "\n$rscu_score";	
			# <>;
		}
	}	

	print OUT "\n";
	print OUT1 "\n";
	
	print OUT2 "\n" if $RSCU_switch;
	
	
	if ($found_gene) {
		system ('CodonW.exe gene_CDS_SEQ.txt -nomenu -silent >NUL 2>&1');
		open (CODONW2, 'gene_CDS_SEQ.blk') or die "Impossibile aprire il file output di CODONW.";
		
		my $CUD_number;
		my $codonW_data;
		while (<CODONW2>) {	
			$codonW_data .= $_; 
			next unless ($_ =~ m|Genetic code|);	#Carica i dati di input e li elabora 1 alla volta
			++$CUD_number;
			my %codon_table_all;
			while ($codonW_data =~ m!([AUCG]+)(?:\s+)?(\d+)\s*([\d\.]+)!g){
				my $triplet = $1;
				my $codon_number = $2;
				$triplet =~ s/U/T/g;
				$codon_table_all{$triplet} = $codon_number;
			}
			print OUT3 $name_list[$CUD_number-1]."	";
			print OUT4 $name_list[$CUD_number-1]."	";
			
			foreach my $AA (sort keys %AA2codons){
				my $ref = $AA2codons{$AA};
				foreach (sort @$ref){		#Da aggiustare perchè così mi fa il sort degli anticodoni e poi li trasforma in anticodoni
				my $triplet = $_;
				print OUT3 "$triplet $codon_table_all{$triplet}	";
				print OUT4 "$codon_table_all{$triplet}	";
				}
			}
			# foreach (sort keys %codon_table_all){
				# print OUT3 "$_ $codon_table_all{$_}	";
			# }

			print OUT3 "\n";
			print OUT4 "\n";
			undef $codonW_data;
		}
		close CODONW2;
	}
	undef $accession_data;
}


chdir;
	# while ($codon_usage =~ m!([AUCG]+)(?:\s+)?([\d\.]+)\(\s*(\d+)\)!g){
# CodonW.exe CDS_SEQ.txt -nomenu -silent -total

close DATA;
close OUT;
close OUT1;
close OUT2 if $RSCU_switch;
close OUT3 if $gene_to_analyze;
close OUT4 if $gene_to_analyze;

sub reverse_complement {
        my $dna = shift;
	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

print "\n\n****FATTO****\n\n";
<>;