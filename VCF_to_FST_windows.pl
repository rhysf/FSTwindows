#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_FASTA;
use read_VCF_lines;
use read_Tab;
use sliding_windows;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -a <Name tab Lineage file> -b <Name tab Type tab Location file> -c <Lineage 1 comparison> -d <Lineage 2 comparison> -r <reference fasta> > out.window.tab\n
Optional: -w\tWindow length [10000]
	  -e\tExclude variants on contig/chromosome []
          -z\tExclude variants not on contig/chromosome []
	  -o\tOutput to stdout (o), or output to file (f) [o]
	  -f\tIf -o f, then filename to output [opt_c-opt_d-FST-wl-opt_w.tab]\n";
our($opt_a, $opt_b, $opt_c, $opt_d, $opt_e, $opt_f, $opt_o, $opt_r, $opt_w, $opt_z);
getopt('abcdeforwz');
die $usage unless ($opt_a && $opt_b && $opt_r && $opt_c && $opt_d);
if(!defined $opt_w) { $opt_w = 10000; }
if(!defined $opt_o) { $opt_o = 'o'; }
if(!defined $opt_f) { $opt_f = ($opt_c . '-' . $opt_d . '-FST-wl-' . $opt_w . '.tab'); }
die "-w needs to be numerical ($opt_w)\n" if (($opt_w) !~ /\d+/);
die "-o needs to be o or f ($opt_o)\n" if (($opt_o ne 'o') && ($opt_o ne 'f'));

# Save Name -> Lineages, and Name->VCF-Location
my $names_to_lineage = tabfile::save_one_to_two_columns($opt_a, 0, 1);
my $names_to_location = tabfile::save_two_to_three_columns($opt_b, 0, 1, 2);
my ($lineage1_bases, $lineage2_bases, $lineage1_and_2_bases) = &save_name_to_location_for_2_lineages($names_to_lineage, $names_to_location, $opt_c, $opt_d);

# Initial calculations for FST
my $r_the_number_of_populations = 2;
my $r_minus_1 = ($r_the_number_of_populations - 1);
my $ni_pop1 = (scalar(keys(%{$lineage1_bases})) * 2);
my $ni_pop2 = (scalar(keys(%{$lineage2_bases})) * 2);
my $sum_of_ni = ($ni_pop1 + $ni_pop2);
my $ni_pop1_squared = ($ni_pop1 ** 2);
my $ni_pop2_squared = ($ni_pop2 ** 2);
my $sum_of_ni_squared_ind = ($ni_pop1_squared + $ni_pop2_squared);
my $n_bar = ($sum_of_ni / $r_the_number_of_populations);
my $start_of_n_sub_c = (1 / $r_minus_1);
my $n_sub_c = ($start_of_n_sub_c * ($sum_of_ni - ($sum_of_ni_squared_ind / $sum_of_ni)));

# Results of initial calculations
warn "r = $r_the_number_of_populations\n";
warn "ni pop1 = $ni_pop1, ni pop2 = $ni_pop2, sum = $sum_of_ni (for just the 2 populations of interested specified by $opt_c and $opt_d)\n";
warn "ni pop1^2 = $ni_pop1_squared, ni pop2^2 = $ni_pop2_squared, sum = $sum_of_ni_squared_ind (2 isolates spec by $opt_c and $opt_d)\n";
warn "n bar = $n_bar\n";
warn "n sub c = $n_sub_c\n";

# Save sequences
my ($sequences, $descriptions, $orderofcontigs) = fastafile::fasta_id_to_seq_hash($opt_r);

# Make the windows
my $want = "FST";
my %isolate_names_fake;
$isolate_names_fake{'name'} = 1;
my ($windows, $total_GC)= slidingwindows::make_windows_for_seq_hash($sequences, $opt_w, \%isolate_names_fake, $want);

# Open output
my $ofh;
if($opt_o eq 'f') { open $ofh, '>', $opt_f or die "Cannot open outfile $opt_f : $!\n"; }

# Save all variants in the entire population
my (%Positions_of_variants);
ALL: foreach my $isolate(keys %{$lineage1_and_2_bases}) {
	# VCF
	my $vcf_file = $$lineage1_and_2_bases{$isolate}{'VCF'};
	#warn "gots $isolate and $vcf_file\n";
	open my $fh, '<', $vcf_file or die "Cannot open $vcf_file : $!\n";
	warn "Reading $vcf_file for $isolate (all)...\n";
	VCF: while (my $line = <$fh>) {
		chomp $line;
		my ($VCF_line) = vcflines::read_VCF_lines($line);
		next VCF if($$VCF_line{'next'} eq 1);
		next VCF if(($$VCF_line{'base_type0'} eq 'deletion') || ($$VCF_line{'base_type0'} eq 'insertion'));
		next VCF if(($$VCF_line{'base_type0'} eq 'het_insertion') || ($$VCF_line{'base_type0'} eq 'het_deletion'));
		next VCF if($$VCF_line{'base_type0'} eq 'ambiguous');

		# new for genome/mito for example
		next VCF if(($opt_e) && ($$VCF_line{'supercontig'} ne $opt_e));
		next VCF if(($opt_z) && ($$VCF_line{'supercontig'} eq $opt_z));

		# Save the bases to scalar
		my ($base1, $base2);
		if($$VCF_line{'base_type0'} eq 'reference') {
			$base1 = $$VCF_line{'reference_VCF_format'};
			$base2 = $$VCF_line{'reference_VCF_format'};
		} elsif($$VCF_line{'base_type0'} eq 'snp') {
			$base1 = $$VCF_line{'consensus_VCF_format'};
			$base2 = $$VCF_line{'consensus_VCF_format'};
		} elsif($$VCF_line{'base_type0'} eq 'heterozygous') {
			($base1, $base2) = vcflines::resolve_ambiguity_char($$VCF_line{'amb_char0'});
		}

		# Save the bases to hash
		if((!defined $base1) || (!defined $base2)) { warn "Ignoring this variant: $line\n"; }
		else {
			$Positions_of_variants{$$VCF_line{'supercontig'}}{$$VCF_line{'position'}}{'1'}{$isolate} = $base1;
			$Positions_of_variants{$$VCF_line{'supercontig'}}{$$VCF_line{'position'}}{'2'}{$isolate} = $base2;
		}
	}
	close $fh;

	# Reference bases in mpileup format (if specified seperately from VCF)
	my $ref_file = $$lineage1_and_2_bases{$isolate}{'REF'};
	next ALL if(!defined $ref_file);
	open IN2, "<$ref_file" or die "Cannot open $ref_file\n";
	REF: while (my $line = <IN2>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($contig, $position, $ref_base) = @bits;
		die "$contig $position is found in VCF and REF for $isolate\n" if(defined $Positions_of_variants{$contig}{$position});
		$Positions_of_variants{$contig}{$position}{'1'}{$isolate} = $ref_base;
		$Positions_of_variants{$contig}{$position}{'2'}{$isolate} = $ref_base;
	}
	close IN2;
}

# Remove bases that are not called in every isolate 
warn "Removing any loci that are not verified in all\n";
foreach my $contig(keys %Positions_of_variants) {
	VARIANTPOS: foreach my $position(keys %{$Positions_of_variants{$contig}}) {
		# Found in everything
		next VARIANTPOS if(scalar(keys(%{$Positions_of_variants{$contig}{$position}{'1'}})) eq scalar(keys(%{$lineage1_and_2_bases})));
		next VARIANTPOS if(scalar(keys(%{$Positions_of_variants{$contig}{$position}{'2'}})) eq scalar(keys(%{$lineage1_and_2_bases})));
		# Not suitable
		delete $Positions_of_variants{$contig}{$position};
	}
}

# Go through each lineage finding inter-population variation
warn "Going through getting inter lineage population variation...\n";
foreach my $contig(keys %{$$windows{'name'}}) {
	warn "Calculating FST values for w.l. $opt_w in $contig\n";
	foreach my $start_window(sort keys %{$$windows{'name'}{$contig}}) {
		foreach my $end_window(sort keys %{$$windows{'name'}{$contig}{$start_window}}) {
			#warn "$start_window -> $end_window\n";
			
			my ($sum_Q2, $sum_Q3, $sum_1_minus_Q3);
			POSITIONS: foreach my $position(keys %{$Positions_of_variants{$contig}}) {
				next POSITIONS if (($position < $start_window) || ($position > $end_window));
				#warn "$position inside window $start_window - $end_window\n";
		
				# Save variants at this allele for each population
				my (%pop1_variants_at_allele, %pop2_variants_at_allele, %pop1_and_2_variants_at_allele);
				foreach my $isolate(keys %{$lineage1_bases}) {
					if(!defined $Positions_of_variants{$contig}{$position}{'1'}{$isolate}) {
						die "Cannot find an expected entry in $contig $position 1 $isolate\n";
					}
					my $variant1 = $Positions_of_variants{$contig}{$position}{'1'}{$isolate};
					my $variant2 = $Positions_of_variants{$contig}{$position}{'2'}{$isolate};
					$pop1_variants_at_allele{$variant1}++;
					$pop1_variants_at_allele{$variant2}++;
					$pop1_and_2_variants_at_allele{$variant1}=1;
					$pop1_and_2_variants_at_allele{$variant2}=1;
				}
				foreach my $isolate(keys %{$lineage2_bases}) {
					my $variant1 = $Positions_of_variants{$contig}{$position}{'1'}{$isolate};
					my $variant2 = $Positions_of_variants{$contig}{$position}{'2'}{$isolate};
					$pop2_variants_at_allele{$variant1}++;
					$pop2_variants_at_allele{$variant2}++;
					$pop1_and_2_variants_at_allele{$variant1}=1;
					$pop1_and_2_variants_at_allele{$variant2}=1;
				}

				# Calculate X and Y (for each allele)
				my ($X, $Y);
				foreach my $variant(keys %pop1_and_2_variants_at_allele) {
					#warn "Variant $variant\n";
					my $frequency_times_ni = 0;
					if(defined $pop1_variants_at_allele{$variant}) {
						my $pop1_tally = $pop1_variants_at_allele{$variant};
						my $pop1_frequency = ($pop1_tally / $ni_pop1);
						$frequency_times_ni += ($pop1_frequency * $ni_pop1);
						my $pop1_frequency_squared = ($pop1_frequency ** 2);
						$X += ($pop1_frequency_squared * $ni_pop1);
						#warn "Pop1 has frequency $pop1_frequency, frequency squared = $pop1_frequency_squared, X now = $X\n";
					}
					if(defined $pop2_variants_at_allele{$variant}) {
						my $pop2_tally = $pop2_variants_at_allele{$variant};
						my $pop2_frequency = ($pop2_tally / $ni_pop2);
						$frequency_times_ni += ($pop2_frequency * $ni_pop2);
						my $pop2_frequency_squared = ($pop2_frequency ** 2);
						$X += ($pop2_frequency_squared * $ni_pop2);
						#warn "Pop2 has frequency $pop2_frequency, frequency squared = $pop2_frequency_squared, X now = $X\n";
					}
					$Y += ($frequency_times_ni ** 2);
				}
				#warn "X and Y for this loci = $X and $Y\n";

				# Calculate Q2 and Q3
				my $Q2 = (($X - $r_the_number_of_populations) / ($r_the_number_of_populations * ($n_bar - 1)));
				#warn "Q2 = $Q2\n";

				my $Q3_fir_fir_part = (1 / (($r_the_number_of_populations * $r_minus_1) * $n_bar * $n_sub_c));
				my $Q3_end_fir_part = ($Y - ((($n_bar * ($n_sub_c - 1)) / ($n_bar - 1)) * $X));
				my $Q3_fir_end_part = (($n_bar - $n_sub_c) / ($n_sub_c * ($n_bar - 1)));
				my $Q3_end_end_part = (1 - ((1 / $r_minus_1) * $X));
				my $Q3 = (($Q3_fir_fir_part * $Q3_end_fir_part) + ($Q3_fir_end_part * $Q3_end_end_part));
				#warn "Q3 parts = $Q3_fir_fir_part $Q3_end_fir_part $Q3_fir_end_part $Q3_end_end_part\nQ3 = $Q3\n";
		
				$sum_Q2 += $Q2;
				$sum_Q3 += $Q3;
				$sum_1_minus_Q3 += (1-$Q3);
				#my $FST_of_each_loci = (($Q2 - $Q3) / (1-$Q3));
				#warn "FST of each loci = $FST_of_each_loci\n";
			}
			#warn "\nsum Q2: $sum_Q2\n";
			#warn "sum Q3: $sum_Q3\n";
			#warn "sum 1 - Q3: $sum_1_minus_Q3\n";
			
			# For multiple loci Q2 and Q3 are summed across loci before calculation
			my $FST = 0;
			if(defined $sum_1_minus_Q3) {
				if($sum_1_minus_Q3 ne 0) { $FST = (($sum_Q2 - $sum_Q3) / $sum_1_minus_Q3); }
			}
			my $rounded_FST = sprintf "%.2f", $FST;
			#warn "FST = $FST -> $rounded_FST\n";
			#$windows{$contig}{$start_window}{$end_window}{'FST'} = $rounded_FST;
			$$windows{'name'}{$contig}{$start_window}{$end_window}{'FST'} = $FST;
		}
	}
}

# Print windows
warn "Printing windows...\n";
if($opt_o eq 'o') { print "Contig\tWindow_Start\tWindow_Stop\tRunning_Position\tFST\n"; }
else { print $ofh "Contig\tWindow_Start\tWindow_Stop\tRunning_Position\tFST\n"; }
my $running_window_length;
foreach my $contigs(@{$orderofcontigs}) {
	STARTS: foreach my $starts(sort {$a <=> $b} keys %{$$windows{'name'}{$contigs}}) {
		STOPS: foreach my $stops(sort {$a <=> $b} keys %{$$windows{'name'}{$contigs}{$starts}}) {
			my $FST = $$windows{'name'}{$contigs}{$starts}{$stops}{'FST'};
			$running_window_length += ($stops - $starts);
			if($opt_o eq 'o') { print "$contigs\t$starts\t$stops\t$running_window_length\t$FST\n"; }
			else { print $ofh "$contigs\t$starts\t$stops\t$running_window_length\t$FST\n"; }
		}
	}
}

sub save_name_to_location_for_2_lineages {
	my ($names2lineage, $names2location, $lin1, $lin2) = @_;

	my (%lineage1_bases, %lineage2_bases);
	NAME: foreach my $name(keys %{$names2location}) {
		my $lineage = $$names2lineage{$name};
		die "lineage not defined for $name : $!\n" if(!defined $lineage);
		next NAME if(($lineage ne $lin1) && ($lineage ne $lin2));
		TYPE: foreach my $type(keys %{$$names2location{$name}}) {
			next TYPE if($type ne 'VCF');
			my $location = $$names2location{$name}{$type};
			if($lineage eq $opt_c) { $lineage1_bases{$name}{'VCF'} = $location; }
			else { $lineage2_bases{$name}{'VCF'} = $location; }
		}
	}
	my %lineage1_and_2_bases = (%lineage1_bases, %lineage2_bases);
	return (\%lineage1_bases, \%lineage2_bases, \%lineage1_and_2_bases);
}
