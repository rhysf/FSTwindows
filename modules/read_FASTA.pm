package fastafile;
use strict;
use Bio::SeqIO;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use FindBin qw($Bin);
use lib "$Bin";
use read_Tab;
use Data::Dumper;

### rfarrer@broadinstitute.org

sub fasta_to_struct {
	my $input = $_[0];
	my %struct;
	$struct{'filename'} = $input;
	warn "fasta_to_struct: saving from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		my $seq = $seq_obj->seq;
		my $desc = $seq_obj->description;

		# Save
		$struct{'seq'}{$id} = $seq;
		$struct{'desc'}{$id} = $desc;
		push @{$struct{'order'}}, $id;
	}
	return \%struct;
}

sub fasta_id_to_seq_hash {
	my $input = $_[0];
	my (%sequences, %descriptions);
	my @order;
	warn "fasta_id_to_seq_hash: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		my $seq = $seq_obj->seq;
		my $desc = $seq_obj->description;
		$sequences{$id}=$seq;
		$descriptions{$id}=$desc;
		push @order, $id;
	}
	return (\%sequences, \%descriptions, \@order);
}

sub fasta_id_to_seq_hash_1_entry {
	my ($input, $want) = @_;
	my (%sequences);
	warn "fasta_id_to_seq_hash_1_entry: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	FASTA: while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
    		my $seq = $seq_obj->seq;
		next FASTA if ($want ne $id);
		$sequences{$id}=$seq;
	}
	return (\%sequences);
}

sub fasta_id_to_seq_length_hash {
	my $input = $_[0];
	my (%lengths);
	warn "fasta_id_to_seq_length_hash: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
    		my $seq = $seq_obj->seq;
		my $length = length($seq);
		$lengths{$id} = $length;
	}
	return (\%lengths);
}

sub fasta_seq_length_to_id_hash {
	my $input = $_[0];
	my (%lengths);
	warn "fasta_seq_length_to_id_hash: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
    		my $seq = $seq_obj->seq;
		my $length = length($seq);
		$lengths{$length} = $id;
	}
	return (\%lengths);
}

sub fasta_struct_seq_length_to_id_hash {
	my $fasta_struct = $_[0];
	my %lengths;
	warn "fasta_struct_seq_length_to_id_hash...\n";
    	foreach my $id(keys %{$$fasta_struct{'seq'}}) {
		my $seq = $$fasta_struct{'seq'}{$id};
		my $length = length($seq);
		$lengths{$length} = $id;
	}
	return (\%lengths);
}

sub fasta_id_to_seq_length_array {
	my $input = $_[0];
	my (@lengths);
	warn "fasta_id_to_seq_length_array: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
    		my $seq = $seq_obj->seq;
		my $length = length($seq);
		push @lengths, length($seq);	
	}
	return (\@lengths);
}

sub fasta_to_total_seq_length {
	my $input = $_[0];
	my ($lengths) = 0;
	warn "fasta_to_total_seq_length: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
    		my $seq = $seq_obj->seq;
		$lengths += length($seq);
	}
	return ($lengths);
}

sub fasta_struct_to_total_length {
	my ($fasta_struct) = $_[0];
	warn "fasta_struct_to_total_length...\n";
	my $total_length = 0;
	FASTA: foreach my $id(keys %{$$fasta_struct{'seq'}}) {
		my $seq = $$fasta_struct{'seq'}{$id};
		my $length = length($seq);
		$total_length += $length;
	}
	return ($total_length);
}

sub fasta_struct_truncate {
	my ($fasta_struct, $truncate_from_5prime, $truncate_from_3prime, $minimum_length) = @_;

	my ($trimmed_count, $skipped_count) = (0, 0);
	my @new_order;
	FASTA: foreach my $id(keys %{$$fasta_struct{'seq'}}) {
		my $seq = $$fasta_struct{'seq'}{$id};
		my $length = length($seq);

		# truncate 5' and 3'
		my $truncated_length = ($length - $truncate_from_3prime);
		my $truncated_seq = substr $seq, $truncate_from_5prime, $truncated_length;
		if($truncated_length < $minimum_length) {
			$skipped_count++;
			delete $$fasta_struct{'seq'}{$id};
			delete $$fasta_struct{'desc'}{$id};
			next FASTA;
		}
		$$fasta_struct{'seq'}{$id} = $truncated_seq;
		$trimmed_count += (length($seq) - length($truncated_seq));
		push @new_order, $id;
	}
	@{$$fasta_struct{'order'}} = @new_order;
	warn "Trimmed $trimmed_count bases.\n";
	warn "Skipped $skipped_count entries for < $minimum_length\n";
	return ($fasta_struct);
}

sub init_fasta_struct {
	my $fasta_struct = $_[0];
	my %new_struct;
	$new_struct{'filename'} = $$fasta_struct{'filename'};
	return \%new_struct
}

sub fasta_struct_append_id {
	my ($fasta_struct, $append_id) = @_;
	my $new_struct = &init_fasta_struct($fasta_struct);

	warn "fasta_struct_append_id: $append_id\n";
	foreach my $id(@{$$fasta_struct{'order'}}) {
		my $new_id = ($id . $append_id);

		# Save
		$$new_struct{'seq'}{$new_id} = $$fasta_struct{'seq'}{$id};
		$$new_struct{'desc'}{$new_id} = $$fasta_struct{'desc'}{$id};
		push @{$$new_struct{'order'}}, $new_id;
	}
	return ($new_struct);
}

sub fasta_struct_replace_id_with_desc {
	my ($fasta_struct, $replace_with_desc_part_number) = @_;
	my $new_struct = &init_fasta_struct($fasta_struct);

	warn "fasta_struct_replace_id_with_desc: $replace_with_desc_part_number\n";
	foreach my $id(@{$$fasta_struct{'order'}}) {
		my @desc_parts = split /\s/, $$fasta_struct{'desc'}{$id};
		my $new_id = $desc_parts[$replace_with_desc_part_number];

		# Save
		$$new_struct{'seq'}{$new_id} = $$fasta_struct{'seq'}{$id};
		$$new_struct{'desc'}{$new_id} = $$fasta_struct{'desc'}{$id};
		push @{$$new_struct{'order'}}, $new_id;
	}
	return ($new_struct);
}

sub fasta_struct_remove_words_from_id {
	my ($fasta_struct, $remove_words) = @_;
	my $new_struct = &init_fasta_struct($fasta_struct);
	my @words = split /,/, $remove_words;

	warn "fasta_struct_remove_words_from_id: $remove_words\n";
	foreach my $id(@{$$fasta_struct{'order'}}) {
		my $new_id = $id;
		foreach(@words) { $new_id =~ s/$_//g; }

		# Save
		$$new_struct{'seq'}{$new_id} = $$fasta_struct{'seq'}{$id};
		$$new_struct{'desc'}{$new_id} = $$fasta_struct{'desc'}{$id};
		push @{$$new_struct{'order'}}, $new_id;
	}
	return ($new_struct);
}

sub fasta_struct_remove_id {
	my ($fasta_struct, $remove_id) = @_;
	my $new_struct = &init_fasta_struct($fasta_struct);

	my $found = 0;
	#warn "fasta_struct_remove_id: $remove_id\n";
	foreach my $id(@{$$fasta_struct{'order'}}) {
		if($id eq $remove_id) {
			$found++;
			next;
		}

		# Save
		$$new_struct{'seq'}{$id} = $$fasta_struct{'seq'}{$id};
		$$new_struct{'desc'}{$id} = $$fasta_struct{'desc'}{$id};
		push @{$$new_struct{'order'}}, $id;
	}
	warn "fasta_struct_remove_id: removed $found ids.\n";
	return ($new_struct);
}

sub fasta_struct_remove_ids_from_text_file {
	my ($fasta_struct, $remove_file) = @_;
	#my $new_struct = &init_fasta_struct($fasta_struct);

	warn "fasta_struct_remove_ids_from_text_file: $remove_file\n";
	#my $remove_ids = tabfile::save_one_column($remove_file, 0);
	my $remove_ids = tabfile::save_columns_to_one_hash($remove_file, 0);
	foreach my $remove_id(keys %{$remove_ids}) { $fasta_struct = &fasta_struct_remove_id($fasta_struct, $remove_id); }
	return ($fasta_struct);
}

sub fasta_struct_remove_description {
	my $fasta_struct = $_[0];
	warn "fasta_struct_remove_description...\n";
	delete $$fasta_struct{'desc'};
	return $fasta_struct;
}

sub fasta_struct_reverse_compliment {
	my ($fasta_struct, $ids) = @_;
	my @reverse_ids = split /,/, $ids;

	foreach my $id(keys %{$$fasta_struct{'seq'}}) {
		foreach(@reverse_ids) {
			if($_ eq $id) {
				 $$fasta_struct{'seq'}{$id} = &reverse_compliment($$fasta_struct{'seq'}{$id});
			}
		}
	}	
	return $fasta_struct;
}

sub fasta_struct_remove_all_but_one_sequence {
	my ($fasta_struct, $keep_sequence) = @_;
	my $new_struct = &init_fasta_struct($fasta_struct);

	warn "fasta_struct_remove_all_but_one_sequence: $keep_sequence\n";
	die "Can't find $keep_sequence\n" unless(defined $$fasta_struct{'seq'}{$keep_sequence});

	# Save
	$$new_struct{'seq'}{$keep_sequence} = $$fasta_struct{'seq'}{$keep_sequence};
	$$new_struct{'desc'}{$keep_sequence} = $$fasta_struct{'desc'}{$keep_sequence};
	push @{$$new_struct{'order'}}, $keep_sequence;
	
	return ($new_struct);
}


sub fasta_hash_reverse_compliment {
	my ($fasta_hash, $ids, $reverse) = @_;
	my %sequences = %{$$fasta_hash};
	if($reverse eq 'y') {
		foreach my $id(keys %{$fasta_hash}) {
			my @reverse_ids = split /\,/, $ids;
			foreach(@reverse_ids) {
				if($_ eq $id) {
					 $sequences{$id} = &reverse_compliment($$fasta_hash{$id});
				}
			}
		}	
	}
	return (\%sequences);
}

sub fasta_struct_reorder {
	my ($fasta_struct, $order_by, $ids, $additional_file) = @_;
	warn "fasta_struct_reorder: $order_by, $ids\n";

	# Reorder n = None/ignore
	if($order_by eq 'n') { return $fasta_struct; }
	delete $$fasta_struct{'order'};

	# Reorder b = from biggest to smallest
	if($order_by eq 'b') { 
		my $lengths = &fasta_struct_seq_length_to_id_hash($fasta_struct);
		foreach my $length(sort {$b<=>$a} keys %{$lengths}) {
			my $id = $$lengths{$length};
			push @{$$fasta_struct{'order'}}, $id;
		}
	}

	# Reorder i = IDs separated by comma
	elsif($order_by eq 'i') { @{$$fasta_struct{'order'}} = split /,/, $ids; }

	# Reorder k = IDS from text file
	elsif($order_by eq 'k') {
		open my $fh, '<', $additional_file or die "Cannot open $additional_file: $!\n";
		while(my $line=<$fh>) {
			chomp $line;
			push @{$$fasta_struct{'order'}}, $line;
		}
		#@{$$fasta_struct{'order'}} = <$fh>; # doesn't chomp
	}

	# Reorder f = IDs from FASTA
	elsif($order_by eq 'f') { 
		my $order_from_fasta = &fasta_id_to_order_array($additional_file); 
		@{$$fasta_struct{'order'}} = @{$order_from_fasta};
	}
	
	else { die "Reorder according to what in $0?: $order_by\n"; }
	return $fasta_struct;
}

sub fasta_hash_print_simple {
	my $fasta_hash = $_[0];
	FASTA: foreach my $id(keys %{$fasta_hash}) {
		my $seq = $$fasta_hash{$id};
		$seq =~ s/(\S{60})/$1\n/g;
		print ">$id\n$seq\n";
	}
	return;
}

sub fasta_hash_print_simple_outfile {
	my ($fasta_hash, $outfile) = @_;
	open my $ofh, '>', $outfile or die "Cannot open $outfile: $!\n";
	FASTA: foreach my $id(keys %{$fasta_hash}) {
		my $seq = $$fasta_hash{$id};
		$seq =~ s/(\S{60})/$1\n/g;
		print $ofh ">$id\n$seq\n";
	}
	close $ofh;
	return;
}

sub fasta_struct_print {
	my ($fasta_struct, $output_format, $data_type, $split_into_files, $outfile_optional_parameter) = @_;
	die "fasta_struct_print: $output_format not recognised\n" unless($output_format =~ m/fasta|tba|fastq|nexus/);
	die "fasta_struct_print: $data_type not recognised for nexus output\n" if(($output_format =~ m/nexus/) && ($data_type !~ m/dna|prot/));

	# Outfile?
	my $ofh;
	if(defined $outfile_optional_parameter) {
		open $ofh, '>', $outfile_optional_parameter or die "Cannot open $outfile_optional_parameter : $!\n";
	}

	# Go through the FASTA struct accordind to the order array
	warn "fasta_struct_print output=$output_format\n";
	my ($sequence_length, $sequence_count, $split_count, $split_name) = (0, 0, 0, 0);
	FASTA: foreach my $id(@{$$fasta_struct{'order'}}) {
		die "no sequence found in FASTA struct for $id\n" if(!defined $$fasta_struct{'seq'}{$id});
		my $seq  = $$fasta_struct{'seq'}{$id};
		my $desc = $$fasta_struct{'desc'}{$id};

		# Nexus format (all are expected to be the same length)
		$sequence_length = length($seq); 
		$sequence_count++;

		# FASTA (base) output
		if($output_format eq 'fasta') { ($split_count, $split_name) = &fasta_struct_print_to_fasta($fasta_struct, $id, $split_into_files, $split_count, $split_name, $ofh); }

		# Threaded-block aligner
		if($output_format eq 'tba') {
			my $new_id= ">$data_type:$id:1:+:$sequence_length";
			if(defined $outfile_optional_parameter) { print $ofh "$new_id\n$seq\n"; }
			else { print "$new_id\n$seq\n"; }
		}

		# FASTQ
		if($output_format eq 'fastq') {
			my $quality;
			for(my $i=0; $i<$sequence_length; $i++) { $quality .= "I"; }
			if(defined $outfile_optional_parameter) { print $ofh '@' . "$id\n$seq\n+\n$quality\n"; }
			else { print '@' . "$id\n$seq\n+\n$quality\n"; }
		}
	}

	# NEXUS output
	if($output_format eq 'nexus') { &print_nexus($$fasta_struct{'seq'}, $sequence_count, $sequence_length, $data_type, $ofh); }
	return 1;
}

sub fasta_summary {
	my ($input, $print_summary_of_full_file, $print_summary_for_each_entry) = @_;
	my %information;
	my @fasta_lengths;
	
	# Print headers for summary per entry
	if($print_summary_for_each_entry ne 'n') {
		# Long
		if($print_summary_for_each_entry eq 'l') { warn "Identity\tDescription\tSequence_length\tA\tC\tT\tG\tN\tU\tGCpercent\n"; }
		# Short
		elsif($print_summary_for_each_entry eq 's') { warn "Identity\tSequence_length\n"; }
		# Unknown
		else { die "fasta_summary in read_FASTA.pm: $print_summary_for_each_entry is not l (long) or s (short)\n"; }
	}

	# Init
	($information{'full_length'}, $information{'n_of_seq'}, $information{'full_length_minus_N'}) = (0, 0, 0);

	#warn "fasta_summary: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
    		my $seq = $seq_obj->seq;
		my $desc = $seq_obj->description;
		if((!defined $desc) || ($desc eq '')) { $desc = 'NA'; }
		my $length = length($seq);
		
		# Base counts per entry (not returned in hash)
		my %base_counts;
		$base_counts{'A'} = ($seq =~ tr/A|a//);
		$base_counts{'C'} = ($seq =~ tr/C|c//);
		$base_counts{'T'} = ($seq =~ tr/T|t//);
		$base_counts{'G'} = ($seq =~ tr/G|g//);
		$base_counts{'N'} = ($seq =~ tr/N|n//);
		$base_counts{'U'} = ($seq =~ tr/U|u//);

		# Summary (returned in hash)
		$information{'description'}{$id}=$desc;
		$information{'full_length'} += $length;
		$information{'lengths'}{$id}=$length;
		$information{'n_of_seq'}++;
		$information{'full_length_minus_N'} += ($information{'lengths'}{$id} - $base_counts{'N'});
		foreach my $type(keys %base_counts) { $information{'base_counts'}{$type} += $base_counts{$type}; }
		push @fasta_lengths, $length;

		# Print summary per entry
		if($print_summary_for_each_entry eq 'l') {
			my $GC = ($base_counts{'G'} + $base_counts{'C'});
			my $GC_percent = sprintf("%.3f", (($GC / $length)*100));
			print "$id\t$desc\t$length\t$base_counts{'A'}\t$base_counts{'C'}\t$base_counts{'T'}\t$base_counts{'G'}\t$base_counts{'N'}\t$base_counts{'U'}\t$GC_percent\n";
		}
		if($print_summary_for_each_entry eq 's') { print "$id\t$length\n"; }
	}
	my $n50_n90 = &calculate_n50_and_n90(\@fasta_lengths);
	%information = (%information, %{$n50_n90});

	# Print summary of whole file
	if($print_summary_of_full_file ne 'n') {
		warn "\nSummary of $input:\n";
		warn "Number of sequences counted:\t$information{'n_of_seq'}\n";
		warn "Length of sequences total:\t$information{'full_length'}\n";
		warn "Length of sequences minus Ns:\t$information{'full_length_minus_N'}\n";
		foreach my $type(keys %{$information{'base_counts'}}) {
			warn "Number of $type:\t$information{'base_counts'}{$type}\n";
		}
		warn "NMAX:\t$information{'NMAX'}\n";
		warn "N50:\t$information{'N50'}\n";
		if(defined $information{'N90'}) { warn "N90:\t$information{'N90'}\n"; }
	}
	return (\%information);
}

sub calculate_n50_and_n90 {
	my $lengths = $_[0];
	my %information;

	# Init
	$information{'NMAX'} = 0;
	$information{'N50'} = 0;
	$information{'N90'} = 0;

	# Calculate N50, N90 and NMAX
	my $seq_length = 0;
	my @fasta_lengths=sort{$b<=>$a} @{$lengths}; 
	foreach(@fasta_lengths) { $seq_length += $_; }
	my ($count,$half)=(0,0);
	COUNT: for(my $i=0; $i<@fasta_lengths; $i++) {
		# Top entry
		if($count eq 0) { 
			$information{'NMAX'} = $fasta_lengths[$i]; 
		}
		$count += $fasta_lengths[$i];
		# Midway
		if(($count >= ($seq_length / 2)) && ($half == 0)) {
			$information{'N50'} = $fasta_lengths[$i];
			$half=$fasta_lengths[$i];
		}
		# 90% through
		elsif($count >= ($seq_length * 0.9)) {
			$information{'N90'} = $fasta_lengths[$i];
			last COUNT;
		}
	}
	return (\%information);
}

sub reverse_compliment { 
	my $sequence = $_[0];
	$sequence = uc ($sequence);

	# Standard
	$sequence =~ tr/[A,C,T,G]/[T,G,A,C]/;

	# IUPAC Ambiguity Codes
	$sequence =~ tr/[Y,R,K,M]/[R,Y,M,K]/;

	my @characters = split(//, $sequence);
	my $reverse_compliment = '';
	for(my $i=scalar(@characters); $i > 0; $i--) { $reverse_compliment .=  $characters[($i - 1)]; }
	return $reverse_compliment;	
}

############# Local subroutines #################
sub print_nexus {
	my ($fasta_hash, $sequence_count, $sequence_length, $data_type, $ofh_optional_parameter) = @_;

	# Make header
	my $nexus_header = "#NEXUS\n";
	$nexus_header .= "begin data;\n";
	$nexus_header .= "dimensions ntax=$sequence_count nchar=$sequence_length;\n";
	$nexus_header .= "format datatype=$data_type interleave=no gap=-;\n";
	$nexus_header .= "matrix\n";

	# Print header
	if(defined $ofh_optional_parameter) { print $ofh_optional_parameter $nexus_header; }
	else { print $nexus_header; }

	# Print Seq
	FASTA: foreach my $id(keys %{$fasta_hash}) {
		my $seq = $$fasta_hash{$id};
		$id =~ s/-/_/g;
		if(defined $ofh_optional_parameter) { print $ofh_optional_parameter "$id\n$seq\n"; }
		else { print "$id\n$seq\n"; }
	}

	# Make and print nexus footer
	my $footer_of_nexus = ";\nend;\n";
	if(defined $ofh_optional_parameter) { print $ofh_optional_parameter $footer_of_nexus; }
	else { print $footer_of_nexus; }
	return 1;
}

sub fasta_struct_print_to_fasta {
	my ($fasta_struct, $id, $split_into_files, $split_count, $split_name, $ofh_optional_parameter) = @_;
	die "$id not found in FASTA structure\n" if(!defined $$fasta_struct{'seq'}{$id});
	my $seq = $$fasta_struct{'seq'}{$id};
	$seq =~ s/(\S{60})/$1\n/g;
	my $desc = "";
	if(defined $$fasta_struct{'desc'}{$id}) { $desc = " $$fasta_struct{'desc'}{$id}"; }
	my $filename = $$fasta_struct{'filename'};
	my $entry = ">$id$desc\n$seq\n";
	my $output_name;

	# Print to standard output
	if($split_into_files eq 'n') { 
		if(defined $ofh_optional_parameter) { print $ofh_optional_parameter "$entry"; }
		else { print "$entry"; }
	}

	# Print to separate files
	else {
		my $ofh2;
		if($split_count eq 0) {
			if($split_into_files eq 1) { $output_name = ($filename . '-' . $id . '.fasta'); }
			else { $output_name = ($filename . '-split-into-' . $split_into_files . '-entries-per-file-pt-' . $split_name . '.fasta'); }
			open $ofh2, '>', $output_name or die "Cannot open $output_name: $!\n";
		}
		print $ofh2 $entry;
		$split_count++;
		if($split_count eq $split_into_files) {
			$split_count = 0;
			close $ofh2;
			$split_name++;
		}
	}
	return ($split_count, $split_name);
}

sub fasta_id_to_order_array {
	my $input = $_[0];
	my @order;
	warn "fasta_id_to_order_array: saving order from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		push @order, $id;
	}
	return (\@order);
}

1;
