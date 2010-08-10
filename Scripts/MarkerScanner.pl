# AMPHORA (version 1.0) An Automated Phylogenomic Inference Pipeline for bacterial sequences. 
# Copyright 2008 by Martin Wu
 
# This file is part of AMPHORA.

# AMPHORA is free software: you may redistribute it and/or modify its under the terms of the 
# GNU General Public License as published by the Free Software Foundation; either version 2 of
# the License, or any later version.

# AMPHORA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details (http://www.gnu.org/licenses/).
 

# For any other inquiries send an Email to Martin Wu
#       mw4yv@virginia.edu
 
# When publishing work that is based on the results from AMPHORA please cite:
# Wu M and Eisen JA. A simple, fast, and accurate method of phylogenomic inference. Genome Biol 2008 Oct 13; 9(10) R151i

use strict;
use lib "bioperl_home_dir";
use Bio::SeqIO;
use Bio::SearchIO;

my (%markerlist, %seq, %Hmmlength) = ();

my $usage = qq~
Usage: $0 fasta-file
~;

my $AMPHORA_home = "AMPHORA_home_dir";
die $usage unless $ARGV[0];

get_marker_list();
read_marker_hmms();

system ("cp $ARGV[0] $$.query");

# Blastp search
system ("BLASTP_path/setdb $$.query >&/dev/null");
system("BLASTP_path/blastp $$.query $AMPHORA_home/Marker/markers.fas -E=0.1 -V=50000 -B=50000 -seqtest -Z=5000> $$.blastp");
get_blast_hits();

# HMM search
system ("HMMER_path/hmmpfam -Z 5000 -E 1e-3 $AMPHORA_home/Marker/markers.swhmm $$.candidate > $$.hmmsearch");		# fix the number of sequences in the database for E-value calculation
get_hmm_hits();

# clean up
system ("rm $$.*");

####################################################################################################################
sub get_marker_list {
	open (IN, "$AMPHORA_home/Marker/marker.list") || die "Can't open $AMPHORA_home/Marker/marker.list";
	while (<IN>) {
		chop;
		/^(\S+)/;
		$markerlist{$1} = 1;
	}
	close IN;
}

sub read_marker_hmms {
	open (IN, "$AMPHORA_home/Marker/markers.swhmm") || die "Can't open $AMPHORA_home/Marker/markers.swhmm\n";
	my $name;
	while (<IN>) {
        chop;
		if (/^NAME\s+(\S+)/) {
			$name = $1;
		}
		elsif (/^LENG\s+(\d+)/) {
			$Hmmlength{$name} = $1;
		}
	}
	close IN;
}

sub get_blast_hits {
	my %hits = ();

	my $in = new Bio::SearchIO('-file' => "$$.blastp");
	while (my $result = $in->next_result) {	
		while( my $hit = $result->next_hit ) {
			$hits{$hit->name()} = 1;
		}
	}
	unless (%hits) {
		system("rm $$.*");
		exit(1) 
	}
	my $seqin = new Bio::SeqIO('-file'=>$ARGV[0]);
	my $seqout = new Bio::SeqIO('-file'=>">$$.candidate",'-format'=>'fasta');

	while (my $seq = $seqin->next_seq) {
		if ($hits{$seq->id}) {
			$seq{$seq->id} = $seq;
			$seqout->write_seq($seq);
		}
	}
}

sub get_hmm_hits {
	my $hmmsearch = new Bio::SearchIO ('-file'=>"$$.hmmsearch", '-format' => 'hmmer');
	my $marker;
	my %hits =();

	RESULT:while ( my $result = $hmmsearch->next_result ){
		my $query = $result->query_name();
		my ($query_match, $hit_match, $best_evalue, $hitname) = ();
		while (my $hit = $result->next_hit() ){
			if ( !(defined $best_evalue) or $best_evalue > $hit->significance) {		
				$best_evalue = $hit->significance;
				$hitname = $hit->name;
	
				$query_match = 0;
				$hit_match = 0;
				while (my $hsp = $hit->next_hsp) {
					$query_match += ($hsp->end('query')-$hsp->start('query') );
					$hit_match += ($hsp->end('hit')-$hsp->start('hit')); 
				}
				$query_match  /= $seq{$query}->length;
				$hit_match /= $Hmmlength{$hit->name()};
			}	
		}
		next RESULT unless $markerlist{$hitname};
		next RESULT if ( ($query_match < 0.7) and ($hit_match < 0.7) );		#ignore the hit if the match is partial
	
		$hits{$hitname}{$query} = 1;	
	}
	unless (%hits) {
		system("rm $$.*");
		exit(1) 
	}
	for my $marker (keys %hits) {
		my $seqout = new Bio::SeqIO('-file'=>">$marker.pep",'-format'=>'fasta');
		for my $seqid (keys %{$hits{$marker}}) {
			$seqout->write_seq($seq{$seqid});
		}
	}
}
