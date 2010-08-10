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
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Getopt::Long;

my $AMPHORA_home = "AMPHORA_home_dir";

my $usage = qq~
Usage:  $0 <options>

Assume the Phylogenetic Marker Database is located at $AMPHORA_home

Options:
	-Partial:	query sequences contain partial genes
	-Trim:		trim the alignment using mask embedded with the marker database
	-Strict:	use a conservative mask
	-Directory	the file directory where sequences are located. Default: current directory
	-Help		print help 
~;

my ($partial, $trim, $help) = undef;
my $strict = undef;
my $dir = ".";

my (%markerlist, %imprint, %mask, @seq) = ();

GetOptions ('Partial'=>\$partial,
			'Trim'=>\$trim,
			'Strict'=>\$strict,
			'Directory=s'=>\$dir,
			'Help'=>\$help) || die "Invalid command line options\n";

die $usage if $help;
die $usage unless ("-e $AMPHORA_home/Marker/marker.list");

get_marker_list();

for my $marker (keys %markerlist) {
	next unless (-e "$marker.pep");
	print STDERR "Aligning $marker ...\n";
	(%imprint, %mask, @seq) = ();
	read_alignment($marker);
	align($marker);
	trim() if $trim;
	output($marker);
	clean();
}



####################################################################################################
sub get_marker_list {
	open (IN, "$AMPHORA_home/Marker/marker.list") || die "Can't open $AMPHORA_home/Marker/marker.list";
	while (<IN>) {
		chop;
		/^(\S+)/;
		$markerlist{$1} = 1;
	}
	close IN;
}
			
sub read_alignment {
	my $marker = shift;
	my ($ID, $seqstart, $name) = undef;
	my %seq = ();
	my $i = 1;

	open (IN, "$AMPHORA_home/Marker/$marker.gde") || die "Can't open $AMPHORA_home/Marker/$marker.gde\n";
	while (<IN>) {
		chop;
		s/"//g;	
		next if /^sequence-ID/i;	
		if (/^name\s+(\S+)/) {
			$name = $1;
			$ID = "HMM_$i";
			$i++;
		}
		if (/^type\s+MASK/) {
			$ID = $name;
		}
		if (/^offset\s+(\d+)/) {
			$seq{$ID} = "-" x $1;
		}
		if (/^sequence\s*(\S*)/) {
			$seqstart=1;
			$seq{$ID} .= $1;
			next;
		}
		if (/^\}/) {
			$seqstart = undef;
		}
		if ($seqstart) {
			s/\s+//;
			s/\./-/g;
			$seq{$ID} .= $_;
		}
	}
	close IN;
	
	if ($strict) {
		$seq{'__mask__'} = $seq{'MASK_strict'};
	}
	else {
		$seq{'__mask__'} = $seq{'MASK_relax'};
	}
	delete $seq{'MASK_strict'};
	delete $seq{'MASK_relax'};
	
	my $seqout = new Bio::SeqIO('-file'=>">$$.seed.fas", '-format'=>'fasta');
	while (my ($key, $value)  = each %seq) {
		next if ($key eq '__mask__');
		$value .= ('-' x (length($seq{'__mask__'})-length($value)));
		my $seq =new Bio::Seq('-seq'=>$value, '-id'=>$key);
		$seqout->write_seq($seq);
	}

	while (my ($key,$value) = each %seq) {
		my $j = 0;
		for ($i=0; $i<length($value); $i++) {
			next if (substr($value, $i, 1) eq '-');
            if ($i < length($seq{'__mask__'})) {
			    $imprint{$key}{$j} = substr($seq{'__mask__'}, $i, 1);
		    }
		    else {
		        $imprint{$key}{$j} = '0';
		    }
			$j++;
		}
	}
}

sub align {
	my $marker = shift;
	if ($partial) {
		system ("HMMER_path/hmmbuild -F -s $$.seed.hmm $$.seed.fas>/dev/null");
	}
	else {
		system ("HMMER_path/hmmbuild -F $$.seed.hmm $$.seed.fas>/dev/null");
	}
	system("HMMER_path/hmmalign -q -o $$.slx --mapali $$.seed.fas $$.seed.hmm $dir/$marker.pep>/dev/null"); 	

	my $in = new Bio::AlignIO('-file'=>"$$.slx");
	my $alignment = $in->next_aln();
	@seq = $alignment->each_seq();
	
	for my $seq (@seq) {
		my ($ID) = ($seq->id() =~ /^(\S+)/);
		$seq->id($ID);
		my $sequence = $seq->seq();
		$sequence =~ s/\./-/g;
		$seq->seq($sequence);
		next unless ($ID =~ /^HMM_/);
		my $j=0;
		for (my $i=0; $i<$seq->length(); $i++) {
			next if (substr($sequence, $i, 1) eq '-');
			$mask{$i} = $imprint{$ID}{$j} if (exists $imprint{$ID}{$j});
			$j++;
		}
	}
}

sub trim {
	for my $seq (@seq) {
		my $i = 0;
		my $sequence = '';
		for (split //,$seq->seq()) {
			$sequence .= $_ if ( (exists $mask{$i}) and ($mask{$i} >= 1) );
			$i++;
		}
		$seq->seq($sequence);
	}
}	

sub output {
	my $marker = shift;
	my $seqout = new Bio::SeqIO('-file'=>">$dir/$marker.aln",'-format'=>'fasta');
	for my $seq (@seq) {
		$seqout->write_seq($seq) unless ($seq->id =~ /^HMM_/);;
	}
}

sub clean {
	system ("rm $$.*");
}
 


