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
use Bio::TreeIO;
use Getopt::Long;
use Bio::AlignIO;
use Cwd;

my $usage = qq~
Description: Given a batch of bacterial sequences, this program will first identify a suite of 
phylogenetic markers from them and then assign each of them a phylotype using a tree-based method 

Usage: $0 <options> sequence-file output-file
        Options:
                -Replicates: number of bootstrap replicates
                -BootstrapCutoff: normalized to 100 replicates (1-100)%  default 70

~;

my $AMPHORA_home = "AMPHORA_home_dir";

my $bootstrap_cutoff = 70;
my $replicates = 100;
my $output = undef;

GetOptions(	'Replicates=i'=>\$replicates,
			'BootstrapCutoff=i'=>\$bootstrap_cutoff) ||die $usage;

die $usage unless ($ARGV[0] and $ARGV[1]);
my $sequence_file = $ARGV[0];
my $output = $ARGV[1];

my $outgroup;

my (%markerlist) = ();

my $working_dir = getcwd();
get_marker_list();
identify_marker();
sequence_align();

mkdir "$working_dir/Trees";
mkdir $$;
chdir $$;

for my $marker (keys %markerlist) {
	my @refseq = ();	
	$outgroup = undef;
	
	next unless (-e "$working_dir/$marker.aln");
	my $refin = new Bio::SeqIO('-file'=>"$AMPHORA_home/Reference/Sequence/$marker.aln");
	while (my $refseq = $refin->next_seq()) {		
		push @refseq, $refseq;
		if ($refseq->id =~ /NC_001263/ or $refseq->id =~ /NC_008025/) {
			$outgroup .= $refseq->id.",";
		}
	}
	$outgroup =~ s/\,$//;
	
	my $queryin = new Bio::SeqIO('-file'=>"$working_dir/$marker.aln");
	while (my $query = $queryin->next_seq()) {
		print STDERR "Assigning ".$query->id()." ...\n";
		my @seq = (@refseq, $query);
		write_sequences(@seq);	
		make_parsimony_trees($marker);		
		add_branch_length($query->id);
		assign_phylotype($query->id);
		cleanup();
	}
}

chdir $working_dir;
system("rm -r $$");

###########################################################################################	

sub get_marker_list {
	open (IN, "$AMPHORA_home/Marker/marker.list") || die "Can't open $AMPHORA_home/Marker/marker.list";
	while (<IN>) {
		chop;
		/^(\S+)/;
		$markerlist{$1} = 1;
	}
	close IN;
}

sub identify_marker {
	print STDERR "Identifying phylogenetic markers...\n";
	my $status = system ("perl $AMPHORA_home/Scripts/MarkerScanner.pl $sequence_file");
	die "No marker sequences were found\n" unless ($status == 0);
}

sub sequence_align {
	print STDERR "Aligning marker sequences...\n";
	system ("perl $AMPHORA_home/Scripts/MarkerAlignTrim.pl -Trim -Partial -Directory $working_dir");
}

sub write_sequences {
	my @seq = @_;
	
	my $seqout = new Bio::SeqIO('-file'=>">$$.fas",'-format'=>'fasta');
	for (@seq) {
		$seqout->write_seq($_);
	}
	
	my $alignin = new Bio::AlignIO('-file'=>"$$.fas");
	my $align = $alignin->next_aln();
	my $alignout = new Bio::AlignIO('-file'=>">$$.phy",'-format'=>'phylip','-idlength'=>50);
	$alignout->write_aln($align);
	$alignout = new Bio::AlignIO('-file'=>">$$.stock",'-format'=>'stockholm');
	$alignout->write_aln($align);
		system ("$AMPHORA_home/bin/seqboot $$.phy $replicates &>/dev/null");
}


sub make_parsimony_trees {
	my $marker = shift;
	system ("$AMPHORA_home/bin/raxmlHPC -f p -t $AMPHORA_home/Reference/Tree/$marker.tre -s $$.phy -m PROTGAMMAWAG -n ref &>/dev/null");
	system ("mv RAxML_parsimonyTree.ref $$.ref.tre");
	for my $i ( 1 .. $replicates) {
		system ("nice -n 15 $AMPHORA_home/bin/raxmlHPC -f p -t $AMPHORA_home/Reference/Tree/$marker.tre -s bs.$i -m PROTGAMMAWAG -n $i &>/dev/null");
	}
	
	system ("cat RAxML_parsimonyTree.* > $$.trees");

	my $treeio = new Bio::TreeIO('-file'=>"$$.ref.tre");
	my $tree = $treeio->next_tree();
	for my $node ($tree->get_nodes()) {
		$node->branch_length(0.5);
	}
	$treeio = new Bio::TreeIO('-file'=>">$$.ref.tre", '-format'=>'newick');
	$treeio->write_tree($tree);
	system ("$AMPHORA_home/bin/raxmlHPC -f b -t $$.ref.tre -z $$.trees -o $outgroup -m PROTGAMMAWAG -n $$ -s $$.phy &>/dev/null");
}

sub add_branch_length {
	my $query_id = shift;
	system ("$AMPHORA_home/bin/quicktree -out m $$.stock > $$.dist");
	system ("perl $AMPHORA_home/Scripts/NJ_CompleteConstraint.pl -Dist $$.dist -Tree RAxML_bipartitions.$$ > $query_id.tree");
}

sub assign_phylotype {
	my $query_id = shift;
	system ("perl $AMPHORA_home/Scripts/AssignPhylotype.pl -Cutoff $bootstrap_cutoff $query_id.tree >> $working_dir/$output");
}

sub cleanup {
	system ("mv *.tree $working_dir/Trees/.");
	system ("rm *");
}
