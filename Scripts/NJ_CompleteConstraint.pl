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
use Getopt::Long;
use lib "bioperl_home_dir";
use Bio::TreeIO;


my $usage = qq~
Given a distance matrix and a tree topology, this program infers the tree branch lengths using a 
constrained NJ algorithm

Usage: $0 <options>
        Options:
                -Distance matrix file: in phylip format 
                -Tree: 
~;

my ($distfile, $treefile) = ();

GetOptions('Distance=s'=>\$distfile, 'Tree=s'=>\$treefile)  || die $usage;
die $usage unless ($distfile and $treefile);

my (%dist, %internal_id, %weight, %tips) = ();

my $treeio = new Bio::TreeIO('-file'=>$treefile);
my $tree = $treeio->next_tree();
for my $leaf ($tree->get_leaf_nodes) {
	my $id = $leaf->internal_id;
	$internal_id{$leaf->id} = $id;
	$tips{$id} = 1;
}

read_distance_file($distfile);

my $node = $tree->get_root_node();
my @children = $node->each_Descendent;
climb($children[0], $children[1]);

$treeio = new Bio::TreeIO('-fh'=>*STDOUT,'-format'=>'newick');
$treeio->write_tree($tree);

sub read_distance_file {
	my $file = shift;
	my $i = 0;
	my @id = ();
	open (IN, $file) || die "Can't open $file\n";
	while (<IN>) {
		chop;
		next if (/^\s+\d/);
		my @l = split;
		$l[0] =~ s/\/.+//g;
		$id[$i] = $l[0];
		for (my $j=0; $j< $i; $j++) {
			my $id1 = $internal_id{$id[$j]};
			my $id2 = $internal_id{$id[$i]};
			$dist{$id1}{$id2} = $l[$j+1];
			$dist{$id2}{$id1} = $l[$j+1];
		}
		$i++;
	}
	close IN;
}

sub climb {
	my ($node1, $node2)  = @_;
	my $n1 = $node1->internal_id;
	my $n2 = $node2->internal_id;
	if (! exists $tips{$n1}) {		
		my @children = $node1->each_Descendent;
		climb($children[0], $children[1]);
	}
	if (! exists $tips{$n2}) {		
		my @children = $node2->each_Descendent;
		climb($children[0], $children[1]);
	}
	if ((exists $tips{$n1}) and (exists $tips{$n2})) {
		nj($node1, $node2);
		return;
	}
}

sub nj {
	my ($node1, $node2) = @_;
	my $n1 = $node1->internal_id;
	my $n2 = $node2->internal_id;

	if ((scalar keys %tips) == 2) {
		$node1->branch_length($dist{$n1}{$n2}/2);
		$node2->branch_length($dist{$n1}{$n2}/2);
		return;
	}

	for my $tip (keys %tips) {
		$weight{$n1} += $dist{$n1}{$tip};
	}

	$weight{$n1} /= ((scalar keys %tips) - 2);	
	for my $tip (keys %tips) {
		$weight{$n2} += $dist{$n2}{$tip};
	}	
	$weight{$n2} /= ((scalar keys %tips) - 2);	
	$node1->branch_length( ($dist{$n1}{$n2}+$weight{$n1}-$weight{$n2})/2);
	$node2->branch_length( ($dist{$n1}{$n2}+$weight{$n2}-$weight{$n1})/2);

	delete $tips{$n1};
	delete $tips{$n2};
	my $joint = $node1->ancestor->internal_id;		
	for my $tip (keys %tips) {
		$dist{$joint}{$tip} =($dist{$n1}{$tip}+$dist{$n2}{$tip}-$dist{$n1}{$n2})/2;
		$dist{$tip}{$joint} = $dist{$joint}{$tip};
	}
	$tips{$joint} = 1;	
}

