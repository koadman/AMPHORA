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
use Bio::Tree::Tree;
use Getopt::Long;
use Bio::Taxon;
use Bio::TreeIO;


my $usage = qq~
Usage:
	$0 <options> tree-file
	
Options:
	Cutoff
	Unrooted
Contact			Martin Wu	mwu\@tigr.org	x7830
~;

my ($taxon, $unrooted, $statistics) = undef;
my $cutoff = 70;

my $AMPHORA_home = "AMPHORA_home_dir";
my $tax_dir = "$AMPHORA_home/Taxonomy";

GetOptions ('Cutoff=s'=>\$cutoff,'Unrooted'=>\$unrooted) || die $usage;

my $treefile = $ARGV[0];
die $usage unless $treefile;

my (%taxonid, %mean, %std_deviation) =();

get_contig_taxonid();

my $taxdb = Bio::DB::Taxonomy->new(	-source   => 'flatfile',
									-nodesfile=> "$tax_dir/nodes.dmp",
									-namesfile=> "$tax_dir/names.dmp",
									-directory=>$tax_dir);

my $tree_functions = new Bio::Tree::Tree();  

my %index = ('species', 1, 'genus',2, 'family', 3, 'order', 4, 'class', 5, 'phylum', 6, 'superkingdom',7);
my %rank = (1,'species',2,'genus',3,'family',4,'order', 5, 'class', 6, 'phylum', 7, 'superkingdom');


open (IN, "$tax_dir/ladder.statistics") || die "Can't open $tax_dir/ladder.statistics";
while (<IN>) {
	chop;
	my ($rank1, $rank2, $mean, $median, $std_deviation) = split /\t/;
	if ($index{$rank1} and $index{$rank2}) {
		$mean{$rank1}{$rank2} = $mean;
		$std_deviation{$rank1}{$rank2} = $std_deviation;
	}
}

my $treeio = new Bio::TreeIO ('-file'=>$treefile);
my $tree = $treeio->next_tree();
$tree = root_tree($tree) if $unrooted;
my ($query_id) = ($treefile =~ /([^\.]+)\.tre/);
my ($outgroup,$outgroup_bootstrap, $lca, $bootstrap) = ();
my $query =$tree->findnode_by_id($query_id);
die "$query_id\tNo query in the tree" unless $query;

my ($node, $dist, $topology_assign, $dist_assign) = ();		# toplogy_assign: taxonomy assignment using only topology
if ($node = $query->ancestor()) {							# dist_assign: taxonomy assignment using topology and distance information
	$outgroup = get_taxonomy($node, $query_id);	
	($outgroup_bootstrap) = ($node->to_string =~ /(\d+)/);
	my $query_dist = get_dist($node);
	$dist += $query_dist;
	LOOP:while (my $ancestor = $node->ancestor()) {
		my $branch_length = $node->branch_length;
		$branch_length = 0 if ($branch_length < 0);
		$dist += $branch_length;
		($bootstrap) = ($ancestor->to_string =~ /(\d+)/);
		if ($bootstrap >= $cutoff) {		
			$lca = get_taxonomy($ancestor, $query_id);
			$topology_assign = $lca unless $topology_assign;
			if ($index{$lca->rank}) {
				my $fraction;
				if ($dist == 0) {
					$fraction = 0;
				}
				else {
					$fraction = ($dist - $query_dist)/$dist;
				}
				$dist_assign = assign($lca, $outgroup, $fraction);
				last LOOP;
			}
		}
		$node = $node->ancestor;
	}
}		

my $assign;
if ($outgroup_bootstrap < ($cutoff-20)) {
	$assign = $topology_assign;
}
elsif (get_depth($topology_assign) > get_depth($dist_assign)) {
	$assign = $topology_assign;
}
else {
	$assign = $dist_assign;
}

$assign = $taxdb->get_taxon(-taxonid => '2') unless $assign; # bacteria

print $query_id, "\t", $assign->scientific_name,"\t",$outgroup->scientific_name,"\tBootstrap:",$bootstrap,"\n";


########################################################################################################################

sub get_contig_taxonid {
	open (IN, "$AMPHORA_home/Reference/contig.taxonid") || die "Can't open $AMPHORA_home/Reference/contig.taxonid\n";
	while (<IN>) {
		my ($contigid, $taxonid) = /^(\S+)\s+(\S+)/;
		$taxonid{$contigid} = $taxonid;
	}
	close IN;
}

sub assign {
	my ($taxon, $outgroup, $fraction) = @_;
	$fraction = 0 if ($fraction < 0);
	my %prob = ();
	my $rank1 = $taxon->rank;
	for my $rank2 (keys %{$mean{$rank1}}) {
		$prob{$rank2} = abs( ($fraction - $mean{$rank1}{$rank2})/$std_deviation{$rank1}{$rank2});
	}
	my @sort = sort {$prob{$a} <=> $prob{$b}} keys %prob;
	my $tree = new Bio::Tree::Tree(-node => $outgroup);
	my $index = $index{$sort[0]};
	my $taxon = undef;
 	while (! ($taxon = $tree->find_node(-rank => $rank{$index}) ) ) {
 		$index ++;
 	}
	return $taxon;
}	
		
sub root_tree {
	my $tree = shift;
	my @nodes;	# root the tree using Thermus thermophilus and Aquifex aeolicus as outgroup
	for my $node ($tree->get_leaf_nodes) {
		my $id = $node->id();
		$id =~ s/'//g;
		$id =~ s/\/.+$//;
		$node->id($id);
		if ($id =~ /NC_001263/) { #Deinococcus radiodurans R1 
			push @nodes, $node;
		}
		elsif ($id =~ /NC_008025/) { #Deinococcus geothermalis DSM 11300 
			push @nodes, $node;
		}
	}
	
	my $root_node = $tree->get_lca(-nodes =>\@nodes);
	$tree->reroot($root_node);
	return $tree;
}


sub get_taxonomy {
	my ($node, $query_id) = @_;
	my @taxon = ();

	for my $descend ($node->get_all_Descendents()) {
		next unless $descend->is_Leaf();
		my $id = $descend->id();	
		$id =~ s/'//g;
		$id =~ s/\/.+$//;	
		next if ($id eq $query_id);
		$id =~ /^[^-]+-(\S+)/;
		push @taxon, $taxdb->get_taxon(-taxonid=>$taxonid{$1});
	}
	if ($#taxon == 0) {
		return $taxon[0];
	}
	else {
		my $tree_functions = new Bio::Tree::Tree();
		return $tree_functions->get_lca(\@taxon);
	}
}

sub get_depth {
        my $depth = 0;
        my $node = shift;
        my @lineage = $tree_functions->get_lineage_nodes($node);
        for (@lineage) {
                $depth++;                                                      
        }
        return $depth;
}

sub get_dist {
	my $node = shift;
	if ($node->is_Leaf) {
		return 0;
	}
	else {
		my @children = $node->each_Descendent;
		return  ($children[0]->branch_length+$children[1]->branch_length+get_dist($children[0])+get_dist($children[1]))/2;
	}
}
