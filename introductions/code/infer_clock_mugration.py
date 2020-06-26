#!/usr/bin/env python3
import subprocess
import shlex
from Bio import Phylo
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import treetime as tt
from treetime import TreeTime
from treetime import wrappers
from treetime.utils import parse_dates
import pandas as pd
import argparse
import ast
import csv


def label_nodes(tre):
	node_counter=0
	for item in tre.get_nonterminals():
		item.confidence = None
		if item.name is None:
			item.name = 'NODE_'+format(node_counter, '07d')
			node_counter+=1
	return(tre)


def make_regions_dict(wi_file, metadata_file):
	wi_regions = pd.read_csv(wi_file)
	wi_tips= {item['Sample Name']: item['wisconsin_catchment_area'] for 
					index, item in wi_regions.iterrows()}
	metadata = pd.read_csv(metadata_file, sep='\t')
	usa_tips = {item['strain']: item['division'] for index, item in metadata[metadata['country'] == 'USA'].iterrows()}
	global_tips = {item['strain']: item['country'] for index, item in metadata[metadata['country'] != 'USA'].iterrows()}
	regions = wi_tips.copy()
	included_keys = set(wi_tips.keys())
	for key, value in usa_tips.items():
		if key not in included_keys:
			regions[key] = value
	included_keys.add(item for item in set(usa_tips.items()))
	for key, value in global_tips.items():
		if key not in included_keys:
			regions[key] = value
	return(regions)


def infer_timetree(tre, aln_file, dates, outgroup, resolve, out_path):
	tre_path = f'{out_path}/ml.newick'
	with open(tre_path, 'w') as out_file:
		Phylo.write(tre, out_file, 'newick', format_branch_length='%1.8f')
	tre = TreeTime(gtr='JC69', tree=tre, precision='auto',
	                 aln=aln_file, verbose=2, dates=dates)
	# from NextStrain Augur refine.py
	# treetime clock filter will mark, but not remove bad tips
	tre.clock_filter(reroot=outgroup, n_iqd=4, plot=False) 
	# remove them explicitly
	leaves = [x for x in tre.tree.get_terminals()]
	for n in leaves:
		if n.bad_branch:
			tre.tree.prune(n)
			print('pruning leaf ', n.name)
	# fix treetime set-up for new tree topology
	tre.prepare_tree()
	print(f'resolve polytomies: {resolve}')
	tre.run(root=outgroup, infer_gtr=True, max_iter=2,
	                branch_length_mode='auto', resolve_polytomies=resolve,
	                time_marginal='assign', vary_rate=0.0004, fixed_clock_rate=0.0008)
	times = pd.DataFrame({'name': [item.name for item in tre.tree.find_clades()],
						  'date': [item.date for item in tre.tree.find_clades()],
						  'lower': [list(tre.get_max_posterior_region(item))[0] for item in tre.tree.find_clades()],
						  'upper': [list(tre.get_max_posterior_region(item))[1] for item in tre.tree.find_clades()]}, 
						  index = range(0, len([item for item in tre.tree.find_clades()])))
	times.to_csv(tre_path.replace('.newick', '_refined_node_times.csv'))
	# Saves refined tree
	with open(tre_path.replace('.newick', '_refined.newick'), 'w') as out_file:
		Phylo.write(tre.tree, out_file, 'newick', format_branch_length='%1.8f')
	tre.branch_length_to_years()
	with open(tre_path.replace('.newick', '_refined_time.newick'), 'w') as out_file:
		Phylo.write(tre.tree, out_file, 'newick', format_branch_length='%1.8f')
	return(tre)


def infer_mugration(tre, aln_file, regions, out_path):
	tre, letter_to_state, state_to_letter = \
		wrappers.reconstruct_discrete_traits(tre, regions, 
			sampling_bias_correction=2.5)
	pd.DataFrame({'letter':list(letter_to_state.keys()), 
				  				'region': list(letter_to_state.values())}).\
								to_csv(f'{out_path}/ml_refined_regions.csv', index=False)
	pd.DataFrame({'name': [item.name for item in tre.tree.find_clades()],
				  'state': [list(item.marginal_profile[0]) for item in tre.tree.find_clades()]}).\
					to_csv(f'{out_path}/ml_refined_node_states.csv', index=False)
	return(tre, letter_to_state, state_to_letter)


def main():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--out_dir', 
						help='output directory', 
						default='results')
	parser.add_argument('--aln_file', 
						help='alignment with nearest neighbors', 
						default='data/subsampled_alignment_neighbors.fasta')
	parser.add_argument('--tree_file', 
						help='tree from which to estimate clock rate',
						default='results/subsampled_alignment_neighbors.fasta.treefile')
	parser.add_argument('--metadata_file',
						help='metadata file',
						default='data/metadata_adjusted.tsv')
	parser.add_argument('--mke_dane_file',
					help='file with list of sequences from Dane and MKE',
					default='data/MKEvsDane.csv')
	parser.add_argument('--outgroup',
					help='Outgroup to use in tree',
					default='Wuhan/Hu-1/2019')
	parser.add_argument('--resolve', dest='resolve', action='store_true')
	parser.add_argument('--no-resolve', dest='resolve', action='store_false')
	parser.set_defaults(resolve=False)
	args = parser.parse_args()
	subprocess.run(shlex.split(f'mkdir {args.out_dir}'))
	regions_dict = make_regions_dict(args.mke_dane_file, args.metadata_file)
	tre = Phylo.read(args.tree_file, 'newick')
	tre.root_with_outgroup(args.outgroup)		# operates in place
	tre = label_nodes(tre)	# labels nodes of tree
	dates = parse_dates(args.metadata_file)	# parses tip dates
	timetree_tre = infer_timetree(tre, args.aln_file, dates, args.outgroup, args.resolve, args.out_dir)
	# infers mugration
	tre, letter_to_state, state_to_letter = \
		infer_mugration(timetree_tre.tree, args.aln_file, regions_dict, args.out_dir)


if __name__ == "__main__":
	main()

