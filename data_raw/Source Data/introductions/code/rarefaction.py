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
import pandas as pd
import argparse
import ast
from collections import Counter
import datetime
import csv
import random
import copy


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


def infer_mugration(tre, aln_file, regions, tre_rep_dir, i):
	out_path = f'{tre_rep_dir}/{i}'
	tre, letter_to_state, state_to_letter = \
		wrappers.reconstruct_discrete_traits(tre, regions, 
			sampling_bias_correction=2.5)
	pd.DataFrame({'letter':list(letter_to_state.keys()), 
				  				'region': list(letter_to_state.values())}).\
								to_csv(f'{out_path}/{i}_refined_regions.csv', index=False)
	pd.DataFrame({'name': [item.name for item in tre.tree.find_clades()],
				  'state': [list(item.marginal_profile[0]) for item in tre.tree.find_clades()]}).\
					to_csv(f'{out_path}/{i}_refined_node_states.csv', index=False)
	return(tre, letter_to_state, state_to_letter)


def estimate_importations(tre, counties, state, tre_rep_dir, i):
	out_path = f'{tre_rep_dir}/{i}'
	subprocess.run(shlex.split(f'mkdir {out_path}'))
	# list of regions and their identifying character
	region_dict = pd.read_csv(f'{out_path}/{i}_refined_regions.csv')
	# Dictionary of each node's state
	node_states = pd.read_csv(f'{out_path}/{i}_refined_node_states.csv')
	node_states['state'] = \
		node_states['state'].apply(lambda k: ast.literal_eval(k))
	node_states['ML_state'] = \
		node_states['state'].apply(lambda k: region_dict.loc[k.index(max(k)), 'region'])
	node_states_dict = {item['name']: item['ML_state'] for index, item in node_states.iterrows()}
	# Dictionary of each node's time
	# formats the unique importation nodes list
	formatted_unique_importation_nodes = []
	formatted_columns = ['destination region', 'source region', 'descendant tips', 
					     'source node', '# descendants', 'midroot time', 'lower midroot time', 'upper midroot time']
	# tabulates the number of introductions into each region
	importation_dict = {}
	for region in counties:
		# List of tree tips from specific region of interest
		region_tre_tips = \
			[item for item in tre.get_terminals() if node_states_dict[item.name] == region]
		# For each tip from the region of interest, get path from root to ip
		region_paths = [tre.get_path(item) for item in region_tre_tips]
		# state of each step in each path
		region_paths_states = [[node_states_dict[item.name] for item in path] for path in region_paths]
		# node which precede the earliest node assigned to ANY of the regions of interest
		importation_node_indices = [min([i for i, item in enumerate(path) if item in counties])-1 for path in region_paths_states]
		importation_nodes = [region_paths[i][item] for i, item in enumerate(importation_node_indices)]
		importation_nodes_dict = {region_paths[i][-1]: item for i, item in enumerate(importation_nodes)}
		importation_nodes = list(set(importation_nodes))
		importation_dict[region] = len(importation_nodes)
	results = pd.DataFrame({'region': list(importation_dict.keys()), 
				  'importations': list(importation_dict.values())})
	results['n'] = i.split('_')[0]
	results['rep'] = i.split('_')[1]
	return(results)


def rarefaction(tre, aln_file, counties, state, all_regions_dict, n_samples, reps, out_path):
	tre_rep_dir = f'{out_path}/rarefaction'
	subprocess.run(shlex.split(f'mkdir {tre_rep_dir}'))
	for n in n_samples:
		for rep in range(reps):
			this_tre = copy.deepcopy(tre)
			i = f'{n}_{rep}'
			subprocess.run(shlex.split(f'mkdir {tre_rep_dir}/{i}'))
			tips_list = [item for item in this_tre.get_terminals() if all_regions_dict[item.name] in counties]
			sample_tips_list = random.sample(tips_list, k=n)
			remove_samples = [item for item in tips_list if item not in sample_tips_list]
			print(f"removing {len(remove_samples)} of {len(tips_list)} Wisconsin tips")
			_ = [this_tre.prune(item) for item in remove_samples]
			this_tre, letter_to_state, state_to_letter = \
				infer_mugration(this_tre, aln_file, all_regions_dict, tre_rep_dir, i)
			these_results = estimate_importations(this_tre.tree, counties, state, tre_rep_dir, i) 
			these_results.to_csv(f'{tre_rep_dir}/{i}/{i}_rarefaction.csv', 
							   index=None)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--tree_file', 
		help='tree to run analyses on',
		default='results/ml_refined_time.newick')
	parser.add_argument('--outgroup',
				help='Outgroup to use in tree',
				default='Wuhan/Hu-1/2019')
	parser.add_argument('--metadata',
		default='data/metadata_adjusted.tsv')
	parser.add_argument('--wi_file',
		default='data/MKEvsDane.csv')
	parser.add_argument('--counties', help='name of counties to analyze', 
		nargs='+', default=('MHD', 'UWHC'))
	parser.add_argument('--state', 
		help='state to analyze, encompasses counties', default='Wisconsin')
	parser.add_argument('--aln',
		default='datasubsampled_alignment_neighbors.fasta')
	parser.add_argument('--outdir', default='results')
	parser.add_argument('--n_samples',
		help='numer of samples from region of interest', nargs='+')
	parser.add_argument('--reps',
		help='number of replicates to do for each n', default=10)
	args = parser.parse_args()
	n_samples = [int(item) for item in args.n_samples]
	print(n_samples)
	all_regions_dict = make_regions_dict(args.wi_file, args.metadata)
	tre = Phylo.read(args.tree_file, 'newick')
	tre.root_with_outgroup(args.outgroup)	# operates in place
	counties = ('MHD', 'UWHC', 'Wisconsin')
	rarefaction(tre, args.aln, args.counties, args.state, all_regions_dict, n_samples, args.reps, args.outdir)


main()
