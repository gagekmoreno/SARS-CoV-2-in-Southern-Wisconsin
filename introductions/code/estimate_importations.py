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
from treetime.utils import numeric_date
import pandas as pd
import argparse
import ast
import numpy as np


def label_nodes(tre):
	node_counter=0
	for item in tre.get_nonterminals():
		if item.name is None:
			item.name = 'NODE_'+format(node_counter, '07d')
			node_counter+=1
	return(tre)


def infer_timetree(tre, aln_file, dates, tre_rep_dir, i, resolve, outgroup):
	out_path = f'{tre_rep_dir}/{i}'
	subprocess.run(shlex.split(f'mkdir {out_path}'))
	# writing as XML because newick writer throws recursive error
	tre_file = f'{out_path}/{i}.xml'
	with open(tre_file, 'w') as out:
		Phylo.write(tre, out, 'phyloxml')
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
    # TODO change max iter to 2 in final version, this only for testing
    # TODO add vary rate parameter
	print(f'resolve polytomies: {resolve}')
	tre.run(root=outgroup, infer_gtr=True, max_iter=2,
                    branch_length_mode='auto', resolve_polytomies=resolve,
                    time_marginal='assign', fixed_clock_rate=0.0008, vary_rate=0.0004)
	times = pd.DataFrame({'name': [item.name for item in tre.tree.find_clades()],
						  'date': [item.date for item in tre.tree.find_clades()],
						  'lower': [list(tre.get_max_posterior_region(item, 0.9))[0] 
						  			for item in tre.tree.find_clades()],
						  'upper': [list(tre.get_max_posterior_region(item, 0.9))[1] 
						  			for item in tre.tree.find_clades()]},
						  index = range(0, len([item for item in tre.tree.find_clades()])))
	times['date'] = pd.to_datetime(times['date']).apply(numeric_date)
	times.to_csv(f'{out_path}/{i}_refined_node_times.csv')
	with open(f'{out_path}/{i}_refined.xml', 'w') as out_file:
		Phylo.write(tre.tree, out_file, 'phyloxml')
	tre.branch_length_to_years()
	with open(f'{out_path}/{i}_refined_time.xml', 'w') as out_file:
		Phylo.write(tre.tree, out_file, 'phyloxml')
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


def estimate_importations(tre, these_regions, tre_rep_dir, i):
	out_path = f'{tre_rep_dir}/{i}'
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
	node_times = pd.read_csv(f'{out_path}/{i}_refined_node_times.csv')
	node_times_dict = {item['name']: (item['lower'], item['date'], item['upper']) for index, item in node_times.iterrows()}
	# formats the unique importation nodes list
	formatted_unique_importation_nodes = []
	formatted_columns = ['destination region', 'source region', 'descendant tips', 
					     'source node', '# descendants', 'midroot time', 'lower midroot time', 'upper midroot time']
	# tabulates the number of introductions into each region
	n_introductions_regions = {}
	for region in these_regions:
		# List of tree tips from specific region of interest
		region_tre_tips = \
			[item for item in tre.get_terminals() if node_states_dict[item.name] == region]
		# For each tip from the region of interest, get path from root to ip
		region_paths = [tre.get_path(item) for item in region_tre_tips]
		# state of each step in each path
		region_paths_states = [[node_states_dict[item.name] for item in path] for path in region_paths]
		# node which precede the earliest node assigned to ANY of the regions of interest
		importation_node_indices = [min([i for i, item in enumerate(path) if item in these_regions])-1 for path in region_paths_states]
		importation_nodes = [region_paths[i][item] for i, item in enumerate(importation_node_indices)]
		importation_nodes_dict = {region_paths[i][-1]: item for i, item in enumerate(importation_nodes)}
		importation_nodes = list(set(importation_nodes))
		n_introductions_regions[region] = len(importation_nodes)
		for importation_node in importation_nodes:
			node_time = node_times_dict[importation_node.name]
			node_descendants = [key.name for key, value in importation_nodes_dict.items() 
								if value == importation_node]
			# sorts by teh median time
			min_descendant_time = sorted([node_times_dict[item] for item in node_descendants], 
										  key=lambda k: k[1])[0]
			this_importation_node = \
					[region, 
					 node_states_dict[importation_node.name], 
					 node_descendants, 
					 importation_node.name,
					 len(node_descendants),
					 (min_descendant_time[1] - node_time[1])/2 + node_time[1],
					 (min_descendant_time[0] - node_time[0])/2 + node_time[0],
					 (min_descendant_time[2] - node_time[2])/2 + node_time[2]]
			# adds output row to larger list
			formatted_unique_importation_nodes.append(this_importation_node)
	# turns formatted list into data frame
	formatted_unique_importation_nodes = \
		pd.DataFrame(formatted_unique_importation_nodes, columns=formatted_columns)
		# saves formatted list
	formatted_unique_importation_nodes.to_csv(f'{out_path}/{i}_refined_importations.csv')
	len(formatted_unique_importation_nodes['source node'])
	n_introductions = len(set(formatted_unique_importation_nodes['source node']))
	print(f'Estimated {n_introductions} introductions into Wisconsin')
	print(f'Of these:')
	for key, value in n_introductions_regions.items():
		print(f'{value} led to tips from {key}')
	return(formatted_unique_importation_nodes)


def format_seq_name(seq_name):
	return(seq_name.replace('hCoV-19/', '').split('|')[0].replace(' ', ''))


# save gtr from states
def main():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--bootstrap_replicate', 
						help='bootstrap replicate to analyze',
						type=int,
						default=0)
	parser.add_argument('--tree_file', 
						help='file which contains set of trees to analyze',
						default='results/subsampled_alignment_neighbors.fasta.ufboot')
	parser.add_argument('--metadata_file',
						help='metadata file',
						default='data/metadata_adjusted.tsv')
	parser.add_argument('--mke_dane_file',
					help='file with list of sequences from Dane and MKE',
					default='data/MKEvsDane.csv')
	parser.add_argument('--aln_file', 
						help='alignment with nearest neighbors', 
						default='data/subsampled_alignment_neighbors.fasta')
	parser.add_argument('--outgroup',
					help='Outgroup to use in tree',
					default='Wuhan/Hu-1/2019')
	parser.add_argument('--resolve', dest='resolve', action='store_true')
	parser.add_argument('--no-resolve', dest='resolve', action='store_false')
	parser.set_defaults(resolve=False)
	args = parser.parse_args()
	dates = parse_dates(args.metadata_file)
	bootstrap_tres = list(Phylo.parse(args.tree_file, 'newick'))
	tre_rep_dir = f'{args.tree_file}_tres'
	print(tre_rep_dir)
	subprocess.run(shlex.split(f'mkdir {tre_rep_dir}'))
	wi_tips = pd.read_csv(args.mke_dane_file)
	wi_tips = {item['Sample Name']: item['wisconsin_catchment_area'] for i, item in wi_tips.iterrows()}
	all_regions_dict = make_regions_dict(args.mke_dane_file, args.metadata_file)
	# Roots bootstrap trees
	i = args.bootstrap_replicate
	print(f'replicate {i}')
	bootstrap_tres[i].root_with_outgroup(args.outgroup)		# operates in place
	print(f'replicate {i} tree rooted on {args.outgroup}')
	bootstrap_tre = label_nodes(bootstrap_tres[i])
	print(f'replicate {i} tree internal nodes labelled')
	timetree_tre = infer_timetree(bootstrap_tre, args.aln_file, dates, tre_rep_dir, i, args.resolve, args.outgroup)
	print(f'replicate {i} timetree inferred')
	tre, letter_to_state, state_to_letter = infer_mugration(timetree_tre.tree, args.aln_file, all_regions_dict, tre_rep_dir, i)
	print(f'replicate {i} mugrations inferred')
	wi_regions = ('Wisconsin', 'UWHC', 'MHD')
	importations = estimate_importations(tre.tree, wi_regions, tre_rep_dir, i) 
	print(f'replicate {i} importations estimated')


if __name__ == "__main__":
	main()


