#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
from calendar import isleap
import numpy as np
import string
import seaborn as sns
import baltic as bt
from treetime.utils import parse_dates
from treetime.utils import numeric_date
import ast


def set_style():
	# Sets plotting style
	# Defines colors based on Nord (https://www.nordtheme.com)
	global greys
	greys = ["#2E3440", "#3b4252", "#434C5E", "#4C566A"]
	global mpl
	# MatPlotLib style
	mpl.rcParams['font.family'] = 'sans-serif'
	#mpl.rcParams['font.sans-serif'] = 'Helvetica'
	#mpl.rcParams['font.weight'] = 'light'
	mpl.rcParams['font.sans-serif'] = 'Helvetica Neue'
	mpl.rcParams['font.weight'] = 'light'
	# Font colors
	mpl.rcParams['text.color'] = greys[0]
	mpl.rcParams['axes.labelcolor'] = greys[0]
	mpl.rcParams['xtick.color'] = greys[0]
	mpl.rcParams['ytick.color'] = greys[0]
	# Font sizes
	mpl.rcParams['figure.titlesize'] = 16
	mpl.rcParams['axes.titlesize'] = 16
	mpl.rcParams['axes.labelsize'] = 16
	mpl.rcParams['xtick.labelsize'] = 16
	mpl.rcParams['ytick.labelsize'] = 16
	# Border colors
	mpl.rcParams['axes.edgecolor'] = greys[0]
	global sns
	# Changes seaborn violinplot edgecolor
	sns.categorical._Old_Violin = sns.categorical._ViolinPlotter
	class _My_ViolinPlotter(sns.categorical._Old_Violin):
		def __init__(self, *args, **kwargs):
			super(_My_ViolinPlotter, self).__init__(*args, **kwargs)
			self.gray=greys[0]
	sns.categorical._ViolinPlotter = _My_ViolinPlotter


def process_results(file, results):
	file_name = file.split('/')[-1].split('.')[0]
	dat = pd.read_csv(file, index_col = 0)
	dat = pd.melt(dat, id_vars=['destination region', 'source region', 'descendant tips', 'source node', '# descendants'], 
		    value_vars=['midroot time', 'lower midroot time', 'upper midroot time'])
	dat['file'] = file_name
	dat = dat.sort_values('value')
	dat = dat.assign(cumsum=dat.groupby(['destination region', 'variable']).cumcount()+1)
	results = results.append(dat)
	return(results.reset_index(drop=True))


# mostly from arviz
# but works directly with a list
# and also returns median
# https://github.com/arviz-devs/arviz/blob/master/arviz/stats/stats.py
def hpd(dat, qs):
    width = qs[2] - qs[0]
    # sorts from smallest to largest
    dat = sorted(dat)
    # length of data
    n = len(dat)
    # number of values we are keeping
    # thus, the index of the begining of 
    # the HPD interval must be <= this value
    # this gives us the tail of the distribution
    interval_idx_inc = int(np.floor(width * n))
    # number of values we are excluding
    # thus, possible number of HPD intervals
    # this gives us the head of the distribution
    n_intervals = n - interval_idx_inc
    # the width of each possible interval
    # for each possible head and tail value, 
    # what is the difference between them
    interval_width = [a_i - b_i for a_i, b_i in 
    				  zip(dat[interval_idx_inc:], 
    				  	  dat[:n_intervals])]
    # find the shortest interval
    min_idx = interval_width.index(min(interval_width))
    hpd_interval = (dat[min_idx], dat[min_idx+interval_idx_inc])
    dat_hpd = [item for item in dat if (item >= hpd_interval[0]) & (item <= hpd_interval[1])]
    dat_mid = np.quantile(dat_hpd, qs[1])
    return((hpd_interval[0], dat_mid, hpd_interval[1]))


def summarize_results(results):
	qs = (0.025, 0.5, 0.975)
	regions = list(set(results['destination region']))
	variables = list(set(results['variable']))
	n_import_data = pd.DataFrame()
	observations = list(set(results['file']))
	# this is hacky, i'm bad at pandas groupby
	# we don't care about introductions leading to general wisconsin tips
	for region in ['MHD', 'UWHC']: 	
		these_results = results[(results['destination region'] == region)]
		time_range = [min(these_results['value']), max(these_results['value'])]
		times = np.arange(np.floor(time_range[0]*100)/100, np.ceil(time_range[1]*100)/100, 1/1000)
		for variable in variables:
			these_results = results[(results['destination region'] == region) & 
									(results['variable'] == variable)]
			chunked_data = pd.DataFrame(index=times)
			chunked_data['region'] = region
			chunked_data['variable'] = variable
			chunked_data['importations'] = [[len(these_results[(these_results['value'] < time) & 
								(these_results['file'] == observation)]) for observation in observations] for time in times]
			chunked_data['importations_summary'] = chunked_data['importations'].apply(hpd, qs=(0.025, 0.5, 0.975))
			chunked_data[['low', 'mid', 'high']] = \
				pd.DataFrame(chunked_data['importations_summary'].tolist(), index=chunked_data.index)
			n_import_data = n_import_data.append(chunked_data)
	return(n_import_data)


def plot_tre(ax, tre, tre_states_dict, inner_colors, outer_colors, names):
	def branch_c_func(k):
		if k.branchType == 'node':
			if 'label' in k.traits.keys():
				name = k.traits['label']
			else:
				return("#575c66")
		else:
			name = k.numName
		if tre_states_dict[name] in inner_colors.keys():
			return(inner_colors[tre_states_dict[name]])
		else:
			return("#575c66")
	def tip_c_func(k):
		if tre_states_dict[k.numName] in inner_colors.keys():
			return(inner_colors[tre_states_dict[k.numName]])
		else:
			return(None)
	def tip_c_o_func(k):
		if tre_states_dict[k.numName] in inner_colors.keys():
			return(outer_colors[tre_states_dict[k.numName]])
		else:
			return(None)
	def tip_s_func(k):
		if tre_states_dict[k.numName] in inner_colors.keys():
			return(50-30*k.height/tre.treeHeight)
		else:
			return(0.0)
	def bw_func(k):
		if k.branchType == 'node':
			n_children = len(k.leaves)
			if n_children > 50: 
				return(1.75)
			else:
				return(0.2)
		else:
			return(0.1)
	tre.plotTree(ax, x_attr=lambda k: k.absoluteTime, 
				 colour_function=branch_c_func, 
				 branchWidth=bw_func,
				 midbranch=True)
	tre.plotPoints(ax, x_attr=lambda k: k.absoluteTime,
	                    colour_function= tip_c_o_func,
	                    size_function=tip_s_func)
	tre.plotPoints(ax, x_attr=lambda k: k.absoluteTime,
	                    colour_function=tip_c_func,
	                    size_function=lambda k: 0.5*tip_s_func(k))
	ax.set_yticks([])
	y_pos = 0.1
	for name in sorted(list(set(names.keys())), reverse=False):
		ax.text(-0.02, y_pos, name, color=outer_colors[names[name]],
		        size=28, transform=ax.transAxes)
		y_pos = y_pos - 0.05
	return(ax)


def plot_results(results, tre, tre_states_dict, outname):
	regions = list(set(results['destination region']))
	light_grey = '#CCCCCC'
	med_grey = '#888888'
	names = {'Milwaukee County': 'MHD', 'Dane County': 'UWHC', 'Wisconsin': 'Wisconsin'}
	inner_colors = {'MHD': '#457B9D', 'UWHC': '#C44536', 'Wisconsin': '#8B668B'}
	outer_colors = {'MHD': '#305882', 'UWHC': '#933831', 'Wisconsin': '#754C78'}
	fig, axs = plt.subplots(3, 1, figsize=(6.4*1.5, 4.8*3), 
							constrained_layout=True, 
							gridspec_kw = {'height_ratios':[5, 1, 1]})
	axs[0] = plot_tre(axs[0], tre, tre_states_dict, inner_colors, outer_colors, names)
	for file in list(set(results['file'])):
		for i, region in enumerate([item for item in regions if item != 'Wisconsin']):
			for variable in list(set(results['variable'])):
				dat = results[(results['file'] == file) & 
							  (results['destination region'] == region) & 
							  (results['variable'] == variable)]
				_ = axs[1].step(dat['value'],
						 dat['cumsum'], alpha=0.2, color=outer_colors[region])
				dat = dat['value'].repeat(dat['# descendants'])
				_ = sns.distplot(dat, kde=True, hist=False, ax=axs[2], 
							 color=outer_colors[region], kde_kws={"alpha": 0.25})
	axs[1].set_ylabel('# of introductions')
	axs[1].set_ylim(0, max(results['cumsum'])+1)
	axs[1].set_yticks([0, 10, 20, 30])
	axs[2].set_ylabel('Descendant sample\ndensity')
	axs[2].set_yticklabels(['' for item in axs[2].get_yticks()])
	axs[2].tick_params(axis='y', which='both', length=0)
	for i,ax in enumerate(axs):
		subplot_label = ax.text(-0.1, 1, string.ascii_uppercase[i], color=greys[0],
                size=24, transform=ax.transAxes, fontweight='bold')
		ax.set_xlim((2019.93, 2020.45))
		ax.set_xticks(pd.to_datetime(['01-01-20', '02-01-20', '03-01-20', '04-01-20', '05-01-20']).\
					  to_series().apply(numeric_date).to_list())
		ax.grid(axis='x', ls='-', color=light_grey, alpha=0.5)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_xlabel(None)
		if i == 0:
			ax.spines['left'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
		if ax == axs[-1]:
			ax.set_xticklabels(['Jan.\n2020', 'Feb.\n2020', 'March\n2020', 'April\n2020', 'May\n2020'])
		else:
			ax.set_xticklabels([])
			ax.tick_params(axis='x', which='both',length=0)
	fig.savefig(f'{outname}.pdf')
	plt.close()


def print_results(results, outdir):
	qs = (0.025, 0.5, 0.975)
	formatted_results = []
	results.groupby('file')['source node'].count()
	region_introductions = {}
	# all WI includes UWHC, MHD, and unspecified WI
	region_introductions['all_wi'] = hpd(results.groupby('file')['source node'].nunique(), qs)
	# to do take as argument
	regions = ('UWHC', 'MHD')
	for region in regions: 
		region_dat = results[results['destination region'] == region]
		region_introductions[region] = hpd(region_dat.groupby('file')['source node'].nunique(), qs)
	# nodes which give rise to tips from BOTH MHD and UWHC
	# first, nodes which give rise to tips from either MHD or UWHC
	either_dat = results[results['destination region'].isin(regions)]
	# leave just midroot time as we're counting nodes and don't care about uncertainty in the time
	either_dat = either_dat[either_dat['variable'] == 'midroot time']
	# groupby file, subset to just source node, and count how many times each node appears
	both_dat = \
		either_dat.groupby('file')['source node'].value_counts()[\
		either_dat.groupby('file')['source node'].value_counts() >= len(regions)]
	region_introductions['both'] = hpd([len(set(index[1])) for index in both_dat.index], qs)
	# probably an easier way to do this with groupby
	both_dat_data = either_dat[either_dat.apply(
		lambda k: k['source node'] in both_dat_nodes[k['file']].index, axis=1)]
	# counts how many descendant tips into each region from shared introduction events
	shared_intro_n = {}
	shared_intro_prop = {}
	for region in regions:
		print(region)
		# how many tips are from this region total?
		region_results = results[results['destination region'] == region]
		region_results['n'] = region_results['descendant tips'].apply(lambda k: len(ast.literal_eval(k)))
		region_n = region_results.groupby('file')['n'].sum()
		region_data = both_dat_data[both_dat_data['destination region'] == region]
		region_data['n'] = region_data['descendant tips'].apply(lambda k: len(ast.literal_eval(k)))
		# number of samples from this region from shared introductions into Wisconsin
		shared_intro_n[region] = hpd(region_data.groupby('file')['n'].sum(), qs)
		# proportion of samples from this region from shared introductions into Wisconsin
		shared_intro_prop[region] = hpd(region_data.groupby('file')['n'].sum()/region_n, qs)
	# todo make loop
	n_all = f'{int(region_introductions["all_wi"][1])} [{region_introductions["all_wi"][0]}, {region_introductions["all_wi"][2]}]'
	n_dane = f'{int(region_introductions["UWHC"][1])} [{region_introductions["UWHC"][0]}, {region_introductions["UWHC"][2]}]' 
	n_mke = f'{int(region_introductions["MHD"][1])} [{region_introductions["MHD"][0]}, {region_introductions["MHD"][2]}]' 
	n_both = f'{int(region_introductions["both"][1])} [{region_introductions["both"][0]}, {region_introductions["both"][2]}]' 
	shared_dane_n = f'{int(shared_intro_n["UWHC"][1])} [{int(shared_intro_n["UWHC"][0])}, {int(shared_intro_n["UWHC"][2])}]'
	print(shared_dane_n)
	shared_dane_prop = f'{round(shared_intro_prop["UWHC"][1], 2)} [{round(shared_intro_prop["UWHC"][0], 2)}, {round(shared_intro_prop["UWHC"][2], 2)}]'
	shared_mke_n = f'{int(shared_intro_n["MHD"][1])} [{int(shared_intro_n["MHD"][0])}, {int(shared_intro_n["MHD"][2])}]'
	print(shared_mke_n)
	shared_mke_prop = f'{round(shared_intro_prop["MHD"][1], 2)} [{round(shared_intro_prop["MHD"][0], 2)}, {round(shared_intro_prop["MHD"][2], 2)}]'
	formatted_results.append([f'We estimate {n_all} introductions into Wisconsin.'])
	formatted_results.append([f'Of these, {n_dane} led to at least one tip from Dane County (UWHC)'])
	formatted_results.append([f'and {n_mke} led to at least one tip from Milwaukee (MHD)'])
	formatted_results.append([f'{n_both} led to at least one tip from Dane County and Milwaukee.'])
	formatted_results.append([f'{shared_dane_n} ({shared_dane_prop}) of the samples from Dane County were from introductions which also led to samples from Milwaukee'])
	formatted_results.append([f'{shared_mke_n} ({shared_mke_prop}) of the samples from Milwaukee were from introductions which also led to samples from Dane County'])
	pd.DataFrame(formatted_results).to_csv(f'{outdir}/results.txt', index=None, header=False)


def import_tre(tree_file, metadata_file):
	tre = bt.loadNewick(tree_file)
	#_ = tre.traverse_tree()
	tre_tip_names = [item.numName for item in tre.Objects if item.branchType == 'leaf']
	metadata_dict = parse_dates(metadata_file)
	metadata_dict = {key: item for key, item in metadata_dict.items() if key in tre_tip_names}
	max_time = max(metadata_dict.values())
	tre.setAbsoluteTime(max_time)
	tre.sortBranches()
	return(tre)
	

def plot_rarefaction_results(rarefaction_results, outname):
	inner_colors = {'MHD': '#457B9D', 'UWHC': '#C44536'}
	outer_colors = {'MHD': '#305882', 'UWHC': '#933831'}
	names = {'MHD': 'Milwaukee County', 'UWHC': 'Dane County'}
	fig, axs = plt.subplots(1, 1, figsize=(6.4, 4.8), 
							constrained_layout=True)
	sns.violinplot(x='n', y='importations', hue='region', 
				   cut=0, data=rarefaction_results, ax=axs, 
				   palette=inner_colors, saturation=1)
	regions = list(set(rarefaction_results['region']))
	y_pos = 0.9
	for region in regions:
		axs.text(0.02, y_pos, names[region], color=outer_colors[region],
				 size=20, transform=axs.transAxes)
		y_pos = y_pos - 0.07
	axs.get_legend().remove()
	axs.set_xlabel('Milwaukee and Dane County samples')
	axs.set_ylabel('# of introductions')
	fig.savefig(f'{outname}.pdf')
	plt.close(fig)


def main():
set_style()
parser = argparse.ArgumentParser()
# input files
parser.add_argument('--introduction_dir', 
					help='estimated introductions from each bootstrap replicate',
					default='results/subsampled_alignment_neighbors.fasta.ufboot_tres/*/*_importations.csv')
parser.add_argument('--outdir', default='results')
parser.add_argument('--outname', default='introductions',
					help='prefix of output files')
parser.add_argument('--metadata_file', default='data/metadata_adjusted.tsv',
					help='metadata file with dates')
parser.add_argument('--tree_file', default='results/ml_refined_time.newick',
					help='newick file to plot')
parser.add_argument('--outgroup',
				help='Outgroup to use in tree',
				default='Wuhan/Hu-1/2019')
parser.add_argument('--tree_states_file', default='results/ml_refined_node_states.csv',
					help='states of nodes/tips in tree')
parser.add_argument('--tree_regions', default='results/ml_refined_regions.csv',
					help='regions dictionary from treetime')
parser.add_argument('--rarefaction_dir', default='results/rarefaction/*/*_rarefaction.csv',
					help='directory with rarefaction results')
args = parser.parse_args()
results = pd.DataFrame()
	for file in glob.glob(args.introduction_dir):
	#for file in args.importation_files:
		results = process_results(file, results)
	results.to_csv(f'{args.outdir}/combined_results.csv', index=None)
	tre = import_tre(args.tree_file, args.metadata_file)
	regions = pd.read_csv(args.tree_regions)
	tre_states = pd.read_csv(args.tree_states_file)
	tre_states['state'] = \
		tre_states['state'].apply(lambda k: ast.literal_eval(k))
	tre_states['ML_state'] = \
		tre_states['state'].apply(lambda k: regions.loc[k.index(max(k)), 'region'])
	tre_states_dict = {item['name']: item['ML_state'] for index, item in tre_states.iterrows()}
	plot_results(results, tre, tre_states_dict, 'figures/tree_timeseries')
	print_results(results, args.outdir)
	rarefaction_results = pd.DataFrame()
	for file in glob.glob(args.rarefaction_dir):
		rarefaction_results = rarefaction_results.append(pd.read_csv(file))
	plot_rarefaction_results(rarefaction_results, f'figures/rarefaction')


if __name__ == "__main__":
	main()
