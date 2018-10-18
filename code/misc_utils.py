#!/usr/bin/env python3
import pandas as pd

def combine_mafs(mafs, output):
	mut_dfs = []
	for file in mafs:
		with open(file, 'r') as f:
			count = 0
			for line in f.readlines():
				if line.startswith('Hugo'):
					break
				count += 1
		mut_df = pd.read_csv(file, sep='\t', header=count)
		mut_dfs.append(mut_df)
	all_mafs = pd.concat(mut_dfs, ignore_index=True)
	all_mafs.to_csv(output, sep='\t', header=True, index=False)
	with open(output, 'r') as f:
		data = f.read()
	with open(output, 'w') as f:
		f.write('#version 2.4\n' + data)

def combine_cnvs(cnvs, samples, output):
	cnv_dfs = []
	for i, file in enumerate(cnvs):
		sample = samples[i]
		try:
			cnv_df = pd.read_csv(file, sep='\t', header=0, usecols=['Hugo_Symbol', 'class', 'depth', 'chromosome', 'start', 'end', 'weight', 'log2', 'ci_lo', 'ci_hi'])
			cnv_df['sample'] = sample
		except:
			cnv_df = pd.DataFrame(columns=['Hugo_Symbol', 'class', 'depth', 'chromosome', 'start', 'end', 'weight', 'log2', 'ci_lo', 'ci_hi', 'sample'])
		cnv_dfs.append(cnv_df)
	all_cnvs = pd.concat(cnv_dfs, ignore_index=True)
	all_cnvs.to_csv(output, sep='\t', header=True, index=False)

def create_mut_mats(muts, cnvs, mut_mat, cnv_mat):
	return 'hello'