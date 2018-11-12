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
			cnv_df = pd.read_csv(file, sep='\t', header=0, usecols=['Hugo_Symbol', 'class', 'depth', 'chromosome', 'start', 'end', 'weight', 'log2']) #, 'ci_lo', 'ci_hi'])
			cnv_df['sample'] = sample
		except:
			cnv_df = pd.DataFrame(columns=['Hugo_Symbol', 'class', 'depth', 'chromosome', 'start', 'end', 'weight', 'log2', 'sample']) #, 'ci_lo', 'ci_hi'])
		cnv_dfs.append(cnv_df)
	all_cnvs = pd.concat(cnv_dfs, ignore_index=True)
	all_cnvs.to_csv(output, sep='\t', header=True, index=False)

def format_pindel(pindel_files, sample_dict, project_dir, all_samples_output):
	pindel_dfs = []
	pindel_dict = {sample: [] for sample in sample_dict}
	for file in pindel_files:
		mut_type = file.split('.filtered.tsv')[0].split('_')[-1]
		def parse(row, pindel_dict):
			sample_count = int(row[28])
			samples = [row[32+x].split()[0] for x in range(0, sample_count)]
			sample_reads = [row[32+x].split()[2] for x in range(0, sample_count)]
			for i, reads in enumerate(sample_reads):
				if reads > 0:
					pindel_dict[samples[i]].append({
					'chr': 'chr' + str(row[8]),
					'start': row[10],
					'end': row[11],
					'gffTag': row[2] + '_' + str(reads)
					})
		
		pindel_df_raw = pd.read_csv(file, sep=None, header=None)
		pindel_df_raw.columns = [i for i in range(1, pindel_df_raw.shape[1]+1)]
		print(pindel_df_raw.head().to_csv())
		pindel_df_raw.apply(parse, axis=1, pindel_dict=pindel_dict)
	for sample in sample_dict:
		sample_df = pd.DataFrame(data=pindel_dict[sample], columns=['chr', 'start', 'end', 'gffTag'])
		with open(sample_dict[sample], 'w') as f:
			f.write('#gffTags\n')
			f.write(sample_df.to_csv(sep='\t', header=False, index=False))
		sample_df['sample'] = sample
		pindel_dfs.append(sample_df)
	all_pindels = pd.concat(cnv_dfs, ignore_index=True)
	all_pindels.to_csv(all_samples_output, sep='\t', header=True, index=False)


def create_mut_mats(muts, cnvs, mut_mat, cnv_mat):
	return 'hello'