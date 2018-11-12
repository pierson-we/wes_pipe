#!/usr/bin/env python3
import pandas as pd
import numpy as np

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

def filter_pindel(pindel_files, sample_dict, project_dir, all_samples_output, min_reads, min_qual, max_inv_length):
	pindel_dfs = []
	pindel_dict = {sample: [] for sample in sample_dict}
	for file in pindel_files:
		mut_type = file.split('.filtered.tsv')[0].split('_')[-1]
		with open(file, 'r') as f:
			lines = f.readlines()
		i = 0
		num_lines = len(lines)
		while i < num_lines:
			while i+1 < num_lines and not lines[i].startswith('#'):
				i += 1
			if i+1 >= num_lines:
				break
			variant_dict = {}
			summary_line = lines[i+1].split()
			if mut_type == 'INV' and int(summary_line[2]) > max_inv_length:
				i += 1
				while i+1 < num_lines and not lines[i].startswith('#'):
					i += 1
				continue
			i += 2
			while i+1 < num_lines and not lines[i+1].startswith('#'):
				i += 1
				read_line = lines[i].split('\t')
				if len(read_line) == 1:
					if read_line[0].startswith('-'):
						while i+1 < num_lines and not lines[i].startswith('#'):
							i += 1
						break
					i += 1
					if i > num_lines:
						break
					read_line += lines[i].split('\t')
				while '' in read_line:
					del read_line[read_line.index('')]
				read_qual = int(read_line[3])
				if read_qual >= min_qual:
					sample = read_line[4]
					if not sample in variant_dict:
						variant_dict[sample] = []
					variant_dict[sample].append(read_qual)
			for sample in variant_dict:
				if len(variant_dict[sample]) >= min_reads:
					info = ['type=%s' % mut_type, 'length=%s' % summary_line[2], 'reads=%s' % len(variant_dict[sample]), 'avg_qual=%s' % np.mean(variant_dict[sample])]
					pindel_dict[sample].append({
						'chr': summary_line[7],
						'start': summary_line[9],
						'end': summary_line[10],
						'gffTag': ';'.join(info)
						})
	for sample in sample_dict:
		sample_df = pd.DataFrame(data=pindel_dict[sample], columns=['chr', 'start', 'end', 'gffTag'])
		with open(sample_dict[sample], 'w') as f:
			f.write('#gffTags\n')
			f.write(sample_df.to_csv(sep='\t', header=False, index=False))
		sample_df['sample'] = sample
		pindel_dfs.append(sample_df)
	all_pindels = pd.concat(pindel_dfs, ignore_index=True)
	all_pindels.to_csv(all_samples_output, sep='\t', header=True, index=False)



	# 	def parse(row, pindel_dict):
	# 		sample_count = int(row[28])
	# 		samples = [row[32+7*x] for x in range(0, sample_count)]
	# 		sample_reads = [row[34+7*x] for x in range(0, sample_count)]
	# 		for i, reads in enumerate(sample_reads):
	# 			if reads >= min_reads:
	# 				info = ['type=%s' % row[2], 'length=%s' % row[3], 'reads=%s' % reads]
	# 				pindel_dict[samples[i]].append({
	# 				'chr': row[8],
	# 				'start': row[10],
	# 				'end': row[11],
	# 				'gffTag': ';'.join(info)
	# 				})
		
	# 	pindel_df_raw = pd.read_csv(file, sep='\s+', header=None)
	# 	pindel_df_raw.columns = [i for i in range(1, pindel_df_raw.shape[1]+1)]
	# 	# print(pindel_df_raw.head().to_csv())
	# 	pindel_df_raw.apply(parse, axis=1, pindel_dict=pindel_dict)
	# for sample in sample_dict:
	# 	sample_df = pd.DataFrame(data=pindel_dict[sample], columns=['chr', 'start', 'end', 'gffTag'])
	# 	with open(sample_dict[sample], 'w') as f:
	# 		f.write('#gffTags\n')
	# 		f.write(sample_df.to_csv(sep='\t', header=False, index=False))
	# 	sample_df['sample'] = sample
	# 	pindel_dfs.append(sample_df)
	# all_pindels = pd.concat(pindel_dfs, ignore_index=True)
	# all_pindels.to_csv(all_samples_output, sep='\t', header=True, index=False)

def format_pindel(pindel_files, sample_dict, project_dir, all_samples_output, min_reads):
	pindel_dfs = []
	pindel_dict = {sample: [] for sample in sample_dict}
	for file in pindel_files:
		mut_type = file.split('.filtered.tsv')[0].split('_')[-1]
		def parse(row, pindel_dict):
			sample_count = int(row[28])
			samples = [row[32+7*x] for x in range(0, sample_count)]
			sample_reads = [row[34+7*x] for x in range(0, sample_count)]
			for i, reads in enumerate(sample_reads):
				if reads >= min_reads:
					info = ['type=%s' % row[2], 'length=%s' % row[3], 'reads=%s' % reads]
					pindel_dict[samples[i]].append({
					'chr': row[8],
					'start': row[10],
					'end': row[11],
					'gffTag': ';'.join(info)
					})
		
		pindel_df_raw = pd.read_csv(file, sep='\s+', header=None)
		pindel_df_raw.columns = [i for i in range(1, pindel_df_raw.shape[1]+1)]
		# print(pindel_df_raw.head().to_csv())
		pindel_df_raw.apply(parse, axis=1, pindel_dict=pindel_dict)
	for sample in sample_dict:
		sample_df = pd.DataFrame(data=pindel_dict[sample], columns=['chr', 'start', 'end', 'gffTag'])
		with open(sample_dict[sample], 'w') as f:
			f.write('#gffTags\n')
			f.write(sample_df.to_csv(sep='\t', header=False, index=False))
		sample_df['sample'] = sample
		pindel_dfs.append(sample_df)
	all_pindels = pd.concat(pindel_dfs, ignore_index=True)
	all_pindels.to_csv(all_samples_output, sep='\t', header=True, index=False)


def create_mut_mats(muts, cnvs, mut_mat, cnv_mat):
	return 'hello'