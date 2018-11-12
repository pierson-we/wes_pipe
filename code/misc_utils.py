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

	def create_pon(row, pon):
		variant_type = row.gffTag.split(';')[0].split('=')[-1]
		variant_length = row.gffTag.split(';')[1].split('=')[-1]
		variant_id = '_'.join([row.start, row.end, variant_type, variant_length])
		if row.chr not in pon:
			pon[row.chr] = []
		if variant_id not in pon[row.chr]:
			pon[row.chr].append(variant_id)

	def filter_pon(row, pon):
		variant_type = row.gffTag.split(';')[0].split('=')[-1]
		variant_length = row.gffTag.split(';')[1].split('=')[-1]
		variant_id = '_'.join([row.start, row.end, variant_type, variant_length])
		if row.chr in pon:
			if variant_id in pon[row.chr]:
				return False
			else:
				return True
		else:
			return True

	pon = {}

	for sample in sample_dict:
		if '_N' in sample:
			sample_df = pd.DataFrame(data=pindel_dict[sample], columns=['chr', 'start', 'end', 'gffTag'])
			with open(sample_dict[sample], 'w') as f:
				f.write('#gffTags\n')
				f.write(sample_df.to_csv(sep='\t', header=False, index=False))
			sample_df['sample'] = sample
			pindel_dfs.append(sample_df)
			sample_df.apply(create_pon, axis=1, pon=pon)

	for sample in sample_dict:
		if '_T' in sample:
			sample_df = pd.DataFrame(data=pindel_dict[sample], columns=['chr', 'start', 'end', 'gffTag'])
			sample_df = sample_df[sample_df.apply(filter_pon, axis=1, pon=pon)]
			with open(sample_dict[sample], 'w') as f:
				f.write('#gffTags\n')
				f.write(sample_df.to_csv(sep='\t', header=False, index=False))
			sample_df['sample'] = sample
			pindel_dfs.append(sample_df)
			sample_df.apply(create_pon, axis=1, pon=pon)

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


def create_mut_mats(mafs, cnvs, pindel, mut_mat_file, cnv_mat_file, mut_counts_file):
	mut_samples = []
	mut_counts = []
	mut_dfs = []
	for file in mafs:
		sample = file.split('/')[-1].split('.')[0]
		with open(file, 'r') as f:
			count = 0
			for line in f.readlines():
				if line.startswith('Hugo'):
					break
				count += 1
		mut_df = pd.read_csv(file, sep='\t', header=count, usecols=['Hugo_Symbol', 'Variant_Classification', 'FILTER', 'dbSNP_RS'])
		mut_df = mut_df[mut_df['FILTER'] == 'PASS']
		mut_counts.append(mut_df.shape[0])
		mut_types = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins', 'Frame_Shift_Ins']
		def filter_mut_types(row, mut_types):
			if row.Variant_Classification in mut_types:
				if row.Variant_Classification == 'Missense_Mutation':
					if row.dbSNP_RS != 'novel':
						return True
					else:
						return False
				else:
					return True
				# return True
			else:
				return False

		mut_df = mut_df[mut_df.apply(filter_mut_types, axis=1, mut_types=mut_types)]
		mut_dfs.append(mut_df)
		mut_samples.append(sample)
	with open(mut_counts_file, 'w') as f:
		f.write('\t'.join(mut_samples))
		f.write('\n')
		f.write('\t'.join([str(x) for x in mut_counts]))

	mut_pindel_dfs = []
	mut_pindel_samples = []
	for file in pindel:
		sample = file.split('/')[-1].split('.')[0]
		pindel_df = pd.read_csv(file, skiprows=1, header=None, sep='\t', names=['chr', 'start', 'end', 'info', 'genes'])
		def parse_pindel(row, new_rows):
			mut_type = row.info.split(';')[0].split('=')[1]
			length = mut_type = row.info.split(';')[1].split('=')[1]
			genes = row.genes.split(';')
			if int(length) > 3 and int(length) % 3 != 0:
				for gene in genes:
					new_rows.append({'Hugo_Symbol': gene, 'Variant_Classification': mut_type, 'FILTER': 'PASS', 'dbSNP_RS': ''})
		pindel_data = []
		pindel_df.apply(parse_pindel, new_rows=pindel_data)
		parsed_pindel_df = pd.DataFrame(pindel_data, columns=['Hugo_Symbol', 'Variant_Classification', 'FILTER', 'dbSNP_RS'])
		mut_df = mut_dfs[mut_samples.index(sample)]
		mut_pindel_df = pd.concat([mut_df, parsed_pindel_df], ignore_index=True)
		mut_pindel_dfs.append(mut_pindel_df)
		mut_pindel_samples.append(sample)

	cnv_samples = []
	cnv_dfs = []

	def filter_ci(row):
		if row['class'] == 'amp':
			if row.log2 > row.ci_lo:
				return True
			else:
				return False
		elif row['class'] == 'del':
			if row.log2 < row.ci_hi:
				return True
			else:
				return False
		else:
			return False
	for file in cnvs:
		sample = file.split('/')[-1].split('.')[0]
		try:
			cnv_df = pd.read_csv(file, sep='\t', header=0, usecols=['Hugo_Symbol', 'class', 'depth', 'chromosome', 'start', 'end', 'weight', 'log2']) #, 'ci_lo', 'ci_hi'])
			cnv_df = cnv_df[cnv_df['class'] != 'wt']
			# cnv_df = cnv_df[cnv_df['depth'] >= 100]
			cnv_df = cnv_df[cnv_df['chromosome'] != 'Y']
			# cnv_df = cnv_df[cnv_df.apply(filter_ci, axis=1)]
			cnv_df['sample'] = sample
			cnv_df['depth_per_kb'] = cnv_df['depth'] / (cnv_df['end'] - cnv_df['start'])*1000
			# min_depth = np.percentile(cnv_df.depth_per_kb, 25)
			# cnv_df = cnv_df[cnv_df['depth_per_kb'] > min_depth]
			# cnv_df = cnv_df[['Hugo_Symbol', 'class']]
		except:
			cnv_df = pd.DataFrame(columns=['Hugo_Symbol', 'class', 'depth', 'chromosome', 'start', 'end', 'weight', 'log2', 'sample', 'depth_per_kb']) # , 'ci_lo', 'ci_hi'])
		# print(sample)
		# print(cnv_df.shape)
		cnv_dfs.append(cnv_df)
		cnv_samples.append(sample)
	all_cnvs = pd.concat(cnv_dfs, ignore_index=True)
	all_cnvs.to_csv(os.path.join(out_dir, 'all_cnvs.tsv'), sep='\t', header=True, index=False)
	min_depth = np.percentile(all_cnvs.depth_per_kb, 25)
	# print('min depth/kb: %s' % min_depth)
	cnv_dfs = [cnv_df[cnv_df['depth_per_kb'] > min_depth] for cnv_df in cnv_dfs]
	assert mut_samples == cnv_samples

	all_genes = []
	for mut_df in mut_dfs:
		genes = mut_df.Hugo_Symbol.unique().tolist()
		all_genes = list(set(all_genes + genes))
	for cnv_df in cnv_dfs:
		genes = cnv_df.Hugo_Symbol.unique().tolist()
		all_genes = list(set(all_genes + genes))

	all_genes = sorted(all_genes)
	
	mut_mat = pd.DataFrame(columns=mut_samples, index=all_genes)
	cnv_mat = pd.DataFrame(columns=cnv_samples, index=all_genes)

	for i, mut_df in enumerate(mut_pindel_dfs):
		for gene in mut_df.Hugo_Symbol.unique().tolist():
			gene_df = mut_df[mut_df['Hugo_Symbol'] == gene]
			variant_types = gene_df.Variant_Classification.unique().tolist()
			if 'D' variant_types or 'INS' in variant_types or 'INV' in variant_types or 'TD' in variant_types:
				mut_mat.loc[gene, mut_pindel_samples[i]] = 5 # SV = 5
			elif 'Nonsense_Mutation' in variant_types or 'Frame_Shift_Ins' in variant_types or 'Frame_Shift_Del' in variant_types:
				mut_mat.loc[gene, mut_pindel_samples[i]] = 1 # nonsense/frameshift = 1
			else:
				mut_mat.loc[gene, mut_samples[i]] = 2 # missense = 2
	mut_mat.fillna(value=0, inplace=True)

	for i, cnv_df in enumerate(cnv_dfs):
		for gene in cnv_df.Hugo_Symbol.unique().tolist():
			variant_types = cnv_df[cnv_df['Hugo_Symbol'] == gene]['class'].unique().tolist()
			# variant_types = gene_df['class'].unique().tolist()
			# if 'Nonsense_Mutation' in variant_types or 'Frame_Shift_Ins' in variant_types or 'Frame_Shift_Del' in variant_types:
			if len(variant_types) == 1:
				if variant_types[0] == 'amp':
					cnv_mat.loc[gene, cnv_samples[i]] = 3 # amplication = 3
				elif variant_types[0] == 'del':
					cnv_mat.loc[gene, cnv_samples[i]] = 4 # deletion = 4
				else:
					print('found a wt straggler')
			else:
				# print(cnv_df[cnv_df['Hugo_Symbol'] == gene])
				print('%s, %s: uh ohhhhhhh' % (cnv_samples[i], gene))
			# else:
			# 	cnv_mat.loc[gene, cnv_samples[i]] = 'del'
	cnv_mat.fillna(value=0, inplace=True)
	# print(cnv_mat.head())
	# print(cnv_mat.shape)

	mut_mat.to_csv(mut_mat_file, header=True, index=True, sep='\t')
	cnv_mat.to_csv(cnv_mat_file, header=True, index=True, sep='\t')