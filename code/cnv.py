#!/usr/bin/env python3
import subprocess
import sys
import pandas as pd
import os
import time
import random
import luigi
import multiprocessing
import bam_processing
import variant_calling
import pipeline_utils
import global_vars


class cnvkit_prep(luigi.Task):
	max_threads = luigi.IntParameter()

	project_dir = luigi.Parameter()
	case_dict = luigi.DictParameter()
	cfg = luigi.DictParameter()

	def requires(self):
		return [bam_processing.index_bam(sample=case_name + '_N', fastq_file=self.case_dict[case_name]['N'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg) for case_name in self.case_dict if self.case_dict[case_name]['N'] != ''] + \
		[bam_processing.index_bam(sample=case_name + '_T', fastq_file=self.case_dict[case_name]['T'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg) for case_name in self.case_dict]

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'ref', 'targets.bed')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'ref', 'antitargets.bed')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'ref', 'access.bed'))]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		cmd = 'python3 %s target %s --annotate %s -o %s' % (self.cfg['cnvkit_location'], self.cfg['library_bed'], self.cfg['refFlat'], os.path.join(self.project_dir, 'output', 'cnvkit', 'ref', 'targets.bed'))
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

		cmd = 'python3 %s access %s -o %s' % (self.cfg['cnvkit_location'], self.cfg['fasta_file'], os.path.join(self.project_dir, 'output', 'cnvkit', 'ref', 'access.bed'))
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

		cmd = 'python3 %s antitarget %s -g %s -o %s' % (self.cfg['cnvkit_location'], os.path.join(self.project_dir, 'output', 'cnvkit', 'ref', 'targets.bed'), os.path.join(self.project_dir, 'output', 'cnvkit', 'ref', 'access.bed'), os.path.join(self.project_dir, 'output', 'cnvkit', 'ref', 'antitargets.bed'))
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

		# cmd = 'python3 %s autobin %s -t %s -g %s --annotate %s' % (cnvkit, bam_base, targets_bed, access_bed, refFlat)
		# cmd = cmd.split(' ')
		# p = subprocess.Popen(cmd)
		# p.wait()

class coverage(luigi.Task):
	case = luigi.Parameter()
	max_threads = luigi.IntParameter()

	project_dir = luigi.Parameter()
	case_dict = luigi.DictParameter()
	cfg = luigi.DictParameter()

	def requires(self):
		requirements = [cnvkit_prep(max_threads=self.max_threads, project_dir=self.project_dir, cfg=self.cfg, case_dict=self.case_dict)]
		if self.case_dict[self.case]['N']:
			return requirements + [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.case_dict[self.case]['T'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), bam_processing.index_bam(sample=self.case + '_N', fastq_file=self.case_dict[self.case]['N'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), variant_calling.mutect_pon(case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]
		else:
			return requirements + [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.case_dict[self.case]['T'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), variant_calling.mutect_pon(case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]
	
	def outputs(self):
		tumor_out = [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'coverage', '%s_T.targetcoverage.cnn' % self.case)),
		luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'coverage', '%s_T.antitargetcoverage.cnn' % self.case))]
		if self.case_dict[self.case]['N'] != '':
			return tumor_out + [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'coverage', '%s_N.targetcoverage.cnn' % self.case)),
			luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'coverage', '%s_N.antitargetcoverage.cnn' % self.case))]
	
	def run(self):
		for output in self.output()():
			pipeline_utils.confirm_path(output.path)
		
		cmd = 'python3 %s coverage %s %s -o %s' % (self.cfg['cnvkit_location'], self.input()[1][0].path, self.input()[0].path, self.output()[0].path)
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

		cmd = 'python3 %s coverage %s %s -o %s' % (self.cfg['cnvkit_location'], self.input()[1][0].path, self.input()[1].path, self.output()[1].path)
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

		if self.case_dict[self.case]['N'] != '':
			cmd = 'python3 %s coverage %s %s -o %s' % (self.cfg['cnvkit_location'], self.input()[2][0].path, self.input()[0].path, self.output()[2].path)
			cmd = cmd.split(' ')
			pipeline_utils.command_call(cmd, self.output())

			cmd = 'python3 %s coverage %s %s -o %s' % (self.cfg['cnvkit_location'], self.input()[2][0].path, self.input()[1].path, self.output()[3].path)
			cmd = cmd.split(' ')
			pipeline_utils.command_call(cmd, self.output())

class reference(luigi.Task):
	max_threads = luigi.IntParameter()

	project_dir = luigi.Parameter()
	case_dict = luigi.DictParameter()
	cfg = luigi.DictParameter()

	def requires(self):
		return [coverage(case=case, max_threads=self.max_threads, project_dir=self.project_dir, cfg=self.cfg, case_dict=self.case_dict) for case in self.case_dict if self.case_dict[case]['N'] != '']
	
	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'reference.cnn'))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		cmd = 'python3 %s reference %s --fasta %s -o %s' % (self.cfg['cnvkit_location'], os.path.join(self.project_dir, 'output', 'cnvkit', 'coverage', '*N.{,anti}targetcoverage.cnn'), self.cfg['fasta_file'], self.output().path)
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

class fix(luigi.Task):
	case = luigi.Parameter()
	max_threads = luigi.IntParameter()

	project_dir = luigi.Parameter()
	case_dict = luigi.DictParameter()
	cfg = luigi.DictParameter()

	def requires(self):
		return [coverage(case=self.case, max_threads=self.max_threads, project_dir=self.project_dir, cfg=self.cfg, case_dict=self.case_dict), reference(max_threads=self.max_threads, project_dir=self.project_dir, cfg=self.cfg, case_dict=self.case_dict)]

	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'coverage', '%s_T.cnr' % self.case))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		cmd = 'python3 %s fix %s %s %s -o %s' % (self.cfg['cnvkit_location'], os.path.join(self.project_dir, 'output', 'cnvkit', 'coverage', '%s_T.targetcoverage.cnn' % self.case), os.path.join(self.project_dir, 'output', 'cnvkit', 'coverage', '%s_T.antitargetcoverage.cnn' % self.case), self.input()[1].path, self.output().path)
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

class refine_cnv(luigi.Task):
	case = luigi.Parameter()
	max_threads = luigi.IntParameter()

	project_dir = luigi.Parameter()
	case_dict = luigi.DictParameter()
	cfg = luigi.DictParameter()

	def requires(self):
		return fix(case=self.case, max_threads=self.max_threads, project_dir=self.project_dir, cfg=self.cfg, case_dict=self.case_dict)

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'segment', '%s_T.segment.cns' % self.case)),
		luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'segment', '%s_T.call.cns' % self.case)),
		luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'segment', '%s_T.scatter.png' % self.case)),
		luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'segment', '%s_T.segmetrics.cns' % self.case)),
		luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'segment', '%s_T.genemetrics.tsv' % self.case)),
		luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'segment', '%s_T.trusted_genes.tsv' % self.case))]
	
	def run(self):
		for output in self.output()():
			pipeline_utils.confirm_path(output.path)

		cmd = 'python3 %s segment %s -m %s --drop-low-coverage -o %s' % (self.cfg['cnvkit_location'], self.input().path, self.cfg['cnv_seg_method'], self.output()[0].path)
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

		cmd = 'python3 %s call %s -i %s_T -m threshold --ci --ampdel --sem -o %s_T.call.cns' % (self.cfg['cnvkit_location'], self.output()[0].path, self.case, self.output()[1].path)
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

		cmd = 'python3 %s scatter %s -s %s -o %s' % (self.cfg['cnvkit_location'], self.input().path, self.output()[1].path, self.output()[2].path)
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

		cmd = 'python3 %s segmetrics %s -s %s --ci --std --mean -o %s' % (self.cfg['cnvkit_location'], self.input().path, self.output()[1].path, self.output()[3].path)
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

		cmd = 'python3 %s genemetrics %s -s %s -t %s -m %s -o %s' % (self.cfg['cnvkit_location'], self.input().path, self.output()[3].path, self.cfg['cnvkit_genemetrics_threshold'], self.cfg['cnvkit_genemetrics_minprobes'], self.output()[4].path)
		cmd = cmd.split(' ')
		cmd = cmd.split(' ')
		pipeline_utils.command_call(cmd, self.output())

		#repurpose command call for pipe command calls
		wait_time = random.uniform(0,3)
		time.sleep(wait_time)
		sys.stdout.flush()
		while not pipeline_utils.add_thread_count(global_vars.thread_file, 1):
			time.sleep(sleep_time)
		
		cmd = 'python3 %s genemetrics %s_T.cnr -s %s' % (self.cfg['cnvkit_location'], self.input().path, self.output()[3].path)
		cmd = cmd.split(' ')
		p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
		# outs, err = p.communicate()
		cmd = 'tail -n+2'
		cmd = cmd.split(' ')
		p2 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=p1.stdout)
		# outs, err = p.communicate()
		cmd = 'cut -f1'
		cmd = cmd.split(' ')
		p3 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=p2.stdout)
		# outs, err = p.communicate()
		cmd = 'sort'
		cmd = cmd.split(' ')
		p4 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=p3.stdout)
		outs, err = p4.communicate()
		with open('%s_segment_genes.txt' % self.case, 'wb') as f:
			f.write(outs)

		cmd = 'python3 %s genemetrics %s' % (self.cfg['cnvkit_location'], self.input().path)
		cmd = cmd.split(' ')
		p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
		# outs, err = p.communicate()
		cmd = 'tail -n+2'
		cmd = cmd.split(' ')
		p2 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=p1.stdout)
		# outs, err = p.communicate()
		cmd = 'cut -f1'
		cmd = cmd.split(' ')
		p3 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=p2.stdout)
		# outs, err = p.communicate()
		cmd = 'sort'
		cmd = cmd.split(' ')
		p4 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=p3.stdout)
		outs, err = p4.communicate()
		with open('%s_ratio_genes.txt' % self.case, 'wb') as f:
			f.write(outs)

		cmd = 'comm -12 %s_ratio_genes.txt %s_segment_genes.txt' % (self.case, self.case)
		cmd = cmd.split(' ')
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
		outs, err = p.communicate()
		with open('%s_trusted_genes.txt' % self.case, 'wb') as f:
			f.write(outs)
		with open('%s_trusted_genes.txt' % self.case, 'r') as f:
			trusted_genes = f.read().split('\n')
		
		while not pipeline_utils.sub_thread_count(global_vars.thread_file, 1):
			time.sleep(sleep_time)

		seg = pd.read_csv(self.output()[4].path, sep='\t', header=0)

		def filter_genes(row, genes):
			if row.gene in genes:
				return True
			else:
				return False
		# # seg = seg[seg.apply(filter_ci, axis=1)]
		# # seg = seg[seg.apply(filter_pten, axis=1)]
		seg = seg[seg.apply(filter_genes, axis=1, genes=trusted_genes)]
		def assign_class(row):
			if row.cn > 2:
				return 'amp'
			elif row.cn < 2:
				return 'del'
			else:
				return 'wt'
		seg['class'] = seg.apply(assign_class, axis=1)
		seg['Tumor_Sample_Barcode'] = '%s_T' % self.case
		seg['FILTER'] = 'PASS'
		seg['Variant_Classification'] = 'SV'
		seg.rename({'gene': 'Hugo_Symbol'}, axis='columns', inplace=True)
		print('postfilter: %s' % seg.shape[0])
		# # print(seg.head().to_string())
		# # print('*****')

		os.remove('%s_segment_genes.txt' % self.case)
		os.remove('%s_ratio_genes.txt' % self.case)
		os.remove('%s_trusted_genes.txt' % self.case)
		
		seg.to_csv(self.output()[5].path, sep='\t', header=True, index=False)