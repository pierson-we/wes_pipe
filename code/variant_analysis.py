#!/usr/bin/env python3
import subprocess
import luigi
import luigi.interface
import os
import sys
import time
import itertools
import gzip
import pipeline_utils
import global_vars
import bam_processing
import variant_calling

class vep(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# vcf_path = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()
	# vcf_path = luigi.Parameter()

	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		# return bam_processing.aggregate_variants(case=self.case, tumor=self.tumor, matched_n=self.matched_n, case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads)
		return fpfilter(case=self.case, tumor=self.tumor, matched_n=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg, case_dict=self.case_dict)

	def output(self):
		# case_dir = os.path.join(self.project_dir, 'output', self.case)
		# vcf_path = os.path.join(case_dir, 'variants')
		# return [luigi.LocalTarget(os.path.join(vcf_path, self.case + '_vep.vcf')), luigi.LocalTarget(os.path.join(vcf_path, self.case + '_vep.vcf_summary.html'))]
		# return luigi.LocalTarget(self.vcf_path.split('fpfilter.vcf')[0] + 'vep.vcf')
		outputs = []
		for fpfilter_vcf in self.input():
			fpfilter_path = fpfilter_vcf.path
			vep_path = fpfilter_path.split('fpfilter.vcf')[0] + 'vep.vcf'
			outputs.append(luigi.LocalTarget(vep_path))
		return outputs

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		
		for fpfilter_vcf in self.input():
			fpfilter_path = fpfilter_vcf.path
			vep_path = fpfilter_path.split('fpfilter.vcf')[0] + 'vep.vcf'

			cmd = ['./packages/ensembl-vep/vep', '-i', fpfilter_path, '-o', vep_path, '--fasta', self.cfg['fasta_file'], '-fork', self.max_threads, '--cache', '--dir_cache', './packages/ensembl-vep/cache', '--protein', '--symbol', '--hgvs', '--force_overwrite', '--check_existing', '--offline', '--vcf'] #, '--buffer_size', '2500']
			pipeline_utils.command_call(cmd, self.output(), threads_needed=self.max_threads)

class fpfilter(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# vcf_path = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()
	# vcf_path = luigi.Parameter()

	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		# return bam_processing.aggregate_variants(case=self.case, tumor=self.tumor, matched_n=self.matched_n, case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads)
		return bam_processing.aggregate_variants(case=self.case, tumor=self.tumor, matched_n=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads, case_dict=self.case_dict, cfg=self.cfg)
		
	def output(self):
		# case_dir = os.path.join(self.project_dir, 'output', self.case)
		# vcf_path = os.path.join(case_dir, 'variants')
		# return luigi.LocalTarget(os.path.join(vcf_path, self.case + '_fpfilter' + '.vcf'))
		# return luigi.LocalTarget(self.vcf_path.split('.vcf')[0] + '_fpfilter.vcf')

		outputs = []
		for variant_caller_output in self.input():
			if not isinstance(variant_caller_output, list):
				vcf_path = variant_caller_output.path.split('.vcf')[0]
				# print(vcf_path)
				outputs.append(luigi.LocalTarget(vcf_path + '_fpfilter' + '.vcf'))
		return outputs

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		
		tumor_bam = os.path.join(self.project_dir, 'output', self.case, 'alignment', self.case + '_T_recalibrated.bam')

		for variant_caller_output in self.input():
			if not isinstance(variant_caller_output, list):
				input_vcf = variant_caller_output.path
				output_vcf = input_vcf.split('.vcf')[0] + '_fpfilter.vcf'
				cmd = ['./packages/fpfilter/fpfilter.pl', '--vcf-file', input_vcf, '--bam-file', tumor_bam, '--reference', self.cfg['fasta_file'], '--sample', self.case + '_T', '--output', output_vcf]
				pipeline_utils.command_call(cmd, [luigi.LocalTarget(output_vcf)])

class msings_baseline(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# vcf_path = luigi.Parameter()
	# case = luigi.Parameter()
	# tumor = luigi.Parameter()
	# matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		# return bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads)
		return [bam_processing.recalibrated_bam(sample=case_name + '_N', fastq_file=self.case_dict[case_name]['N'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg) for case_name in self.case_dict if self.case_dict[case_name]['N'] != '']

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'normal_bams.txt')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'MSI_BASELINE.txt'))] \
		+ list(itertools.chain(*[[luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'msings', 'baseline', case_name + '_N_recalibrated', case_name + '_N_recalibrated.%s' % file_ext)) for file_ext in ['mpileup', 'msi_output', 'msi.txt']] for case_name in self.case_dict if self.case_dict[case_name]['N'] != '']))
	
	def run(self):
		print(self.output())
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		# cmd = ['./packages/msings/scripts/create_intervals.sh', './packages/MANTIS/b37_exome_microsatellites.bed']
		# pipeline_utils.command_call(cmd, self.output())

		normal_bams_file = os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'normal_bams.txt')
		with open(normal_bams_file, 'w') as f:
			normal_bams_list = [os.path.join(self.project_dir, 'output', case_name, 'alignment', case_name + '_N_recalibrated.bam') for case_name in self.case_dict if self.case_dict[case_name]['N'] != '']
			f.write('\n'.join(normal_bams_list))
		cmd = ['./packages/msings/scripts/create_baseline.sh', normal_bams_file, './packages/MANTIS/b37_1000_ms_loci_annotated_mono.msi_intervals', './packages/MANTIS/b37_1000_ms_loci_annotated_mono.bed', self.cfg['fasta_file'], os.path.join(self.project_dir, 'output', 'msings', 'baseline')]
		pipeline_utils.command_call(cmd, self.output())

class msi(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	vcf_path = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), bam_processing.recalibrated_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), msings_baseline(project_dir=self.project_dir, max_threads=self.max_threads, case_dict=self.case_dict, cfg=self.cfg)]
		else:
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), msings_baseline(project_dir=self.project_dir, max_threads=self.max_threads, case_dict=self.case_dict, cfg=self.cfg)]
			# return [msings_baseline(project_dir=self.project_dir, max_threads=self.max_threads, case_dict=self.case_dict)] #, vcf_path=self.vcf_path, case=self.case, tumor=self.tumor, matched_n=self.matched_n)]

	def output(self):
		if self.matched_n != '':
			return [luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msi.txt')), luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msi.txt.status')), luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msi.kmer_counts.txt')), luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msi.kmer_counts_filtered.txt')), luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msings.txt'))]
		else:
			return [luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msings.txt'))]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		if self.matched_n != '':
			cmd = ['python3', './packages/MANTIS/mantis.py', '-b', './packages/MANTIS/b37_1000_ms_loci_annotated.sorted.bed', '--genome', self.cfg['fasta_file'], '-t', self.input()[0][0].path, '-n', self.input()[1][0].path, '-mrq', '20.0', '-mlq', '25.0', '-mlc', '20', '-mrr', '1', '-o', self.output()[0].path]
			pipeline_utils.command_call(cmd, self.output())
		# else:
		# tumor_bams_file = os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'tumor_bams.txt')
		# with open(tumor_bams_file, 'w') as f:
		# 	tumor_bams_list = [os.path.join(self.project_dir, 'output', self.case, 'alignment', case_name + '_T_recalibrated.bam') for case_name in self.case_dict if self.case_dict[case_name]['N'] == '']
		# 	f.write('\n'.join(tumor_bams_list))
	
		cmd = ['./packages/msings/scripts/run_msings_single_sample.sh', self.input()[0][0].path, './packages/MANTIS/b37_1000_ms_loci_annotated_mono.msi_intervals', './packages/MANTIS/b37_1000_ms_loci_annotated_mono.bed', self.cfg['fasta_file'], os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'MSI_BASELINE.txt'), os.path.join(self.project_dir, 'output', 'msings', 'tumor')]
		pipeline_utils.command_call(cmd, self.output())
		os.rename(os.path.join(self.project_dir, 'output', 'msings', 'tumor', self.case + '_T_recalibrated', self.case + '_T_recalibrated.MSI_Analysis.txt'), os.path.join(self.vcf_path, self.case + '_msings.txt'))
		# cmd = ['echo', '"mSINGS still needs to be set up for tumor-only samples"', '>', self.output()[0].path] # this will be a pain to get up and going: https://bitbucket.org/uwlabmed/msings/src/8269e0e01acfc5e01d0de9d63ffc1e399996ce8a/Recommendations_for_custom_assays?at=master&fileviewer=file-view-default
		# pipeline_utils.command_call(cmd, self.output())

