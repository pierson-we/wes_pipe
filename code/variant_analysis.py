#!/usr/bin/env python3
import subprocess
import luigi
import luigi.interface
import os
import sys
import time
import itertools
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

	fasta_file = luigi.Parameter()

	def requires(self):
		return bam_processing.aggregate_variants(case=self.case, tumor=self.tumor, matched_n=self.matched_n, case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads)

	def output(self):
		case_dir = os.path.join(self.project_dir, 'output', self.case)
		vcf_path = os.path.join(case_dir, 'variants')
		return [luigi.LocalTarget(os.path.join(vcf_path, self.case + '_vep.vcf')), luigi.LocalTarget(os.path.join(vcf_path, self.case + '_vep.vcf_summary.html'))]
	
	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		cmd = ['./packages/ensembl-vep/vep', '-i', self.input()[0].path, '-o', self.output()[0].path, '--fasta', self.fasta_file, '--cache', '--dir_cache', './packages/ensembl-vep/cache', '--protein', '--symbol', '--hgvs', '--force_overwrite', '--check_existing', '--offline', '--buffer_size', '2500']
		pipeline_utils.command_call(cmd, self.output())

class fpfilter(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# vcf_path = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()

	fasta_file = luigi.Parameter()

	def requires(self):
		return bam_processing.aggregate_variants(case=self.case, tumor=self.tumor, matched_n=self.matched_n, case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads)

	def output(self):
		case_dir = os.path.join(self.project_dir, 'output', self.case)
		vcf_path = os.path.join(case_dir, 'variants')
		return luigi.LocalTarget(os.path.join(vcf_path, self.case + '_snvs' + '.fpfilter'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		fpfilter_path = os.path.join(self.project_dir, 'output', self.case, 'fpfilter')
		pipeline_utils.confirm_path(fpfilter_path)
		cmd = ['mkdir', fpfilter_path]
		pipeline_utils.command_call(cmd, [self.output()])
		snvs_var = os.path.join(fpfilter_path, 'snvs.var')
		cmd = ['perl', '-ane', '''\'print join("\t",@F[0,1,1])."\n" unless(m/^#/)\'''', self.input()[0].path, '>', snvs_var]
		pipeline_utils.command_call(cmd, [self.output()])
		snvs_readcount = os.path.join(fpfilter_path, 'snvs.readcount')
		tumor_bam = os.path.join(self.project_dir, 'output', self.case, 'alignment', self.case + '_T_recalibrated.bam')
		cmd = ['./packages/fpfilter/bam-readcount-master/build/bin/bam-readcount', '-q1', '-b15', '-w1', '-l', snvs_var, '-f', self.fasta_file, tumor_bam, '>', snvs_readcount]
		pipeline_utils.command_call(cmd, [self.output()])
		cmd = ['./packages/fpfilter/fpfilter.pl', '--var-file', self.input()[0].path, '--readcount-file', snvs_readcount, '--output-file', self.output().path]
		pipeline_utils.command_call(cmd, [self.output()])

class msings_baseline(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	vcf_path = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()

	fasta_file = luigi.Parameter()

	def requires(self):
		return bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads)

	def output(self):
		return self.input() + [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'normal_bams.txt')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'MSI_BASELINE.txt'))] \
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
		cmd = ['source', './packages/msings/scripts/create_baseline.sh', normal_bams_file, './packages/MANTIS/b37_exome_microsatellites.msi_intervals', './packages/MANTIS/b37_exome_microsatellites.bed', self.fasta_file, os.path.join(self.project_dir, 'output', 'msings', 'baseline')]
		pipeline_utils.command_call(cmd, self.output())

class msi(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	vcf_path = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()

	fasta_file = luigi.Parameter()

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads), bam_processing.recalibrated_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads)]
		else:
			return [msings_baseline(project_dir=self.project_dir, vcf_path=self.vcf_path, case=self.case, tumor=self.tumor, matched_n=self.matched_n, max_threads=self.max_threads, case_dict=self.case_dict)]

	def output(self):
		if self.matched_n != '':
			return [luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msi.txt')), luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msi.txt.status')), luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msi.kmer_counts.txt')), luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msi.kmer_counts_filtered.txt'))]
		else:
			return [luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msi.txt'))]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		if self.matched_n != '':
			cmd = ['python3', './packages/MANTIS/mantis.py', '-b', './packages/MANTIS/b37_exome_microsatellites.bed', '--genome', self.fasta_file, '-t', self.input()[0][0].path, '-n', self.input()[1][0].path, '--threads', self.max_threads, '-mrq', '20.0', '-mlq', '25.0', '-mlc', '20', '-mrr', '1', '-o', self.output()[0].path]
			pipeline_utils.command_call(cmd, self.output(), threads_needed=self.max_threads)
		# else:
		# tumor_bams_file = os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'tumor_bams.txt')
		# with open(tumor_bams_file, 'w') as f:
		# 	tumor_bams_list = [os.path.join(self.project_dir, 'output', self.case, 'alignment', case_name + '_T_recalibrated.bam') for case_name in self.case_dict if self.case_dict[case_name]['N'] == '']
		# 	f.write('\n'.join(tumor_bams_list))
	
		cmd = ['source', './packages/msings/scripts/run_msings_single_sample.sh', self.input()[0][0].path, './packages/MANTIS/b37_exome_microsatellites.msi_intervals', './packages/MANTIS/b37_exome_microsatellites.bed', self.fasta_file, os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'MSI_BASELINE.txt'), os.path.join(self.project_dir, 'output', 'msings', 'tumor')]
		pipeline_utils.command_call(cmd, self.output())
		os.rename(os.path.join(self.project_dir, 'output', 'msings', 'tumor', self.case + '_T_recalibrated', self.case + '_T_recalibrated.MSI_Analysis.txt'), os.path.join(self.vcf_path, self.case + '_msi.txt'))
		# cmd = ['echo', '"mSINGS still needs to be set up for tumor-only samples"', '>', self.output()[0].path] # this will be a pain to get up and going: https://bitbucket.org/uwlabmed/msings/src/8269e0e01acfc5e01d0de9d63ffc1e399996ce8a/Recommendations_for_custom_assays?at=master&fileviewer=file-view-default
		# pipeline_utils.command_call(cmd, self.output())

