#!/usr/bin/env python3
import subprocess
import luigi
import luigi.interface
import os
import sys
import time
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
		return luigi.LocalTarget(os.path.join(vcf_path, self.case + '_vep' + '.vcf'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		cmd = ['./package/ensembl-vep/vep', '-i', self.input()[0].path, '-o', self.output().path, '--fasta', self.fasta_file, '--cache', '--dir_cache', './packages/ensembl-vep/cache', '--protein', '--symbol', '--hgvs', '--force_overwrite', '--check_existing', '--offline']
		pipeline_utils.command_call(cmd, [self.output()])


class msi(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	vcf_path = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()

	fasta_file = luigi.Parameter()

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads), bam_processing.recalibrated_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads)]
		else:
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads)]

	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_msi' + '.txt'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		if self.matched_n != '':
			cmd = ['python', './packages/MANTIS/mantis.py', '--bed_file', './packages/MANTIS/b37_exome_microsatellites.bed', '--genome', self.fasta_file, '-t', self.input()[0].path, '-n', self.input()[1].path, '--threads', self.max_threads, '-mrq', '20.0', '-mlq', '25.0', '-mlc', '20', '-mrr', '1', '-o', self.output().path]
			pipeline_utils.command_call(cmd, [self.output()], threads_needed=self.max_threads)
		else:
			cmd = ['echo', '"mSINGS still needs to be set up for tumor-only samples"', '>', self.output().path] # this will be a pain to get up and going: https://bitbucket.org/uwlabmed/msings/src/8269e0e01acfc5e01d0de9d63ffc1e399996ce8a/Recommendations_for_custom_assays?at=master&fileviewer=file-view-default
			pipeline_utils.command_call(cmd, [self.output()])