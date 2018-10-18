#!/usr/bin/env python3
import subprocess
import luigi
import luigi.interface
import os
import sys
import time
import random
import gzip
import pipeline_utils
import global_vars
import variant_analysis
from variant_calling import *

# ***
# GATK Resource Bundle Files
# https://software.broadinstitute.org/gatk/download/bundle
# ***

# @luigi.Task.event_handler(luigi.Event.FAILURE)
# def error_handling(exception):
# 	print('Current working files at time of interruption:')
# 	print(global_vars.working_files)
# 	# print(cwd)
# 	os.chdir(global_vars.cwd)
# 	# for file in global_vars.working_files:
# 	# 	os.remove(file)
# 	raise exception

class genome_index(luigi.Task):
	max_threads = luigi.IntParameter()
	# fasta_file = luigi.Parameter()
	# threads = luigi.Parameter()
	# base_name = luigi.Parameter()
	# bowtie_build_location = luigi.Parameter()

	# fasta_dir = os.path.join(*luigi.Parameter().task_value('genome_index', 'fasta_file').split('/')[:-1])

	cfg = luigi.DictParameter()

	def output(self):
		# fasta_dir = os.path.join(*self.cfg['fasta_file'].split('/')[:-1])
		fasta_dir = '/'.join(self.cfg['fasta_file'].split('/')[:-1])
		return luigi.LocalTarget(os.path.join(fasta_dir, 'index', self.cfg['base_name'] + '.1.bt2'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		# cwd = os.getcwd()
		# if cwd.split('/')[-1] != 'wes_pipe':
		# print(cwd)
		# os.chdir(self.fasta_dir)
		# if not os.path.exists('./index'):
		# 	os.mkdir('./index')
		# os.chdir('./index')
		# cmd = [cwd + self.cfg['bowtie_location'] + 'bowtie2-build', '--threads=%s' % self.max_threads, self.cfg['fasta_file'], self.cfg['base_name']]
		# fasta_dir = os.path.join(self.cfg['fasta_file'].split('/')[:-1])
		fasta_dir = '/'.join(self.cfg['fasta_file'].split('/')[:-1])

		cmd = [self.cfg['bowtie_build_location'], '--threads=%s' % self.max_threads, self.cfg['fasta_file'], os.path.join(fasta_dir, 'index', self.cfg['base_name'])]
		pipeline_utils.command_call(cmd, [self.output()], threads_needed=self.max_threads, sleep_time=0.1)
		# subprocess.call([self.cfg['bowtie_location'] + 'bowtie2-build', '--threads=%s' % self.max_threads, self.cfg['fasta_file'], self.cfg['base_name']], stdout=subprocess.PIPE)
		os.chdir(global_vars.cwd)

class samtools_index(luigi.Task):
	max_threads = luigi.IntParameter()
	# fasta_file = luigi.Parameter()
	# samtools_location = luigi.Parameter()

	cfg = luigi.DictParameter()

	def output(self):
		return luigi.LocalTarget(self.cfg['fasta_file'] + '.fai')
	
	def run(self):
		# cmd = [self.cfg['samtools_location'], 'faidx', '--nthreads=%s' % self.max_threads, self.cfg['fasta_file']]
		cmd = [self.cfg['samtools_location'], 'faidx', '--nthreads=%s' % self.max_threads, self.cfg['fasta_file']]
		pipeline_utils.command_call(cmd, [self.output()], threads_needed=self.max_threads, sleep_time=0.1)
		# samtools faidx Homo_sapiens_assembly18.fasta

class picard_index(luigi.Task):
	# fasta_file = luigi.Parameter()
	# picard_location = luigi.Parameter()

	cfg = luigi.DictParameter()

	def output(self):
		return luigi.LocalTarget(self.cfg['fasta_file'] + '.dict')
	
	def run(self):
		cmd = ['java', '-jar', self.cfg['picard_location'], 'CreateSequenceDictionary', 'R=%s' % self.cfg['fasta_file'], 'O=%s' % self.output().path]
		pipeline_utils.command_call(cmd, [self.output()], sleep_time=0.1)

class bwa_index(luigi.Task):
	cfg = luigi.DictParameter()

	def output(self):
		return [luigi.LocalTarget(self.cfg['fasta_file'] + '.amb'), luigi.LocalTarget(self.cfg['fasta_file'] + '.ann'), luigi.LocalTarget(self.cfg['fasta_file'] + '.bwt'), luigi.LocalTarget(self.cfg['fasta_file'] + '.pac'), luigi.LocalTarget(self.cfg['fasta_file'] + '.sa')]
	
	def run(self):
		cmd = ['bwa', 'index', '-a', 'bwtsw', self.cfg['fasta_file']]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.1)

class trim(luigi.Task):
	fastq_file = luigi.Parameter()
	project_dir = luigi.Parameter()
	# trim_location = luigi.Parameter()
	sample = luigi.Parameter()

	cfg = luigi.DictParameter()

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.fastq_file.split('\t')[0].split('/')[-1].split('.')[0] + '_val_1.fq.gz')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.fastq_file.split('\t')[1].split('/')[-1].split('.')[0] + '_val_2.fq.gz'))]
	
	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		cmd = [self.cfg['trim_location'], '--paired', self.fastq_file.split('\t')[0], self.fastq_file.split('\t')[1], '-o', os.path.join(self.project_dir, 'output', self.sample[:-2])]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.05)

class fastqc(luigi.Task):
	sample = luigi.Parameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	both_fastq_files = luigi.Parameter()
	# fastqc_location = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return trim(fastq_file=self.both_fastq_files, sample=self.sample, project_dir=self.project_dir, cfg=self.cfg)

	def output(self):
		return [luigi.LocalTarget(self.fastq_file), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'fastqc', self.fastq_file.split('/')[-1].split('.')[0] + '_fastqc.html'))]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		pipeline_utils.confirm_path(self.output()[1].path)
		cmd = [self.cfg['fastqc_location'], '--outdir=%s' % os.path.join(self.project_dir, 'output', self.sample[:-2], 'fastqc'), os.path.join(self.project_dir, 'output', self.sample[:-2], self.fastq_file)]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.1)

class fastqc_launch(luigi.Task):
	sample = luigi.Parameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return [trim(fastq_file=self.fastq_file, sample=self.sample, project_dir=self.project_dir, cfg=self.cfg), 
		fastqc(fastq_file=os.path.join(self.project_dir, 'output', self.sample[:-2], self.fastq_file.split('\t')[0].split('/')[-1].split('.')[0] + '_val_1.fq.gz'), sample=self.sample, project_dir=self.project_dir, both_fastq_files=self.fastq_file, cfg=self.cfg),
		fastqc(fastq_file=os.path.join(self.project_dir, 'output', self.sample[:-2], self.fastq_file.split('\t')[1].split('/')[-1].split('.')[0] + '_val_2.fq.gz'), sample=self.sample, project_dir=self.project_dir, both_fastq_files=self.fastq_file, cfg=self.cfg)]

	def output(self):
		return [self.input()[1], self.input()[2]]

	def run(self):
		pass

class bowtie(luigi.Task):
	max_threads = luigi.IntParameter()
	fastq_file = luigi.Parameter()
	# sam_file = luigi.Parameter()
	# threads = luigi.Parameter()
	project_dir = luigi.Parameter()
	# bowtie_location = luigi.Parameter()
	# samtools_location = luigi.Parameter()
	# fasta_file = luigi.Parameter()
	# base_name = luigi.Parameter()
	sample = luigi.Parameter()

	cfg = luigi.DictParameter()

	# def on_failure(self, exception):
	# 	error_handling(self, exception)

	# fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	cfg = luigi.DictParameter()

	def requires(self):
		return [genome_index(max_threads=self.max_threads, cfg=self.cfg), #threads=self.threads, base_name=self.cfg['base_name'], fasta_path=self.fasta_path)
		samtools_index(max_threads=self.max_threads, cfg=self.cfg),
		picard_index(cfg=self.cfg),
		fastqc_launch(fastq_file=self.fastq_file, sample=self.sample, project_dir=self.project_dir, cfg=self.cfg)]
		#fastqc(fastq_file=self.fastq_file.split('\t')[1], sample=self.sample, project_dir=self.project_dir)]

	def output(self):
		# print(os.path.join(os.path.join('/', *self.cfg['fasta_file'].split('/')[:-1]), 'index', self.sample + '_raw.sam'))
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_raw.bam'))

	def run(self):
		# try:
		# cwd = os.getcwd()
		# print(cwd)

		# os.chdir(os.path.join(os.path.join(*self.cfg['fasta_file'].split('/')[:-1]), 'index'))
		pipeline_utils.confirm_path(self.output().path)
		# os.chdir(os.path.join(self.fasta_dir, 'index'))
		# print(os.getcwd())
		# cmd = [os.path.join(cwd, self.cfg['bowtie_location'], 'bowtie2'), '-x', self.cfg['base_name'], '--threads=%s' % self.max_threads, '-U', self.fastq_file, '-S', self.sample + '_raw.sam']
		fasta_dir = os.path.join(*self.cfg['fasta_file'].split('/')[:-1])

		cmd = [self.cfg['bowtie_location'], '-x', os.path.join(fasta_dir, 'index', self.cfg['base_name']), '-1', self.input()[-1][0][0].path, '-2', self.input()[-1][1][0].path, '-p', self.max_threads, '--very-sensitive-local' '|', self.cfg['samtools_location'], 'view', '-bh', '-', '>', self.output().path]
		pipeline_utils.command_call(cmd, [self.output()], cwd=os.getcwd(), threads_needed=self.max_threads, sleep_time=0.2)

		# for input_file in self.input()[-1][0]:
		# 	input_file.remove()

		# os.chdir(global_vars.cwd)
		# print(os.getcwd())
		# except KeyboardInterrupt:
		# 	pipeline_utils.error_handling(KeyboardInterrupt)
		# 	os.chdir(cwd)
		# 	self.output().remove()
		# 	sys.exit()

class bwa(luigi.Task):
	max_threads = luigi.IntParameter()
	fastq_file = luigi.Parameter()
	# sam_file = luigi.Parameter()
	# threads = luigi.Parameter()
	project_dir = luigi.Parameter()
	# bowtie_location = luigi.Parameter()
	# samtools_location = luigi.Parameter()
	# fasta_file = luigi.Parameter()
	# base_name = luigi.Parameter()
	sample = luigi.Parameter()

	cfg = luigi.DictParameter()

	# def on_failure(self, exception):
	# 	error_handling(self, exception)

	# fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	def requires(self):
		return [bwa_index(cfg=self.cfg), #threads=self.threads, base_name=self.cfg['base_name'], fasta_path=self.fasta_path)
		samtools_index(max_threads=self.max_threads, cfg=self.cfg),
		picard_index(cfg=self.cfg),
		fastqc_launch(fastq_file=self.fastq_file, sample=self.sample, project_dir=self.project_dir, cfg=self.cfg)]
		#fastqc(fastq_file=self.fastq_file.split('\t')[1], sample=self.sample, project_dir=self.project_dir)]

	def output(self):
		# print(os.path.join(os.path.join('/', *self.cfg['fasta_file'].split('/')[:-1]), 'index', self.sample + '_raw.sam'))
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_raw.bam'))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		r1 = self.fastq_file.split('\t')[0]
		r2 = self.fastq_file.split('\t')[1]
		sai1 = self.sample + '_1.sai'
		sai2 = self.sample + '_2.sai'
		cmd = ['bwa aln %s %s > %s' % (self.cfg['fasta_file'], r1, sai1)]
		pipeline_utils.command_call(cmd, [self.output()], threads_needed=self.max_threads, sleep_time=0.2)
		cmd = ['bwa aln %s %s > %s' % (self.cfg['fasta_file'], r2, sai2)]
		pipeline_utils.command_call(cmd, [self.output()], threads_needed=self.max_threads, sleep_time=0.2)
		cmd = ['bwa sampe %s %s %s | samtools view -bh | samtools sort -o %s' % (self.cfg['fasta_file'], sai1, sai2, self.output().path)]
		

# ~45 mins
class add_read_groups(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	# picard_location = luigi.Parameter()
	# library_prep = luigi.Parameter()
	# platform = luigi.Parameter()

	# input_bam = luigi.Parameter()
	# output_bam = luigi.Parameter()
	sample = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return bowtie(fastq_file=self.fastq_file, sample=self.sample, max_threads=self.max_threads, project_dir=self.project_dir, cfg=self.cfg)

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_preprocessed.bam')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_preprocessed.bai'))]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		cmd = ['java', '-Xmx8g', '-Xms8g', '-XX:+UseSerialGC', '-Djava.io.tmpdir=%s' % self.cfg['tmp_dir'], '-jar', self.cfg['picard_location'], 'AddOrReplaceReadGroups', 'I=%s' % self.input().path, 'O=%s' % self.output()[0].path, 'CREATE_INDEX=true', 'SORT_ORDER=coordinate', 'RGID=%s' % self.sample, 'RGLB=%s' % self.cfg['library_prep'], 'RGPL=%s' % self.cfg['platform'], 'RGPU=%s' % self.sample + '_barcode', 'RGSM=%s' % self.sample]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.3)
		# self.input().remove()

# https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
# ~35 minutes
class mark_duplicates(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	# picard_location = luigi.Parameter()

	# input_bam = luigi.Parameter() #replace all of the input bams with the outputs of the previous
	# output_bam = luigi.Parameter()
	# metrics_file = luigi.Parameter()
	sample = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return add_read_groups(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads, cfg=self.cfg)

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_marked_dups.bam')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'tmp', self.sample + '_marked_dups_metrics.txt'))]

	def run(self):
		pipeline_utils.confirm_path(self.output()[0].path)
		pipeline_utils.confirm_path(self.output()[1].path)
		cmd = ['java', '-Xmx8g', '-Xms8g', '-XX:+UseSerialGC', '-Djava.io.tmpdir=%s' % self.cfg['tmp_dir'], '-jar', self.cfg['picard_location'], 'MarkDuplicates', 'I=%s' % self.input()[0].path, 'O=%s' % self.output()[0].path, 'M=%s' % self.output()[1].path, 'CREATE_INDEX=true', 'REMOVE_DUPLICATES=true', 'ASSUME_SORT_ORDER=coordinate']
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.4)
		# for input_file in self.input():
		# 	input_file.remove()

class index_bam(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	# samtools_location = luigi.Parameter()

	sample = luigi.Parameter()

	# input_bam = luigi.Parameter()
	# recalibration_report = luigi.Parameter()
	# output_bam = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return mark_duplicates(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads, cfg=self.cfg)

	def output(self):

		return [self.input()[0], luigi.LocalTarget(self.input()[0].path + '.bai')]

	def run(self):
		pipeline_utils.confirm_path(self.output()[1].path)
		# cmd = [os.getcwd() + '/' + self.cfg['samtools_location'], 'index', '-b', self.input()[0].path]
		cmd = [self.cfg['samtools_location'], 'index', '-b', self.input()[0].path]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.5)
		# self.input()[1].remove()

# ~20 mins w/2 cores
# This might be useless: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5048557/
class realigner_target(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	# gatk3_location = luigi.Parameter()

	sample = luigi.Parameter()
	# fasta_file = luigi.Parameter()
	# mills = luigi.Parameter()
	# kg = luigi.Parameter()

	# input_bam = luigi.Parameter()
	# output_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		# return index_bam(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads, cfg=self.cfg)
		return mark_duplicates(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads, cfg=self.cfg)

	def output(self):
		return self.input() + [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_forIndelRealigner.intervals'))]

	def run(self):
		pipeline_utils.confirm_path(self.output()[2].path)
		cmd = ['java', '-Xmx8g', '-Xms8g', '-XX:+UseSerialGC', '-Djava.io.tmpdir=%s' % self.cfg['tmp_dir'], '-jar', self.cfg['gatk3_location'], '-T', 'RealignerTargetCreator', '-nt', str(self.max_threads), '-R', self.cfg['fasta_file'], '-I', self.input()[0].path, '--known', self.cfg['mills'], '--known', self.cfg['kg'], '-o', self.output()[2].path]
		pipeline_utils.command_call(cmd, self.output(), threads_needed=self.max_threads, sleep_time=0.6)
		# for input_file in self.input():
		# 	input_file.remove()

# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php
# ~45 mins (no multiprocessing available)
# This might be useless: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5048557/
class indel_realignment(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	# gatk3_location = luigi.Parameter()

	# fasta_file = luigi.Parameter()
	# mills = luigi.Parameter()
	# kg = luigi.Parameter()
	sample = luigi.Parameter()
	
	# input_bam = luigi.Parameter()
	# output_bam = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return realigner_target(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads, cfg=self.cfg)

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_realigned.bam')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_realigned.bai'))]

	def run(self):
		pipeline_utils.confirm_path(self.output()[0].path)
		pipeline_utils.confirm_path(self.output()[1].path)
		cmd = ['java', '-Xmx8g', '-Xms8g', '-XX:+UseSerialGC', '-Djava.io.tmpdir=%s' % self.cfg['tmp_dir'], '-jar', self.cfg['gatk3_location'], '-T', 'IndelRealigner', '-R', self.cfg['fasta_file'], '-I', self.input()[0].path, '-known', self.cfg['mills'], '-known', self.cfg['kg'], '-targetIntervals', self.input()[2].path, '-o', self.output()[0].path]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.7)
		# for input_file in self.input():
		# 	input_file.remove()

# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
# ~45 mins (no multiprocessing available)
# This might be useless: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5048557/
class bqsr(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	# gatk4_location = luigi.Parameter()

	# fasta_file = luigi.Parameter()
	# mills = luigi.Parameter()
	# kg = luigi.Parameter()
	# dbsnp = luigi.Parameter()
	sample = luigi.Parameter()
	
	# input_bam = luigi.Parameter()
	# output_table = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return mark_duplicates(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads, cfg=self.cfg)

	def output(self):
		return self.input() + [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_recalibrated.table'))]

	def run(self):
		pipeline_utils.confirm_path(self.output()[2].path)
		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx8g -Xms8g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'BaseRecalibrator', '-R', self.cfg['fasta_file'], '-I', self.input()[0].path, '--known-sites', self.cfg['known_vcf'], '--known-sites', self.cfg['mills'], '--known-sites', self.cfg['kg'], '-O',  self.output()[2].path]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.8)
		# for input_file in self.input():
		# 	input_file.remove()
		# self.input().remove()

# ~75 mins (no multiprocessing available)
class recalibrated_bam(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	# gatk4_location = luigi.Parameter()

	# fasta_file = luigi.Parameter()
	sample = luigi.Parameter()

	cfg = luigi.DictParameter()

	# input_bam = luigi.Parameter()
	# recalibration_report = luigi.Parameter()
	# output_bam = luigi.Parameter()

	def requires(self):
		return bqsr(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads, cfg=self.cfg)

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_recalibrated.bam')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_recalibrated.bam.bai'))]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx8g -Xms8g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'ApplyBQSR', '-R', self.cfg['fasta_file'], '-I', self.input()[0].path, '--bqsr-recal-file', self.input()[2].path, '-O',  self.output()[0].path]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.9)
		cmd = ['mv', luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_recalibrated.bai')).path, self.output()[1].path]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.9)
		# for input_file in self.input():
		# 	input_file.remove()

class aggregate_variants(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()

	cfg = luigi.DictParameter()
	
	def requires(self):
		case_dir = os.path.join(self.project_dir, 'output', self.case)
		vcf_path = os.path.join(case_dir, 'variants')

		kwargs = {
		'project_dir': self.project_dir,
		# 'case_dir': case_dir,
		'vcf_path': vcf_path,
		'case': self.case,
		'tumor': self.tumor,
		'matched_n': self.matched_n,
		'max_threads': self.max_threads,
		'cfg': self.cfg
		}
		return [filter_mutect(case_dict=self.case_dict, **kwargs),
		vardict(**kwargs)
		# freebayes(**kwargs),
		# scalpel_export(**kwargs),
		#variant_analysis.msi(case_dict=self.case_dict, **kwargs)
		# yield run_variant_caller(caller='VarDict', **kwargs)
		# yield run_variant_caller(caller='FreeBayes', **kwargs)
		# yield run_variant_caller(caller='VarScan', **kwargs)

		# yield run_variant_caller(caller='Scalpel', **kwargs)

		# yield run_variant_caller(caller='CNVkit', **kwargs)

		# yield run_variant_caller(caller='LUMPY', **kwargs)
		# yield run_variant_caller(caller='DELLY', **kwargs)
		# yield run_variant_caller(caller='WHAM', **kwargs)
		]

	def output(self):
		# outputs = []
		# for variant_caller_output in self.input():
		# 	if not isinstance(variant_caller_output, list):
		# 		vcf_path = variant_caller_output.path.split('.vcf')[0]
		# 		# print(vcf_path)
		# 		for vcf_filter in ['fpfilter', 'vep']:
		# 			outputs.append(luigi.LocalTarget(vcf_path + '_' + vcf_filter + '.vcf'))
		# 	else:
		# 		outputs += variant_caller_output
		# return outputs
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'all_samples.vcf'))

	def run(self):
		input_vcfs = []
		for variant_caller in self.input():
			if not isinstance(variant_caller, list):
				vcf_path = variant_caller.path
				if '.gz' in vcf_path:
					with gzip.open(vcf_path, 'rb') as f:
						with open(vcf_path.split('.gz')[0], 'wb') as new_f:
							new_f.write(f.read())
					vcf_path = vcf_path.split('.gz')[0]
				input_vcfs.append(vcf_path)
		cmd = ['java', '-jar', self.cfg['gatk3_location'], '-T', 'CombineVariants', '-R', self.cfg['fasta_file', '-o', self.output().path]] + ['--variant']
		
	# 			variant_analysis.vep(max_threads=self.max_threads, project_dir=self.project_dir, case=self.case, tumor=self.tumor, matched_n=self.matched_n, vcf_path=vcf_path, cfg=self.cfg)

# class filter_variants(luigi.Task):
# 	max_threads = luigi.IntParameter()
# 	project_dir = luigi.Parameter()
# 	case = luigi.Parameter()
# 	tumor = luigi.Parameter()
# 	matched_n = luigi.Parameter()
# 	case_dict = luigi.DictParameter()

# 	cfg = luigi.DictParameter()
	
# 	def requires(self):
# 		# requirements = []
# 		# for variant_caller in self.input():
# 		# 	if not isinstance(variant_caller, list):
# 		# 		vcf_path = variant_caller.path
# 		# 		if '.gz' in vcf_path:
# 		# 			with gzip.open(vcf_path, 'rb') as f:
# 		# 				with open(vcf_path.split('.gz')[0], 'wb') as new_f:
# 		# 					new_f.write(f.read())
# 		# 			vcf_path = vcf_path.split('.gz')[0]
# 		# 		requirements.append(variant_analysis.vep(max_threads=self.max_threads, project_dir=self.project_dir, case=self.case, tumor=self.tumor, matched_n=self.matched_n, vcf_path=vcf_path, cfg=self.cfg))

# 		return aggregate_variants(case=self.case, tumor=self.tumor, matched_n=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads, case_dict=self.case_dict, cfg=self.cfg)
# 		# + requirements
		
# 		# yield run_variant_caller(caller='VarDict', **kwargs)
# 		# yield run_variant_caller(caller='FreeBayes', **kwargs)
# 		# yield run_variant_caller(caller='VarScan', **kwargs)

# 		# yield run_variant_caller(caller='Scalpel', **kwargs)

# 		# yield run_variant_caller(caller='CNVkit', **kwargs)

# 		# yield run_variant_caller(caller='LUMPY', **kwargs)
# 		# yield run_variant_caller(caller='DELLY', **kwargs)
# 		# yield run_variant_caller(caller='WHAM', **kwargs)

# 	def output(self):
# 		# outputs = []
# 		# for variant_caller_output in self.input():
# 		# 	if not isinstance(variant_caller_output, list):
# 		# 		vcf_path = variant_caller_output.path.split('.vcf')[0]
# 		# 		# print(vcf_path)
# 		# 		for vcf_filter in ['fpfilter', 'vep']:
# 		# 			outputs.append(luigi.LocalTarget(vcf_path + '_' + vcf_filter + '.vcf'))
# 		# 	else:
# 		# 		outputs += variant_caller_output
# 		# return outputs

# 		# outputs = []
# 		# for variant_caller in self.input():
# 		# 	if not isinstance(variant_caller, list):
# 		# 		vcf_path = variant_caller.path
# 		# 		if '.gz' in vcf_path:
# 		# 			with gzip.open(vcf_path, 'rb') as f:
# 		# 				with open(vcf_path.split('.gz')[0], 'wb') as new_f:
# 		# 					new_f.write(f.read())
# 		# 			vcf_path = vcf_path.split('.gz')[0]
# 		# 		outputs.append(variant_analysis.vep(max_threads=self.max_threads, project_dir=self.project_dir, case=self.case, tumor=self.tumor, matched_n=self.matched_n, vcf_path=vcf_path, cfg=self.cfg))

# 		return luigi.LocalTarget()

# 	def run(self):
# 		for output in self.output():
# 			pipeline_utils.confirm_path(output.path)
# 		for variant_caller in self.input():
# 			if not isinstance(variant_caller, list):
# 				vcf_path = variant_caller.path
# 				if '.gz' in vcf_path:
# 					with gzip.open(vcf_path, 'rb') as f:
# 						with open(vcf_path.split('.gz')[0], 'wb') as new_f:
# 							new_f.write(f.read())
# 					vcf_path = vcf_path.split('.gz')[0]
# 				variant_analysis.vep(max_threads=self.max_threads, project_dir=self.project_dir, case=self.case, tumor=self.tumor, matched_n=self.matched_n, vcf_path=vcf_path, cfg=self.cfg, input_vcf=vcf_path, output_vcf=self.output.path ,output=self.output())


	

	
	
