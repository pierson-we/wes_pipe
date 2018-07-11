#!/usr/bin/env python3
import subprocess
import luigi
import luigi.interface
import os
import sys
import time
import random
import pipeline_utils
import global_vars
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
	fasta_file = luigi.Parameter()
	# threads = luigi.Parameter()
	base_name = luigi.Parameter()
	bowtie_location = luigi.Parameter()

	fasta_dir = os.path.join(*luigi.Parameter().task_value('genome_index', 'fasta_file').split('/')[:-1])

	def output(self):
		return luigi.LocalTarget(os.path.join(self.fasta_dir, 'index', self.base_name + '.1.bt2'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		# cwd = os.getcwd()
		# if cwd.split('/')[-1] != 'wes_pipe':
		# print(cwd)
		# os.chdir(self.fasta_dir)
		# if not os.path.exists('./index'):
		# 	os.mkdir('./index')
		# os.chdir('./index')
		# cmd = [cwd + self.bowtie_location + 'bowtie2-build', '--threads=%s' % self.max_threads, self.fasta_file, self.base_name]
		cmd = [self.bowtie_location, '--threads=%s' % self.max_threads, self.fasta_file, os.path.join(self.fasta_dir, 'index', self.base_name)]
		pipeline_utils.command_call(cmd, [self.output()], threads_needed=self.max_threads, sleep_time=0.1)
		# subprocess.call([self.bowtie_location + 'bowtie2-build', '--threads=%s' % self.max_threads, self.fasta_file, self.base_name], stdout=subprocess.PIPE)
		os.chdir(global_vars.cwd)

class samtools_index(luigi.Task):
	max_threads = luigi.IntParameter()
	fasta_file = luigi.Parameter()
	samtools_location = luigi.Parameter()

	def output(self):
		return luigi.LocalTarget(self.fasta_file + '.fai')
	
	def run(self):
		# cmd = [self.samtools_location, 'faidx', '--nthreads=%s' % self.max_threads, self.fasta_file]
		cmd = [self.samtools_location, 'faidx', '--nthreads=%s' % self.max_threads, self.fasta_file]
		pipeline_utils.command_call(cmd, [self.output()], threads_needed=self.max_threads, sleep_time=0.1)
		# samtools faidx Homo_sapiens_assembly18.fasta

class picard_index(luigi.Task):
	fasta_file = luigi.Parameter()
	picard_location = luigi.Parameter()

	def output(self):
		return luigi.LocalTarget(self.fasta_file + '.dict')
	
	def run(self):
		cmd = ['java', '-jar', self.picard_location, 'CreateSequenceDictionary', 'R=%s' % self.fasta_file, 'O=%s' % self.output().path]
		pipeline_utils.command_call(cmd, [self.output()], sleep_time=0.1)

class trim(luigi.Task):
	fastq_file = luigi.Parameter()
	project_dir = luigi.Parameter()
	trim_location = luigi.Parameter()
	sample = luigi.Parameter()

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.fastq_file.split('\t')[0].split('/')[-1].split('.')[0] + '_val_1.fq.gz')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.fastq_file.split('\t')[1].split('/')[-1].split('.')[0] + '_val_2.fq.gz'))]
	
	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		cmd = [self.trim_location, '--paired', self.fastq_file.split('\t')[0], self.fastq_file.split('\t')[1], '-o', os.path.join(self.project_dir, 'output', self.sample[:-2])]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.05)

class fastqc(luigi.Task):
	sample = luigi.Parameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	both_fastq_files = luigi.Parameter()
	fastqc_location = luigi.Parameter()

	def requires(self):
		return trim(fastq_file=self.both_fastq_files, sample=self.sample, project_dir=self.project_dir)

	def output(self):
		return [luigi.LocalTarget(self.fastq_file), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'fastqc', self.fastq_file.split('/')[-1].split('.')[0] + '_fastqc.html'))]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		pipeline_utils.confirm_path(self.output()[1].path)
		cmd = [self.fastqc_location, '--outdir=%s' % os.path.join(self.project_dir, 'output', self.sample[:-2], 'fastqc'), os.path.join(self.project_dir, 'output', self.sample[:-2], self.fastq_file)]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.1)

class fastqc_launch(luigi.Task):
	sample = luigi.Parameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()

	def requires(self):
		return [trim(fastq_file=self.fastq_file, sample=self.sample, project_dir=self.project_dir), 
		fastqc(fastq_file=os.path.join(self.project_dir, 'output', self.sample[:-2], self.fastq_file.split('\t')[0].split('/')[-1].split('.')[0] + '_val_1.fq.gz'), sample=self.sample, project_dir=self.project_dir, both_fastq_files=self.fastq_file),
		fastqc(fastq_file=os.path.join(self.project_dir, 'output', self.sample[:-2], self.fastq_file.split('\t')[1].split('/')[-1].split('.')[0] + '_val_2.fq.gz'), sample=self.sample, project_dir=self.project_dir, both_fastq_files=self.fastq_file)]

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
	bowtie_location = luigi.Parameter()
	samtools_location = luigi.Parameter()
	fasta_file = luigi.Parameter()
	base_name = luigi.Parameter()
	sample = luigi.Parameter()

	fasta_dir = os.path.join(*luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	# def on_failure(self, exception):
	# 	error_handling(self, exception)

	# fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	def requires(self):
		return [genome_index(max_threads=self.max_threads), #threads=self.threads, base_name=self.base_name, fasta_path=self.fasta_path)
		samtools_index(max_threads=self.max_threads),
		picard_index(),
		fastqc_launch(fastq_file=self.fastq_file, sample=self.sample, project_dir=self.project_dir)]
		#fastqc(fastq_file=self.fastq_file.split('\t')[1], sample=self.sample, project_dir=self.project_dir)]

	def output(self):
		# print(os.path.join(os.path.join('/', *self.fasta_file.split('/')[:-1]), 'index', self.sample + '_raw.sam'))
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_raw.bam'))

	def run(self):
		# try:
		# cwd = os.getcwd()
		# print(cwd)

		# os.chdir(os.path.join(os.path.join(*self.fasta_file.split('/')[:-1]), 'index'))
		pipeline_utils.confirm_path(self.output().path)
		# os.chdir(os.path.join(self.fasta_dir, 'index'))
		# print(os.getcwd())
		# cmd = [os.path.join(cwd, self.bowtie_location, 'bowtie2'), '-x', self.base_name, '--threads=%s' % self.max_threads, '-U', self.fastq_file, '-S', self.sample + '_raw.sam']
		cmd = [self.bowtie_location, '-x', os.path.join(self.fasta_dir, 'index', self.base_name), '-1', self.input()[-1][0][0].path, '-2', self.input()[-1][1][0].path, '-p', self.max_threads, '|', self.samtools_location, 'view', '-b', '-', '>', self.output().path]
		pipeline_utils.command_call(cmd, [self.output()], cwd=os.getcwd(), threads_needed=self.max_threads, sleep_time=0.2)

		# os.chdir(global_vars.cwd)
		# print(os.getcwd())
		# except KeyboardInterrupt:
		# 	pipeline_utils.error_handling(KeyboardInterrupt)
		# 	os.chdir(cwd)
		# 	self.output().remove()
		# 	sys.exit()

class add_read_groups(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	picard_location = luigi.Parameter()
	library_prep = luigi.Parameter()
	platform = luigi.Parameter()

	# input_bam = luigi.Parameter()
	# output_bam = luigi.Parameter()
	sample = luigi.Parameter()

	def requires(self):
		return bowtie(fastq_file=self.fastq_file, sample=self.sample, max_threads=self.max_threads, project_dir=self.project_dir)

	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_preprocessed.bam'))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		cmd = ['java', '-jar', self.picard_location, 'AddOrReplaceReadGroups', 'I=%s' % self.input().path, 'O=%s' % self.output().path, 'SORT_ORDER=coordinate', 'RGID=%s' % self.sample, 'RGLB=%s' % self.library_prep, 'RGPL=%s' % self.platform, 'RGPU=%s' % self.sample + '_barcode', 'RGSM=%s' % self.sample]
		pipeline_utils.command_call(cmd, [self.output()], sleep_time=0.3)
		self.input().remove()

# https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
class mark_duplicates(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	picard_location = luigi.Parameter()

	# input_bam = luigi.Parameter() #replace all of the input bams with the outputs of the previous
	# output_bam = luigi.Parameter()
	# metrics_file = luigi.Parameter()
	sample = luigi.Parameter()

	def requires(self):
		return add_read_groups(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads)

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_marked_dups.bam')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'tmp', self.sample + '_marked_dups_metrics.txt'))]

	def run(self):
		pipeline_utils.confirm_path(self.output()[0].path)
		pipeline_utils.confirm_path(self.output()[1].path)
		cmd = ['java', '-jar', self.picard_location, 'MarkDuplicatesWithMateCigar', 'I=%s' % self.input().path, 'O=%s' % self.output()[0].path, 'M=%s' % self.output()[1].path]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.4)
		# self.input().remove()

class index_bam(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	samtools_location = luigi.Parameter()

	sample = luigi.Parameter()

	# input_bam = luigi.Parameter()
	# recalibration_report = luigi.Parameter()
	# output_bam = luigi.Parameter()

	def requires(self):
		return mark_duplicates(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads)

	def output(self):

		return [self.input()[0], luigi.LocalTarget(self.input()[0].path + '.bai')]

	def run(self):
		pipeline_utils.confirm_path(self.output()[1].path)
		# cmd = [os.getcwd() + '/' + self.samtools_location, 'index', '-b', self.input()[0].path]
		cmd = [self.samtools_location, 'index', '-b', self.input()[0].path]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.5)
		self.input()[1].remove()

# ~20 mins w/2 cores
class realigner_target(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	gatk3_location = luigi.Parameter()

	sample = luigi.Parameter()
	fasta_file = luigi.Parameter()
	mills = luigi.Parameter()
	kg = luigi.Parameter()

	# input_bam = luigi.Parameter()
	# output_file = luigi.Parameter()

	def requires(self):
		return index_bam(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads)

	def output(self):
		return self.input() + [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_forIndelRealigner.intervals'))]

	def run(self):
		pipeline_utils.confirm_path(self.output()[2].path)
		cmd = ['java', '-Xmx8g', '-XX:+UseSerialGC', '-jar', self.gatk3_location, '-T', 'RealignerTargetCreator', '-nt', str(self.max_threads), '-R', self.fasta_file, '-I', self.input()[0].path, '--known', self.mills, '--known', self.kg, '-o', self.output()[2].path]
		pipeline_utils.command_call(cmd, self.output(), threads_needed=self.max_threads, sleep_time=0.6)
		# for input_file in self.input():
		# 	input_file.remove()

# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php
# ~45 mins (no multiprocessing available)
class indel_realignment(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	gatk3_location = luigi.Parameter()

	fasta_file = luigi.Parameter()
	mills = luigi.Parameter()
	kg = luigi.Parameter()
	sample = luigi.Parameter()
	
	# input_bam = luigi.Parameter()
	# output_bam = luigi.Parameter()

	def requires(self):
		return realigner_target(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads)

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_realigned.bam')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_realigned.bai'))]

	def run(self):
		pipeline_utils.confirm_path(self.output()[0].path)
		pipeline_utils.confirm_path(self.output()[1].path)
		cmd = ['java', '-Xmx8g', '-XX:+UseSerialGC', '-jar', self.gatk3_location, '-T', 'IndelRealigner', '-R', self.fasta_file, '-I', self.input()[0].path, '-known', self.mills, '-known', self.kg, '-targetIntervals', self.input()[2].path, '-o', self.output()[0].path]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.7)
		for input_file in self.input():
			input_file.remove()

# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
# ~45 mins (no multiprocessing available)
class bqsr(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	gatk4_location = luigi.Parameter()

	fasta_file = luigi.Parameter()
	mills = luigi.Parameter()
	kg = luigi.Parameter()
	dbsnp = luigi.Parameter()
	sample = luigi.Parameter()
	
	# input_bam = luigi.Parameter()
	# output_table = luigi.Parameter()

	def requires(self):
		return indel_realignment(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads)

	def output(self):
		return self.input() + [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_recalibrated.table'))]

	def run(self):
		pipeline_utils.confirm_path(self.output()[2].path)
		cmd = [self.gatk4_location, '--java-options', '"-Xmx8g -XX:+UseSerialGC"', 'BaseRecalibrator', '-R', self.fasta_file, '-I', self.input()[0].path, '--known-sites', self.dbsnp, '--known-sites', self.mills, '--known-sites', self.kg, '-O',  self.output()[2].path]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.8)
		# self.input().remove()

# ~75 mins (no multiprocessing available)
class recalibrated_bam(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	gatk4_location = luigi.Parameter()

	fasta_file = luigi.Parameter()
	sample = luigi.Parameter()

	# input_bam = luigi.Parameter()
	# recalibration_report = luigi.Parameter()
	# output_bam = luigi.Parameter()

	def requires(self):
		return bqsr(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads)

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_recalibrated.bam')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], 'alignment', self.sample + '_recalibrated.bam'))]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		cmd = [self.gatk4_location, 'ApplyBQSR', '-R', self.fasta_file, '-I', self.input()[0].path, '--bqsr-recal-file', self.input()[2].path, '-O',  self.output()[0].path]
		pipeline_utils.command_call(cmd, self.output(), sleep_time=0.9)
		for input_file in self.input():
			input_file.remove()

class run_variant_caller(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	picard_location = luigi.Parameter()
	gatk3_location = luigi.Parameter()
	gatk4_location = luigi.Parameter()
	vardict_location = luigi.Parameter()
	library_bed = luigi.Parameter()

	caller = luigi.Parameter()
	case = luigi.Parameter()
	case_dir = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	fasta_file = luigi.Parameter()
	fasta_dir = os.path.join(*luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])
	
	def requires(self):
		# tumor_sam = os.path.join(self.fasta_dir, 'index', self.case + 'T' + '.sam')
		# normal_sam = os.path.join(self.fasta_dir, 'index', self.case + 'N' + '.sam')

		# tumor_bam = os.path.join(self.case_dir, 'alignment', self.case + 'T' + '.bam')
		# normal_bam = os.path.join(self.case_dir, 'alignment', self.case + 'N' + '.bam')
		if self.matched_n != '':
			# tumor_sam = os.path.join(self.fasta_dir, 'index', self.case + 'T' + '.sam')
			# normal_sam = os.path.join(self.fasta_dir, 'index', self.case + 'N' + '.sam')

			# tumor_bam = os.path.join(self.case_dir, 'alignment', self.case + 'T' + '.bam')
			# normal_bam = os.path.join(self.case_dir, 'alignment', self.case + 'N' + '.bam')

			# yield convert_bam(case_dir=self.case_dir, sample=self.case + 'T', fastq_file=self.tumor, sam_file=tumor_sam, bam_file=tumor_bam)
			# yield convert_bam(case_dir=self.case_dir, sample=self.case + 'N', fastq_file=self.matched_n, sam_file=normal_sam, bam_file=normal_bam)
			return [recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=int(self.max_threads/2)), recalibrated_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=int(self.max_threads/2))]
		else:
			# tumor_sam = os.path.join(self.fasta_dir, 'index', self.case + 'T' + '.sam')
			# tumor_bam = os.path.join(self.case_dir, 'alignment', self.case + 'T' + '.bam')
			# yield convert_bam(case_dir=self.case_dir, sample=self.case + 'T', fastq_file=self.tumor, sam_file=tumor_sam, bam_file=tumor_bam)
			return [recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads)]

	def output(self):
		vcf_file = os.path.join(self.vcf_path, self.caller, '.vcf')
		return luigi.LocalTarget(vcf_file)

	def run(self):
		tumor_bam = self.input()[0][0].path #os.path.join(self.case_dir, 'alignment', self.case + 'T' + '.bam')
		if self.matched_n:
			normal_bam = self.input()[1][0].path #os.path.join(self.case_dir, 'alignment', self.case + 'T' + '.bam')
		else:
			normal_bam = ''

		vcf_file = os.path.join(self.vcf_path, self.caller, '.vcf')
		
		caller_dict = {
		'MuTect': {'tumor_only': 'java -jar %s -T Mutect2 -R %s -I:tumor %s --dbsnp %s -o %s' % (self.gatk4_location, self.fasta_file, tumor_bam, self.known_vcf, vcf_file),
			'matched_n': 'java -jar %s -T Mutect2 -R %s -I:tumor %s -I:normal %s --dbsnp %s -o %s' % (self.gatk4_location, self.fasta_file, tumor_bam, normal_bam, self.known_vcf, vcf_file)},
		'VarDict': {'tumor_only': 'AF_THR="0.01" \
vardict -G %s -f $AF_THR -N %s -b %s -c 1 -S 2 -E 3 -g 4 -x 100 %s | teststrandbias.R | var2vcf_valid.pl -N %s -E -f $AF_THR' % (self.fasta_file, self.case + 'T', tumor_bam, self.library_bed, self.case + 'T'),
			'matched_n': 'AF_THR="0.01" \
vardict -G %s -f $AF_THR -N %s -b "%s|%s" -c 1 -S 2 -E 3 -g 4 -x 100 %s | testsomatic.R | var2vcf_paired.pl -N "%s|%s" -f $AF_THR' % (self.fasta_file, self.case + 'T', tumor_bam, normal_bam, self.library_bed, self.case + 'T', self.case + 'N'),},
		'FreeBayes': {'tumor_only': '',
			'matched_n': ''},
		'VarScan': {'tumor_only': '',
			'matched_n': ''},
		'Scalpel': {'tumor_only': '',
			'matched_n': ''},
		'LUMPY': {'tumor_only': '',
			'matched_n': ''},
		'DELLY': {'tumor_only': '',
			'matched_n': ''},
		'WHAM': {'tumor_only': '',
			'matched_n': ''},
		'CNVkit': {'tumor_only': '',
			'matched_n': ''}
		}
		if self.matched_n:
			subprocess.call([caller_dict[self.caller]['matched_n']], stdout=subprocess.PIPE)
		else:
			subprocess.call([caller_dict[self.caller]['tumor_only']], stdout=subprocess.PIPE)

class aggregate_variants(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()


	
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
		'max_threads': self.max_threads
		}
		yield filter_mutect(**kwargs)
		yield scalpel_export(**kwargs)
		yield vardict(**kwargs)
		yield cnvkit(case_dict=self.case_dict, **kwargs)
		# yield run_variant_caller(caller='VarDict', **kwargs)
		# yield run_variant_caller(caller='FreeBayes', **kwargs)
		# yield run_variant_caller(caller='VarScan', **kwargs)

		# yield run_variant_caller(caller='Scalpel', **kwargs)

		# yield run_variant_caller(caller='CNVkit', **kwargs)

		# yield run_variant_caller(caller='LUMPY', **kwargs)
		# yield run_variant_caller(caller='DELLY', **kwargs)
		# yield run_variant_caller(caller='WHAM', **kwargs)

	def output(self):
		return self.input()

class cases(luigi.Task):
	sample_dict = luigi.DictParameter()
	project_dir = luigi.Parameter()
	sample_threads = luigi.IntParameter()

	def requires(self):
		return [aggregate_variants(case=case, tumor=self.sample_dict[case]['T'], matched_n=self.sample_dict[case]['N'], project_dir=self.project_dir, max_threads=self.sample_threads, case_dict=self.sample_dict) for case in self.sample_dict]
		
		# # global global_max_threads, thread_count
		
		# global_vars.global_max_threads = self.max_threads
		# global_vars.thread_file = os.path.join(os.getcwd(), 'thread_count_temp_%s.txt' % timestamp)
		# print(global_vars.thread_file)
		# pipeline_utils.init_thread_file(global_vars.thread_file)
		# global_vars.working_files = os.path.join(os.getcwd(), 'working_files_%s.pkl' % timestamp)
		# pipeline_utils.init_working_files(global_vars.working_files)
		# global_vars.cwd = os.getcwd()

		# sample_dict = {}
		# # try:
		# for sample in os.listdir(self.sample_dir):
		# 	if os.path.isdir(os.path.join(self.sample_dir, sample)):
		# 		sample_dict[sample] = {'T': '', 'N': ''}
		# 		tumor_fastq = os.path.join(self.sample_dir, sample, 'tumor', os.listdir(os.path.join(self.sample_dir, sample, 'tumor'))[0]) + '\t' + os.path.join(self.sample_dir, sample, 'tumor', os.listdir(os.path.join(self.sample_dir, sample, 'tumor'))[1])
		# 		sample_dict[sample]['T'] = tumor_fastq
		# 		if os.path.exists(os.path.join(self.sample_dir, sample, 'normal')):
		# 			if len(os.listdir(os.path.join(self.sample_dir, sample, 'normal'))) > 0:
		# 				normal_fastq = os.path.join(self.sample_dir, sample, 'normal', os.listdir(os.path.join(self.sample_dir, sample, 'normal'))[0]) + '\t' + os.path.join(self.sample_dir, sample, 'normal', os.listdir(os.path.join(self.sample_dir, sample, 'normal'))[1])
		# 				sample_dict[sample]['N'] = normal_fastq
		# # except:
		# # 	raise ValueError("Error in parsing fastq directory.")
		# # print('\n\n\n\n')
		# # print(self.sample_dir)
		# # print(sample_dict)
		# if self.threads_per_sample:
		# 	sample_threads = self.threads_per_sample
		# else:
		# 	sample_threads = max(1, self.max_threads//len(sample_dict.keys()))
		# for case in sample_dict:
		# 	tumor = sample_dict[case]['T']
		# 	matched_n = sample_dict[case]['N']
		# 	yield aggregate_variants(case=case, tumor=tumor, matched_n=matched_n, project_dir=self.project_dir, max_threads=sample_threads, case_dict=sample_dict)

def run_pipeline(args):
	timestamp = str(int(time.time()))

	global_vars.global_max_threads = args.max_threads
	global_vars.thread_file = os.path.join(os.getcwd(), 'thread_count_temp_%s.txt' % timestamp)
	print(global_vars.thread_file)
	pipeline_utils.init_thread_file(global_vars.thread_file)
	global_vars.working_files = os.path.join(os.getcwd(), 'working_files_%s.pkl' % timestamp)
	pipeline_utils.init_working_files(global_vars.working_files)
	global_vars.cwd = os.getcwd()

	sample_dict = {}
	# try:
	for sample in os.listdir(args.sample_dir):
		if os.path.isdir(os.path.join(args.sample_dir, sample)):
			sample_dict[sample] = {'T': '', 'N': ''}
			tumor_list = [filename for filename in os.listdir(os.path.join(args.sample_dir, sample, 'tumor')) if 'fastq' in filename]
			tumor_fastq = os.path.join(args.sample_dir, sample, 'tumor', tumor_list[0]) + '\t' + os.path.join(args.sample_dir, sample, 'tumor', tumor_list[1])
			sample_dict[sample]['T'] = tumor_fastq
			if os.path.exists(os.path.join(args.sample_dir, sample, 'normal')):
				if len(os.listdir(os.path.join(args.sample_dir, sample, 'normal'))) > 0:
					normal_list = [filename for filename in os.listdir(os.path.join(args.sample_dir, sample, 'normal')) if 'fastq' in filename]
					normal_fastq = os.path.join(args.sample_dir, sample, 'normal', normal_list[0]) + '\t' + os.path.join(args.sample_dir, sample, 'normal', normal_list[1])
					sample_dict[sample]['N'] = normal_fastq
	# except:
	# 	raise ValueError("Error in parsing fastq directory.")
	# print('\n\n\n\n')
	# print(self.sample_dir)
	# print(sample_dict)
	if args.threads_per_sample:
		sample_threads = args.threads_per_sample
	else:
		sample_threads = max(1, args.max_threads//len(sample_dict.keys()))

	worker_scheduler_factory = luigi.interface._WorkerSchedulerFactory()
	# for case in sample_dict:
	# 	tumor = sample_dict[case]['T']
	# 	matched_n = sample_dict[case]['N']
	# luigi.build([bam_processing.aggregate_variants(case=case, tumor=sample_dict[case]['T'], matched_n=sample_dict[case]['N'], project_dir=args.project_dir, max_threads=sample_threads, case_dict=sample_dict) for case in sample_dict], workers=args.workers, local_scheduler=args.local_scheduler, worker_scheduler_factory=worker_scheduler_factory) #, scheduler_port=int(args.port)) # workers=sample_threads
	luigi.build([cases(sample_dict=sample_dict, project_dir=args.project_dir, sample_threads=sample_threads)], workers=args.workers, local_scheduler=args.local_scheduler, worker_scheduler_factory=worker_scheduler_factory) #, scheduler_port=int(args.port)) # workers=sample_threads

		# [(max_threads=args.max_threads, project_dir=args.project_dir, sample_dir=args.sample_dir, threads_per_sample=args.threads_per_sample, timestamp=timestamp)], workers=args.workers, local_scheduler=args.local_scheduler)
		# yield aggregate_variants(case=case, tumor=tumor, matched_n=matched_n, project_dir=self.project_dir, max_threads=sample_threads, case_dict=sample_dict)



	
	
