#!/usr/bin/env python3
import subprocess
import luigi
import luigi.interface
import os
import sys
import time
import pipeline_utils
import bam_processing

class mutect(luigi.Task):
	max_threads = luigi.IntParameter()
	matched_n = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()
	project_dir = luigi.Parameter()

	gatk4_location = luigi.Parameter()
	gatk3_location = luigi.Parameter()
	known_vcf = luigi.Parameter()
	germline_resource = luigi.Parameter()
	fasta_file = luigi.Parameter()
	fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=int(self.max_threads/2)), bam_processing.recalibrated_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=int(self.max_threads/2))]
		else:
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads)]


	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_mutect' + '.vcf.gz'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		if self.matched_n:
			cmd = [self.gatk4_location, 'Mutect2', '-R', self.fasta_file, '-I', self.input()[0][0].path, '-tumor', self.case + '_T', '-I', self.input()[1][0].path, '-normal', self.case + '_N', '--germline-resource', self.germline_resource, '-O', self.output().path]
		else:
			cmd = [self.gatk4_location, 'Mutect2', '-R', self.fasta_file, '-I', self.input()[0][0].path, '-tumor', self.case + '_T', '--germline-resource', self.germline_resource, '-O', self.output().path]
		pipeline_utils.command_call(cmd, [self.output()], threads_needed=4)

class filter_mutect(luigi.Task):
	max_threads = luigi.IntParameter()
	matched_n = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()
	project_dir = luigi.Parameter()

	gatk4_location = luigi.Parameter()
	fasta_file = luigi.Parameter()
	fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	def requires(self):
		if self.matched_n != '':
			return mutect(project_dir=self.project_dir, vcf_path=self.vcf_path, case=self.case, tumor=self.tumor, matched_n=self.matched_n, max_threads=self.max_threads)
		else:
			return mutect(project_dir=self.project_dir, vcf_path=self.vcf_path, case=self.case, tumor=self.tumor, matched_n=self.matched_n, max_threads=self.max_threads)


	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_mutect_filtered' + '.vcf.gz'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		cmd = [self.gatk4_location, 'FilterMutectCalls', '-V', self.input().path, '-O', self.output().path]
		pipeline_utils.command_call(cmd, [self.output()], sleep_time=1.1)
		self.input().remove()


class scalpel_discovery(luigi.Task):
	max_threads = luigi.IntParameter()
	matched_n = luigi.Parameter()
	project_dir = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	library_bed = luigi.Parameter()
	fasta_file = luigi.Parameter()

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=int(self.max_threads/2)), bam_processing.recalibrated_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=int(self.max_threads/2))]
		else:
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads)]


	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, 'scalpel', 'variants.db.dir'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		if self.matched_n:
			cmd = ['./packages/scalpel-0.5.4/scalpel-discovery', '--somatic', '--normal', self.input()[1][0].path, '--tumor', self.input()[0][0].path, '--bed', self.library_bed, '--ref', self.fasta_file, '--two-pass', '--dir', os.path.join(self.vcf_path, 'scalpel'), '--numprocs', str(self.max_threads)]
		else:
			cmd = ['./packages/scalpel-0.5.4/scalpel-discovery', '--single', '--bam', self.input()[0][0].path, '--bed', self.library_bed, '--ref', self.fasta_file, '--dir', os.path.join(self.vcf_path, 'scalpel'), '--numprocs', str(self.max_threads)]
		pipeline_utils.command_call(cmd, [self.output()], threads_needed=self.max_threads)

class scalpel_export(luigi.Task):
	max_threads = luigi.IntParameter()
	matched_n = luigi.Parameter()
	project_dir = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	library_bed = luigi.Parameter()
	fasta_file = luigi.Parameter()
	fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	def requires(self):
		# if self.matched_n != '':
		return scalpel_discovery(case=self.case, tumor=self.tumor, matched_n=self.matched_n, project_dir=self.project_dir, vcf_path=self.vcf_path, max_threads=self.max_threads) #, scalpel_discovery(case=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=int(self.max_threads/2))]
		# else:
		# 	return scalpel_discovery(case=self.case + '_T', tumor=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads)


	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_scalpel' + '.vcf'))
	
	def run(self):
		if self.matched_n:
			cmd = ['./scalpel-0.5.4/scalpel-export', '--somatic', '--db', self.input().path[:-4], '--bed', self.library_bed, '--ref', self.fasta_file]
		else:
			cmd = ['./scalpel-0.5.4/scalpel-export', '--single', '--db', self.input().path[:-4], '--bed', self.library_bed, '--ref', self.fasta_file]
		pipeline_utils.command_call(cmd, [self.output()], sleep_time=1.1)

class freebayes(luigi.Task):
	max_threads = luigi.IntParameter()
	matched_n = luigi.Parameter()
	project_dir = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	library_bed = luigi.Parameter()
	fasta_file = luigi.Parameter()

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=int(self.max_threads/2)), bam_processing.recalibrated_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=int(self.max_threads/2))]
		else:
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads)]


	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_freebayes' + '.vcf'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		if self.matched_n:
			cmd = ['./packages/freebayes/scripts/freebayes-parallel', '<(fasta_generate_regions.py %s.fai 100000)' % self.fasta_file, self.max_threads, '-f', self.fasta_file, '-t', self.library_bed, '--pooled-continuous', '--pooled-discrete', '-F', '0.01', '-C', '2', self.input()[0][0].path, self.input()[1][0].path, '>', os.path.join(self.vcf_path, 'freebayes')]
		else:
			cmd = ['./packages/freebayes/scripts/freebayes-parallel', '<(fasta_generate_regions.py %s.fai 100000)' % self.fasta_file, self.max_threads, '-f', self.fasta_file, '-t', self.library_bed, '--pooled-continuous', '--pooled-discrete', '-F', '0.01', '-C', '2', self.input()[0][0].path, '>', os.path.join(self.vcf_path, 'freebayes')]
		pipeline_utils.command_call(cmd, [self.output()])

class vardict(luigi.Task):
	max_threads = luigi.IntParameter()
	matched_n = luigi.Parameter()
	project_dir = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	library_bed = luigi.Parameter()
	fasta_file = luigi.Parameter()

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=int(self.max_threads/2)), bam_processing.recalibrated_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=int(self.max_threads/2))]
		else:
			return [bam_processing.recalibrated_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads)]


	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_freebayes' + '.vcf'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		if self.matched_n:
			cmd = ['./packagers//VarDictJava/build/install/VarDict/bin/VarDict', '-G', self.fasta_file, '-f', '0.01', '-N', self.case + '_T', '-b', '"%s|%s"' % (self.input()[0][0].path, self.input()[1][0].path), '-z', '-F', '-c', '1', '-S', '2', '-E', '3', '-g', '4', self.library_bed, '|', './packages/VarDictJava/VarDict/testsomatic.R', '|', './packages/VarDictJava/VarDict/var2vcf_paired.pl', '-N', '"%s|%s"' % (self.case + '_T', self.case + '_N'), '-f', '0.01', '>%s' % os.path.join(self.vcf_path, 'vardict')]
		else:
			cmd = ['./packages//VarDictJava/build/install/VarDict/bin/VarDict', '-G', self.fasta_file, '-f', '0.01', '-N', self.case + '_T', '-b', self.input()[0][0].path, '-z', '-c', '1', '-S', '2', '-E', '3', '-g', '4', self.library_bed, '|', './packages/VarDictJava/VarDict/teststrandbias.R', '|', './packages/VarDictJava/VarDict/var2vcf_valid.pl', '-N', self.case + '_T', 'E', '-f', '0.01', '>%s' % os.path.join(self.vcf_path, 'vardict')]
		pipeline_utils.command_call(cmd, [self.output()])
