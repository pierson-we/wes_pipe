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
			cmd = ['java', '-jar', self.gatk4_location, '-T', 'Mutect2', '-R', self.fasta_file, '-I:tumor', self.input()[0].path, '-I:normal', self.input()[1].path, '--dbsnp', self.known_vcf, '-o', self.output().path]
		else:
			cmd = [self.gatk4_location, 'Mutect2', '-R', self.fasta_file, '-I', self.input()[0].path, '-tumor', self.case + '_T', '--germline-resource', self.germline_resource, '-O', self.output().path]
		pipeline_utils.command_call(cmd, threads_needed=4)

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
		pipeline_utils.command_call(cmd)


class scalpel_discovery(luigi.Task):
	max_threads = luigi.Parameter()
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
			cmd = ['./scalpel-0.5.4/scalpel-discovery', '--somatic', '--normal', self.input()[1].path, '--tumor', self.input()[0].path, '--bed', self.library_bed, '--ref', self.fasta_file, '--two-pass', '--dir', os.path.join(self.vcf_path, 'scalpel'), '--numprocs', str(self.max_threads)]
		else:
			cmd = ['./scalpel-0.5.4/scalpel-discovery', '--single', '--bam', self.input()[0].path, '--bed', self.library_bed, '--ref', self.fasta_file, '--dir', os.path.join(self.vcf_path, 'scalpel'), '--numprocs', str(self.max_threads)]
		pipeline_utils.command_call(cmd)

class scalpel_export(luigi.Task):
	max_threads = luigi.Parameter()
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
		pipeline_utils.command_call(cmd)

class freebayes(luigi.Task):
	max_threads = luigi.Parameter()
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
			cmd = ['./scalpel-0.5.4/scalpel-discovery', '--somatic', '--normal', self.input()[1].path, '--tumor', self.input()[0].path, '--bed', self.library_bed, '--ref', self.fasta_file, '--two-pass', '--dir', os.path.join(self.vcf_path, 'scalpel'), '--numprocs', str(self.max_threads)]
		else:
			cmd = ['./scalpel-0.5.4/scalpel-discovery', '--single', '--bam', self.input()[0].path, '--bed', self.library_bed, '--ref', self.fasta_file, '--dir', os.path.join(self.vcf_path, 'scalpel'), '--numprocs', str(self.max_threads)]
		pipeline_utils.command_call(cmd)
# 'VarDict': {'tumor_only': 'AF_THR="0.01" \
# vardict -G %s -f $AF_THR -N %s -b %s -c 1 -S 2 -E 3 -g 4 -x 100 %s | teststrandbias.R | var2vcf_valid.pl -N %s -E -f $AF_THR' % (self.fasta_file, self.case + 'T', tumor_bam, self.library_bed, self.case + 'T'),
# 			'matched_n': 'AF_THR="0.01" \
# vardict -G %s -f $AF_THR -N %s -b "%s|%s" -c 1 -S 2 -E 3 -g 4 -x 100 %s | testsomatic.R | var2vcf_paired.pl -N "%s|%s" -f $AF_THR' % (self.fasta_file, self.case + 'T', tumor_bam, normal_bam, self.library_bed, self.case + 'T', self.case + 'N'),},
# 		