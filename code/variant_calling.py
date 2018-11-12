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
import bam_processing
import misc_utils

class mutect_single_normal(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	sample = luigi.Parameter()
	fastq_file = luigi.Parameter()
	# tumor = luigi.Parameter()
	# matched_n = luigi.Parameter()
	# vcf_path = luigi.Parameter()
	# case_dict = luigi.DictParameter()

	cfg = luigi.DictParameter()

	def requires(self):
		# if self.matched_n != '':
		# 	return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads), bam_processing.index_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads)]
		# else:
		return bam_processing.index_bam(sample=self.sample, fastq_file=self.fastq_file, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)


	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'mutect', self.sample + '.vcf.gz')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'mutect', self.sample + '.vcf.gz.tbi'))]
		
	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		# if self.matched_n:
		# 	cmd = ['./packages/VarDictJava/build/install/VarDict/bin/VarDict', '-G', self.cfg['fasta_file'], '-f', '0.01', '-N', self.case + '_T', '-b', '"%s|%s"' % (self.input()[0][0].path, self.input()[1][0].path), '-z', '-F', '-c', '1', '-S', '2', '-E', '3', '-g', '4', self.cfg['library_bed'], '|', './packages/VarDictJava/VarDict/testsomatic.R', '|', './packages/VarDictJava/VarDict/var2vcf_paired.pl', '-N', '"%s|%s"' % (self.case + '_T', self.case + '_N'), '-f', '0.01', '>%s' % os.path.join(self.vcf_path, 'vardict')]
		# else:
		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx8g -Xms8g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'Mutect2', '-R', self.cfg['fasta_file'], '-I', self.input()[0].path, '-tumor', self.sample, '-L', self.cfg['library_bed'], '--native-pair-hmm-threads', self.max_threads, '-O', self.output()[0].path]
		pipeline_utils.command_call(cmd, self.output(), threads_needed=self.max_threads)

class mutect_pon(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# case = luigi.Parameter()
	# tumor = luigi.Parameter()
	# matched_n = luigi.Parameter()
	# vcf_path = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		# if self.matched_n != '':
		# 	return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads), bam_processing.index_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads)]
		# else:
		return [mutect_single_normal(sample=case_name + '_N', fastq_file=self.case_dict[case_name]['N'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg) for case_name in self.case_dict if self.case_dict[case_name]['N'] != '']


	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'mutect', 'pon.vcf.gz'))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		# if self.matched_n:
		# 	cmd = ['./packages/VarDictJava/build/install/VarDict/bin/VarDict', '-G', self.cfg['fasta_file'], '-f', '0.01', '-N', self.case + '_T', '-b', '"%s|%s"' % (self.input()[0][0].path, self.input()[1][0].path), '-z', '-F', '-c', '1', '-S', '2', '-E', '3', '-g', '4', self.cfg['library_bed'], '|', './packages/VarDictJava/VarDict/testsomatic.R', '|', './packages/VarDictJava/VarDict/var2vcf_paired.pl', '-N', '"%s|%s"' % (self.case + '_T', self.case + '_N'), '-f', '0.01', '>%s' % os.path.join(self.vcf_path, 'vardict')]
		# else:
		cmd = [self.cfg['gatk4_location'], 'CreateSomaticPanelOfNormals']
		for normal_vcf in self.input():
			cmd.append('--vcfs')
			cmd.append(normal_vcf[0].path)
		cmd.append('--output')
		cmd.append(self.output().path)
		pipeline_utils.command_call(cmd, [self.output()])

class mutect(luigi.Task):
	max_threads = luigi.IntParameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()
	project_dir = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# gatk4_location = luigi.Parameter()
	# gatk3_location = luigi.Parameter()
	# dbsnp = luigi.Parameter()
	# germline_resource = luigi.Parameter()
	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()
	# fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	cfg = luigi.DictParameter()

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), bam_processing.index_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), mutect_pon(case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]
		else:
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), mutect_pon(case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]


	def output(self):
		return [luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_mutect.vcf.gz')), luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_mutect.vcf.gz.tbi'))]
	
	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		if self.matched_n:
			cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx8g -Xms8g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'Mutect2', '-R', self.cfg['fasta_file'], '-I', self.input()[0][0].path, '-tumor', self.case + '_T', '-I', self.input()[1][0].path, '-normal', self.case + '_N', '--germline-resource', self.cfg['germline_resource'], '--af-of-alleles-not-in-resource', '0.0000025', '-L', self.cfg['library_bed'], '-pon', self.input()[-1].path, '--native-pair-hmm-threads', self.max_threads, '-O', self.output()[0].path]
		else:
			cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx8g -Xms8g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'Mutect2', '-R', self.cfg['fasta_file'], '-I', self.input()[0][0].path, '-tumor', self.case + '_T', '--germline-resource', self.cfg['germline_resource'], '--af-of-alleles-not-in-resource', '0.0000025', '-L', self.cfg['library_bed'], '-pon', self.input()[-1].path, '--native-pair-hmm-threads', self.max_threads, '-O', self.output()[0].path]
		pipeline_utils.command_call(cmd, self.output(), threads_needed=self.max_threads)

class filter_mutect(luigi.Task):
	max_threads = luigi.IntParameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()
	project_dir = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# gatk4_location = luigi.Parameter()
	# fasta_file = luigi.Parameter()
	# fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	cfg = luigi.DictParameter()

	def requires(self):
		return mutect(project_dir=self.project_dir, vcf_path=self.vcf_path, case=self.case, tumor=self.tumor, matched_n=self.matched_n, max_threads=self.max_threads, case_dict=self.case_dict, cfg=self.cfg)

	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_mutect_filtered' + '.vcf.gz'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx8g -Xms8g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'FilterMutectCalls', '-V', self.input()[0].path, '-O', self.output().path]
		pipeline_utils.command_call(cmd, [self.output()], sleep_time=1.1)
		# for input_file in self.input():
		# 	input_file.remove()


class scalpel_discovery(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), bam_processing.index_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]
		else:
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]


	def output(self):
		if self.matched_n != '':
			return luigi.LocalTarget(os.path.join(self.vcf_path, 'scalpel', 'twopass', 'somatic.db.dir'))
		else:
			return luigi.LocalTarget(os.path.join(self.vcf_path, 'scalpel', 'variants.db.dir'))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		if self.matched_n:
			cmd = ['./packages/scalpel-0.5.4/scalpel-discovery', '--somatic', '--normal', self.input()[1][0].path, '--tumor', self.input()[0][0].path, '--bed', self.cfg['library_bed'], '--ref', self.cfg['fasta_file'], '--two-pass', '--dir', os.path.join(self.vcf_path, 'scalpel'), '--numprocs', str(self.max_threads)]
		else:
			cmd = ['./packages/scalpel-0.5.4/scalpel-discovery', '--single', '--bam', self.input()[0][0].path, '--bed', self.cfg['library_bed'], '--ref', self.cfg['fasta_file'], '--dir', os.path.join(self.vcf_path, 'scalpel'), '--numprocs', str(self.max_threads)]
		pipeline_utils.command_call(cmd, [self.output()], threads_needed=self.max_threads)

class scalpel_export(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()
	# fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	cfg = luigi.DictParameter()

	def requires(self):
		# if self.matched_n != '':
		return scalpel_discovery(case=self.case, tumor=self.tumor, matched_n=self.matched_n, project_dir=self.project_dir, vcf_path=self.vcf_path, max_threads=self.max_threads, cfg=self.cfg) #, scalpel_discovery(case=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=int(self.max_threads/2))]
		# else:
		# 	return scalpel_discovery(case=self.case + '_T', tumor=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads)


	def output(self):
		if self.matched_n != '':
			return luigi.LocalTarget(os.path.join(self.vcf_path, 'scalpel', 'twopass', 'somatic.indel.vcf'))
		else:
			return luigi.LocalTarget(os.path.join(self.vcf_path, 'scalpel', 'variants.indel.vcf'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		if self.matched_n:
			cmd = ['./packages/scalpel-0.5.4/scalpel-export', '--somatic', '--db', self.input().path[:-4], '--bed', self.cfg['library_bed'], '--ref', self.cfg['fasta_file']]
		else:
			cmd = ['./packages/scalpel-0.5.4/scalpel-export', '--single', '--db', self.input().path[:-4], '--bed', self.cfg['library_bed'], '--ref', self.cfg['fasta_file']]
		pipeline_utils.command_call(cmd, [self.output()], sleep_time=1.1)

# not yet tested - need to install GNU Parallel on cluster... but might be able to run local install http://git.savannah.gnu.org/cgit/parallel.git/tree/README
class freebayes(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), bam_processing.index_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]
		else:
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]


	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_freebayes' + '.vcf'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		if self.matched_n:
			cmd = ['./packages/freebayes/scripts/freebayes-parallel', '''<(/bin/bash python3 ./packages/freebayes/scripts/fasta_generate_regions.py''', '%s.fai' % self.cfg['fasta_file'], '''100000)''', self.max_threads, '-f', self.cfg['fasta_file'], '-t', self.cfg['library_bed'], '--pooled-continuous', '--pooled-discrete', '-F', '0.01', '-C', '2', self.input()[0][0].path, self.input()[1][0].path, '>', os.path.join(self.vcf_path, 'freebayes.vcf')]
		else:
			cmd = ['./packages/freebayes/scripts/freebayes-parallel', '''<(/bin/bash python3 ./packages/freebayes/scripts/fasta_generate_regions.py''', '%s.fai' % self.cfg['fasta_file'], '''100000)''', self.max_threads, '-f', self.cfg['fasta_file'], '-t', self.cfg['library_bed'], '--pooled-continuous', '--pooled-discrete', '-F', '0.01', '-C', '2', self.input()[0][0].path, '>', os.path.join(self.vcf_path, 'freebayes.vcf')]
		pipeline_utils.command_call(cmd, [self.output()])

class vardict(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), bam_processing.index_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]
		else:
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]


	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_vardict' + '.vcf'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		if self.matched_n:
			cmd = ['./packages/VarDictJava/build/install/VarDict/bin/VarDict', '-G', self.cfg['fasta_file'], '-f', '0.01', '-N', self.case + '_T', '-b', '"%s|%s"' % (self.input()[0][0].path, self.input()[1][0].path), '-th', self.max_threads, '-z', '-c', '1', '-S', '2', '-E', '3', '-g', '4', self.cfg['library_bed'], '|', './packages/VarDictJava/VarDict/testsomatic.R', '|', './packages/VarDictJava/VarDict/var2vcf_paired.pl', '-N', '"%s|%s"' % (self.case + '_T', self.case + '_N'), '-f', '0.01', '>%s' % self.output().path]
		else:
			cmd = ['./packages/VarDictJava/build/install/VarDict/bin/VarDict', '-G', self.cfg['fasta_file'], '-f', '0.01', '-N', self.case + '_T', '-b', self.input()[0][0].path, '-th', self.max_threads, '-z', '-c', '1', '-S', '2', '-E', '3', '-g', '4', self.cfg['library_bed'], '|', './packages/VarDictJava/VarDict/teststrandbias.R', '|', './packages/VarDictJava/VarDict/var2vcf_valid.pl', '-N', self.case + '_T', '-E', '-f', '0.01', '>%s' % self.output().path]
		pipeline_utils.command_call(cmd, [self.output()], threads_needed=self.max_threads)

class sort_vardict(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return [vardict(max_threads=self.max_threads, project_dir=self.project_dir, case=self.case, tumor=self.tumor, matched_n=self.matched_n, vcf_path=self.vcf_path, cfg=self.cfg),
		bam_processing.picard_index(cfg=self.cfg)]

	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_vardict_sorted' + '.vcf'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		cmd = ['java', '-jar', self.cfg['picard_location'], 'SortVcf', 'I=%s' % self.input()[0].path, 'O=%s' % self.output().path, 'SEQUENCE_DICTIONARY=%s' % self.input()[1].path]
		pipeline_utils.command_call(cmd, [self.output()], threads_needed=self.max_threads)

# this will be pretty annoying to get up and going
class varscan(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		if self.matched_n != '':
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), bam_processing.index_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]
		else:
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]


	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_varscan' + '.vcf'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		if self.matched_n:
			cmd = ['./packages/VarDictJava/build/install/VarDict/bin/VarDict', '-G', self.cfg['fasta_file'], '-f', '0.01', '-N', self.case + '_T', '-b', '"%s|%s"' % (self.input()[0][0].path, self.input()[1][0].path), '-z', '-F', '-c', '1', '-S', '2', '-E', '3', '-g', '4', self.cfg['library_bed'], '|', './packages/VarDictJava/VarDict/testsomatic.R', '|', './packages/VarDictJava/VarDict/var2vcf_paired.pl', '-N', '"%s|%s"' % (self.case + '_T', self.case + '_N'), '-f', '0.01', '>%s' % os.path.join(self.vcf_path, 'vardict')]
		else:
			cmd = ['./packages/VarDictJava/build/install/VarDict/bin/VarDict', '-G', self.cfg['fasta_file'], '-f', '0.01', '-N', self.case + '_T', '-b', self.input()[0][0].path, '-z', '-c', '1', '-S', '2', '-E', '3', '-g', '4', self.cfg['library_bed'], '|', './packages/VarDictJava/VarDict/teststrandbias.R', '|', './packages/VarDictJava/VarDict/var2vcf_valid.pl', '-N', self.case + '_T', 'E', '-f', '0.01', '>%s' % os.path.join(self.vcf_path, 'vardict')]
		pipeline_utils.command_call(cmd, [self.output()])

class pindel(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# case = luigi.Parameter()
	# tumor = luigi.Parameter()
	# matched_n = luigi.Parameter()
	# vcf_path = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return [bam_processing.index_bam(sample=case_name + '_N', fastq_file=self.case_dict[case_name]['N'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg) for case_name in self.case_dict if self.case_dict[case_name]['N'] != ''] \
		+ [bam_processing.index_bam(sample=case_name + '_T', fastq_file=self.case_dict[case_name]['T'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg) for case_name in self.case_dict]

	def output(self):
		pindel_files = ['_D', '_SI', '_TD', '_INV'] #, '_LI', '_BP',
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'pindel', 'pindel_all_samples' + ext)) for ext in pindel_files]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		with open('___pindel_bams___.txt', 'w') as f:
			for input_bam in self.input():
				case = input_bam[0].path.split('/')[-1].split('_')[0]
				if '_N' in input_bam[0].path:
					f.write('%s %s %s\n' % (input_bam[0].path, self.cfg['insert_size'], case + '_N'))
				else:
					f.write('%s %s %s\n' % (input_bam[0].path, self.cfg['insert_size'], case + '_T'))
		cmd = ['./packages/pindel/pindel', '-f', self.cfg['fasta_file'], '-i', '___pindel_bams___.txt', '-T', self.max_threads, '-c', 'ALL', '-o', os.path.join(self.project_dir, 'pindel', 'pindel_all_samples')]
		pipeline_utils.command_call(cmd, self.output(), threads_needed=self.max_threads)

		os.remove('___pindel_bams___.txt')

class filter_pindel(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# case = luigi.Parameter()
	# tumor = luigi.Parameter()
	# matched_n = luigi.Parameter()
	# vcf_path = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return pindel(case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', case_name, 'variants', case_name + '_T.pindel.bed')) for case_name in self.case_dict] + \
		[luigi.LocalTarget(os.path.join(self.project_dir, 'output', case_name, 'variants', case_name + '_N.pindel.bed')) for case_name in self.case_dict if self.case_dict[case_name]['N'] != ''] + \
		[luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'all_samples', 'all_samples_pindel.tsv'))]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)

		wait_time = random.uniform(0,3)
		time.sleep(wait_time)
		sys.stdout.flush()
		while not pipeline_utils.add_thread_count(global_vars.thread_file, 1):
			time.sleep(1.2)

		sample_dict = {output.path.split('/')[-1].split('.pindel.bed')[0]: output.path for output in self.output()[:-1]}

		misc_utils.filter_pindel(pindel_files=[input_file.path for input_file in self.input()], sample_dict=sample_dict, project_dir=self.project_dir, all_samples_output=self.output()[-1].path, min_reads=self.cfg['pindel_min_reads'], min_qual=self.cfg['pindel_min_qual'], max_inv_length=self.cfg['pindel_max_inv_length'])

		while not pipeline_utils.sub_thread_count(global_vars.thread_file, 1):
			time.sleep(1.2)

class annotate_pindel(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# case = luigi.Parameter()
	# tumor = luigi.Parameter()
	# matched_n = luigi.Parameter()
	# vcf_path = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return filter_pindel(case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)

	def output(self):
		outputs = []
		for input_file in self.input()[:-1]:
			output_file = input_file.path.split('.bed')[0] + '_final.bed'
			outputs.append(luigi.LocalTarget(output_file))
		return outputs

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)

		wait_time = random.uniform(0,3)
		time.sleep(wait_time)
		sys.stdout.flush()
		while not pipeline_utils.add_thread_count(global_vars.thread_file, 1):
			time.sleep(1.2)

		for i, input_file in enumerate(self.input()[:-1]):
			cmd = 'sort-bed %s' % input_file.path
			# print(cmd)
			p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

			cmd = ['bedtools', 'intersect', '-wa', '-u', '-sorted', '-a', 'stdin', '-b', self.cfg['exons_bed']]
			cmd = ' '.join(cmd)
			p2 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=p1.stdout, shell=True)
			print(cmd)

			cmd = ["bedmap", "--echo", "--echo-map-id-uniq", "--delim", r"'\t'", "-", self.cfg['genmap']]
			cmd = " ".join(cmd)
			# print(cmd)
			p3 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=p2.stdout, shell=True)

			outs, err = p3.communicate()
			with open(self.output()[i].path, 'wb') as f:
				f.write(str.encode('#gffTags\n'))
				f.write(outs)

		while not pipeline_utils.sub_thread_count(global_vars.thread_file, 1):
			time.sleep(1.2)

class filter_pindel_old(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# case = luigi.Parameter()
	# tumor = luigi.Parameter()
	# matched_n = luigi.Parameter()
	# vcf_path = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return pindel(case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)

	def output(self):
		return [luigi.LocalTarget(input_file.path + '.filtered.tsv') for input_file in self.input()]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)

		wait_time = random.uniform(0,3)
		time.sleep(wait_time)
		sys.stdout.flush()
		while not pipeline_utils.add_thread_count(global_vars.thread_file, 1):
			time.sleep(1.2)

		for i, input_file in enumerate(self.input()):
			cmd = 'grep "ChrID" %s' % input_file.path #| awk '$17 >= 3' > $file_out
			p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
			# outs, err = p.communicate()
			cmd = "awk '$17>=%s'" % self.cfg['pindel_min_reads']
			p2 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=p1.stdout, shell=True)
			# outs, err = p.communicate()
			outs, err = p2.communicate()
			with open(self.output()[i].path, 'wb') as f:
				f.write(outs)

		while not pipeline_utils.sub_thread_count(global_vars.thread_file, 1):
			time.sleep(1.2)

class parse_pindel(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# case = luigi.Parameter()
	# tumor = luigi.Parameter()
	# matched_n = luigi.Parameter()
	# vcf_path = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return filter_pindel(case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', case_name, 'variants', case_name + '_T.pindel.bed')) for case_name in self.case_dict] + \
		[luigi.LocalTarget(os.path.join(self.project_dir, 'output', case_name, 'variants', case_name + '_N.pindel.bed')) for case_name in self.case_dict if self.case_dict[case_name]['N'] != ''] + \
		[luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'all_samples', 'all_samples_pindel.tsv'))]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)

		wait_time = random.uniform(0,3)
		time.sleep(wait_time)
		sys.stdout.flush()
		while not pipeline_utils.add_thread_count(global_vars.thread_file, 1):
			time.sleep(1.2)

		sample_dict = {output.path.split('/')[-1].split('.pindel.bed')[0]: output.path for output in self.output()[:-1]}

		misc_utils.format_pindel(pindel_files=[input_file.path for input_file in self.input()], sample_dict=sample_dict, project_dir=self.project_dir, all_samples_output=self.output()[-1].path, min_reads=self.cfg['pindel_min_reads'])

		while not pipeline_utils.sub_thread_count(global_vars.thread_file, 1):
			time.sleep(1.2)

class pindel2vcf(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# case = luigi.Parameter()
	# tumor = luigi.Parameter()
	# matched_n = luigi.Parameter()
	# vcf_path = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return pindel(case_dict=self.case_dict, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)

	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'pindel', 'pindel_all_samples.vcf'))

	def run(self):
		pindel_input = '_'.join(self.input()[0].path.split('_')[:-1])
		pipeline_utils.confirm_path(self.output().path)
		cmd = ['./packages/pindel/pindel2vcf', '-r', self.cfg['fasta_file'], '-G', '-R', self.cfg['base_name'], '-d', 'idk', '-P', pindel_input, '-v', self.output().path]
		pipeline_utils.command_call(cmd, [self.output()])

class msisensor(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)


	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, self.case + '_T.msisensor'))
	
	def run(self):
		pipeline_utils.confirm_path(self.output().path)

		cmd = ['./packages/msisensor/binary/msisensor.linux', 'msi', '-d', './packages/msisensor/microsatellites.list', '-t', self.input()[0].path, '-e', self.cfg['library_bed'], '-o', self.output().path] # , '-b', self.max_threads
		pipeline_utils.command_call(cmd, [self.output()]) # , threads_needed=self.max_threads)


class cnvkit(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# case = luigi.Parameter()
	# tumor = luigi.Parameter()
	# matched_n = luigi.Parameter()
	# vcf_path = luigi.Parameter()
	case_dict = luigi.DictParameter()

	# library_bed = luigi.Parameter()
	# fasta_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		# if self.matched_n != '':
		# 	return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads), bam_processing.index_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads)]
		# else:
		return [bam_processing.index_bam(sample=case_name + '_N', fastq_file=self.case_dict[case_name]['N'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg) for case_name in self.case_dict if self.case_dict[case_name]['N'] != ''] \
		+ [bam_processing.index_bam(sample=case_name + '_T', fastq_file=self.case_dict[case_name]['T'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg) for case_name in self.case_dict]


	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'reference', 'reference.cnn'))] \
		+ [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'variants', '%s_T_recalibrated.cnr' % case_name)) for case_name in self.case_dict] \
		+ [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'variants', '%s_T_recalibrated.cns' % case_name)) for case_name in self.case_dict] \
		+ [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'variants', '%s_T_recalibrated.targetcoverage.cnn' % case_name)) for case_name in self.case_dict] \
		+ [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'variants', '%s_T_recalibrated.antitargetcoverage.cnn' % case_name)) for case_name in self.case_dict] \
		+ [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'variants', '%s_N_recalibrated.targetcoverage.cnn' % case_name)) for case_name in self.case_dict if self.case_dict[case_name]['N'] != ''] \
		+ [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'cnvkit', 'variants', '%s_N_recalibrated.antitargetcoverage.cnn' % case_name)) for case_name in self.case_dict if self.case_dict[case_name]['N'] != '']
	
	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		# if self.matched_n:
		# 	cmd = ['./packages/VarDictJava/build/install/VarDict/bin/VarDict', '-G', self.cfg['fasta_file'], '-f', '0.01', '-N', self.case + '_T', '-b', '"%s|%s"' % (self.input()[0][0].path, self.input()[1][0].path), '-z', '-F', '-c', '1', '-S', '2', '-E', '3', '-g', '4', self.cfg['library_bed'], '|', './packages/VarDictJava/VarDict/testsomatic.R', '|', './packages/VarDictJava/VarDict/var2vcf_paired.pl', '-N', '"%s|%s"' % (self.case + '_T', self.case + '_N'), '-f', '0.01', '>%s' % os.path.join(self.vcf_path, 'vardict')]
		# else:
		cmd = ['python3', './packages/cnvkit/cnvkit.py', 'batch', os.path.join(self.project_dir, 'output', '*', 'alignment', '*T*recalibrated.bam'), '--normal', os.path.join(self.project_dir, 'output', '*', 'alignment', '*N*recalibrated.bam'), '--targets', self.cfg['library_bed'], '--fasta', self.cfg['fasta_file'], '--output-reference', self.output()[0].path, '--output-dir', os.path.join(self.project_dir, 'output', 'cnvkit', 'variants'), '--diagram', '--scatter', '--rlibpath', './packages/R', '--annotate', './packages/cnvkit/data/refFlat_b37.txt', '--drop-low-coverage', '-p', self.max_threads]
		pipeline_utils.command_call(cmd, self.output(), threads_needed=self.max_threads)

