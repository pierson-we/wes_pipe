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

class haplotype_caller(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	sample = luigi.Parameter()
	fastq_file = luigi.Parameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return bam_processing.index_bam(sample=self.sample, fastq_file=self.fastq_file, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)


	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'haplotype_caller', self.sample + '.g.vcf.gz')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'mutect', self.sample + '.g.vcf.gz.tbi'))]
		
	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx4g -Xms4g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'HaplotypeCaller', '-R', self.cfg['fasta_file'], '-I', self.input()[0].path, '-L', self.cfg['library_bed'], '--native-pair-hmm-threads', '1', '-ERC', 'GVCF', '-G', 'Standard', '-G', 'AS_Standard', '-O', self.output()[0].path]
		pipeline_utils.command_call(cmd, self.output())

class consolidate_gvcfs(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case_dict = luigi.DictParameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return [haplotype_caller(sample=case_name + '_N', fastq_file=self.case_dict[case_name]['N'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg) for case_name in self.case_dict if self.case_dict[case_name]['N'] != '']


	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'haplotype_caller', 'genomicsdb'))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		with open('gvcf_sample_map', 'w') as f:
			for gvcf in self.input():
				sample = gvcf[0].path.split('/')[-1].split('.g.vcf')[0]
				f.write('%s\t%s\n' % (sample, gvcf[0].path))
		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx8g -Xms8g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'GenomicsDBImport', '--genomicsdb-workspace-path', os.path.join(self.project_dir, 'output', 'haplotype_caller', 'genomicsdb'), '--batch-size', '50', '-L', self.cfg['library_bed'], '--sample-name-map', 'gvcf_sample_map', '--tmp-dir=%s' % self.cfg['tmp_dir'], '--reader-threads', str(self.max_threads)]
		pipeline_utils.command_call(cmd, [self.output()])

		os.remove('gvcf_sample_map')

class genotype_gvcfs(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case_dict = luigi.DictParameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return consolidate_gvcfs(project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg, case_dict=self.case_dict)


	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'haplotype_caller', 'all_germline.vcf.gz'))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		
		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx4g -Xms4g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'GenotypeGVCFs', '-R', self.cfg['fasta_file'], '-V', 'gendb://' + os.path.join(self.project_dir, 'output', 'haplotype_caller', 'genomicsdb'), '-O', self.output().path, '--tmp-dir=%s' % self.cfg['tmp_dir']]
		pipeline_utils.command_call(cmd, [self.output()])

class filter_snps(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case_dict = luigi.DictParameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return genotype_gvcfs(project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg, case_dict=self.case_dict)


	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'haplotype_caller', 'all_germline_snps_filtered.vcf.gz'))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		
		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx4g -Xms4g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'SelectVariants', '-R', self.cfg['fasta_file'], '-V', self.input().path, '--select-type-to-include', 'SNP', '-O', os.path.join(self.project_dir, 'output', 'haplotype_caller', 'all_germline_snps.vcf.gz')]
		pipeline_utils.command_call(cmd, [self.output()])

		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx4g -Xms4g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'VariantFiltration', '-R', self.cfg['fasta_file'], '-V', os.path.join(self.project_dir, 'output', 'haplotype_caller', 'all_germline_snps.vcf.gz'), '-O', self.output().path, '--filter-expression', '"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"', '--filter-name', '"generic_germline_snp_filter']
		pipeline_utils.command_call(cmd, [self.output()])

		os.remove(os.path.join(self.project_dir, 'output', 'haplotype_caller', 'all_germline_snps.vcf.gz'))

class filter_indels(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case_dict = luigi.DictParameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return genotype_gvcfs(project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg, case_dict=self.case_dict)


	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'haplotype_caller', 'all_germline_indels_filtered.vcf.gz'))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		
		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx4g -Xms4g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'SelectVariants', '-R', self.cfg['fasta_file'], '-V', self.input().path, '--select-type-to-include', 'INDEL', '-O', os.path.join(self.project_dir, 'output', 'haplotype_caller', 'all_germline_indels.vcf.gz')]
		pipeline_utils.command_call(cmd, [self.output()])

		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx4g -Xms4g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'VariantFiltration', '-R', self.cfg['fasta_file'], '-V', os.path.join(self.project_dir, 'output', 'haplotype_caller', 'all_germline_indels.vcf.gz'), '-O', self.output().path, '--filter-expression', '"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"', '--filter-name', '"generic_germline_indel_filter']
		pipeline_utils.command_call(cmd, [self.output()])

		os.remove(os.path.join(self.project_dir, 'output', 'haplotype_caller', 'all_germline_indels.vcf.gz'))

class filter_germline(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case_dict = luigi.DictParameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return [filter_snps(project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg, case_dict=self.case_dict), filter_indels(project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg, case_dict=self.case_dict)]


	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'haplotype_caller', 'all_germline.vcf.gz'))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)
		
		cmd = ['java', '-Xmxg', '-jar', self.cfg['gatk3_location'], '-T', 'CombineVariants', '-R', self.cfg['fasta_file'], '--variant', self.input()[0].path, '--variant', self.input()[1].path, '-o', self.output()[1].path, '-genotypeMergeOptions', 'UNIQUIFY']
		pipeline_utils.command_call(cmd, [self.output()])

