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
import misc_utils
import global_vars
import bam_processing
import variant_calling
import cnv
import germline

class somatic_vcf_intersection(luigi.Task):
	max_threads = luigi.IntParameter()

	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	vcf_path = luigi.Parameter()
	project_dir = luigi.Parameter()
	case_dict = luigi.DictParameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return [variant_calling.filter_mutect(project_dir=self.project_dir, vcf_path=self.vcf_path, case=self.case, tumor=self.tumor, matched_n=self.matched_n, max_threads=self.max_threads, case_dict=self.case_dict, cfg=self.cfg),
		variant_calling.sort_vardict(project_dir=self.project_dir, vcf_path=self.vcf_path, case=self.case, tumor=self.tumor, matched_n=self.matched_n, max_threads=self.max_threads, cfg=self.cfg),
		variant_calling.freebayes(project_dir=self.project_dir, vcf_path=self.vcf_path, case=self.case, tumor=self.tumor, matched_n=self.matched_n, max_threads=self.max_threads, cfg=self.cfg)] 

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.vcf_path, '%s_T_intersect.vcf' % self.case)), luigi.LocalTarget(os.path.join(self.vcf_path, '%s_T_union.vcf' % self.case))]

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		cmd = ['java', '-Xmx2g', '-jar', self.cfg['gatk3_location'], '-T', 'CombineVariants', '-R', self.cfg['fasta_file'], '--variant:mutect', self.input()[0].path, '--variant:vardict', self.input()[1].path, '--variant:freebayes', self.input()[2].path, '-o', self.output()[1].path, '-genotypeMergeOptions', 'PRIORITIZE', '-priority', 'mutect,vardict,freebayes', '--minimumN', '2']
		pipeline_utils.command_call(cmd, self.output())

		cmd = [self.cfg['gatk4_location'], '--java-options', '"-Xmx2g -Xms2g -XX:+UseSerialGC -Djava.io.tmpdir=%s"' % self.cfg['tmp_dir'], 'SelectVariants', '-R', self.cfg['fasta_file'], '-V', self.output()[1].path, '-O', self.output()[0].path, '-select', """'set == "Intersection";'"""]
		pipeline_utils.command_call(cmd, self.output())

class vep(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	# vcf_path = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()

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

	vcf_path = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return [somatic_vcf_intersection(case=self.case, tumor=self.tumor, matched_n=self.matched_n, vcf_path=self.vcf_path, project_dir=self.project_dir, max_threads=self.max_threads, case_dict=self.case_dict, cfg=self.cfg),
		bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg)]
		
	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, '%s_T_fpfilter.vcf' % self.case))

	def run(self):
		pipeline_utils.confirm_path(self.output().path)

		cmd = ['./packages/fpfilter/fpfilter.pl', '--vcf-file', self.input()[0][0].path, '--bam-file', self.input()[1][0].path, '--reference', self.cfg['fasta_file'], '--sample', self.case + '_T', '--output', self.output().path]
		pipeline_utils.command_call(cmd, [self.output()])

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
		# return bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads)
		return [bam_processing.index_bam(sample=case_name + '_N', fastq_file=self.case_dict[case_name]['N'], project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg) for case_name in self.case_dict if self.case_dict[case_name]['N'] != '']

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'normal_bams.txt')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'MSI_BASELINE.txt'))] \
		+ list(itertools.chain(*[[luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'msings', 'baseline', case_name + '_N_recalibrated', case_name + '_N_recalibrated.%s' % file_ext)) for file_ext in ['mpileup', 'msi_output', 'msi.txt']] for case_name in self.case_dict if self.case_dict[case_name]['N'] != '']))
	
	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		# cmd = ['./packages/msings/scripts/create_intervals.sh', './packages/MANTIS/b37_exome_microsatellites.bed']
		# pipeline_utils.command_call(cmd, self.output())

		normal_bams_file = os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'normal_bams.txt')
		with open(normal_bams_file, 'w') as f:
			normal_bams_list = [os.path.join(self.project_dir, 'output', case_name, 'alignment', case_name + '_N_recalibrated.bam') for case_name in self.case_dict if self.case_dict[case_name]['N'] != '']
			f.write('\n'.join(normal_bams_list))
		cmd = ['./packages/msings/scripts/create_baseline.sh', normal_bams_file, './packages/msings/doc/mSINGS_TCGA.msi_intervals', './packages/msings/doc/mSINGS_TCGA.bed', self.cfg['fasta_file'], os.path.join(self.project_dir, 'output', 'msings', 'baseline')]
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
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), bam_processing.index_bam(sample=self.case + '_N', fastq_file=self.matched_n, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), msings_baseline(project_dir=self.project_dir, max_threads=self.max_threads, case_dict=self.case_dict, cfg=self.cfg)]
		else:
			return [bam_processing.index_bam(sample=self.case + '_T', fastq_file=self.tumor, project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg), msings_baseline(project_dir=self.project_dir, max_threads=self.max_threads, case_dict=self.case_dict, cfg=self.cfg)]
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
			cmd = ['python3', './packages/MANTIS/mantis.py', '-b', './packages/msings/doc/mSINGS_TCGA_MANTIS.bed', '--genome', self.cfg['fasta_file'], '-t', self.input()[0][0].path, '-n', self.input()[1][0].path, '-mrq', '20.0', '-mlq', '25.0', '-mlc', '20', '-mrr', '1', '-o', self.output()[0].path]
			pipeline_utils.command_call(cmd, self.output())
		# else:
		# tumor_bams_file = os.path.join(self.project_dir, 'output', 'msings', 'baseline', 'tumor_bams.txt')
		# with open(tumor_bams_file, 'w') as f:
		# 	tumor_bams_list = [os.path.join(self.project_dir, 'output', self.case, 'alignment', case_name + '_T_recalibrated.bam') for case_name in self.case_dict if self.case_dict[case_name]['N'] == '']
		# 	f.write('\n'.join(tumor_bams_list))
	
		cmd = ['./packages/msings/scripts/run_msings_single_sample.sh', self.input()[0][0].path, './packages/msings/doc/mSINGS_TCGA.msi_intervals', './packages/msings/doc/mSINGS_TCGA.bed', self.cfg['fasta_file'], './packages/msings/doc/mSINGS_TCGA.baseline', os.path.join(self.project_dir, 'output', 'msings', 'tumor')]
		pipeline_utils.command_call(cmd, self.output())
		os.rename(os.path.join(self.project_dir, 'output', 'msings', 'tumor', self.case + '_T_recalibrated', self.case + '_T_recalibrated.MSI_Analysis.txt'), os.path.join(self.vcf_path, self.case + '_msings.txt'))
		# cmd = ['echo', '"mSINGS still needs to be set up for tumor-only samples"', '>', self.output()[0].path] # this will be a pain to get up and going: https://bitbucket.org/uwlabmed/msings/src/8269e0e01acfc5e01d0de9d63ffc1e399996ce8a/Recommendations_for_custom_assays?at=master&fileviewer=file-view-default
		# pipeline_utils.command_call(cmd, self.output())

class vcf2maf(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	vcf_path = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()
	case_dict = luigi.DictParameter()

	cfg = luigi.DictParameter()

	def requires(self):
		return fpfilter(max_threads=self.max_threads, project_dir=self.project_dir, vcf_path=self.vcf_path, case=self.case, tumor=self.tumor, matched_n=self.matched_n, case_dict=self.case_dict, cfg=self.cfg)

	def output(self):
		return luigi.LocalTarget(os.path.join(self.vcf_path, '%s_T.maf' % self.case))

	def run(self):
		if self.matched_n != '':
			cmd = ['./packages/misc/mskcc-vcf2maf-decbf60/vcf2maf.pl', '--input-vcf', self.input().path, '--output-maf', self.output().path, '--tumor-id', self.case + '_T', '--vcf-tumor-id', self.case + '_T', '--normal-id', self.case + '_N', '--vcf-normal-id', self.case + '_N', '--vep-path', './packages/ensembl-vep', '--vep-data', './packages/ensembl-vep/cache', '--ref-fasta', self.cfg['fasta_file'], '--species', 'homo_sapiens', '--ncbi-build', 'GRCh38', '--cache-version', '94', '--filter-vcf', '0']
		else:
			cmd = ['./packages/misc/mskcc-vcf2maf-decbf60/vcf2maf.pl', '--input-vcf', self.input().path, '--output-maf', self.output().path, '--tumor-id', self.case + '_T', '--vcf-tumor-id', self.case + '_T', '--vep-path', './packages/ensembl-vep', '--vep-data', './packages/ensembl-vep/cache', '--ref-fasta', self.cfg['fasta_file'], '--species', 'homo_sapiens', '--ncbi-build', 'GRCh38', '--cache-version', '94', '--filter-vcf', '0']
		pipeline_utils.command_call(cmd, [self.output()])

class combine_mafs(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case_dict = luigi.DictParameter()
	cfg = luigi.DictParameter()

	def requires(self):
		return [vcf2maf(max_threads=self.max_threads, project_dir=self.project_dir, case=case, tumor=self.case_dict[case]['T'], matched_n=self.case_dict[case]['N'], case_dict=self.case_dict, cfg=self.cfg, vcf_path=os.path.join(self.project_dir, 'output', case, 'variants')) for case in self.case_dict]

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'all_samples', 'all_samples.maf'))] + self.input()

	def run(self):
		pipeline_utils.confirm_path(self.output()[0].path)
		misc_utils.combine_mafs([maf_input.path for maf_input in self.input()], self.output()[0].path)

class combine_cnvs(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case_dict = luigi.DictParameter()
	cfg = luigi.DictParameter()

	def requires(self):
		return [cnv.refine_cnv(case=case, max_threads=self.max_threads, project_dir=self.project_dir, cfg=self.cfg, case_dict=self.case_dict) for case in self.case_dict]

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'all_samples', 'all_samples_cnv.tsv'))] + self.input()

	def run(self):
		pipeline_utils.confirm_path(self.output()[0].path)
		misc_utils.combine_cnvs([cnv_input[-1].path for cnv_input in self.input()], [case for case in self.case_dict], self.output()[0].path)

class create_mut_mats(luigi.Task):
	max_threads = luigi.IntParameter()
	project_dir = luigi.Parameter()

	case_dict = luigi.DictParameter()
	cfg = luigi.DictParameter()

	def requires(self):
		return [combine_mafs(max_threads=self.max_threads, project_dir=self.project_dir, cfg=self.cfg, case_dict=self.case_dict), 
		combine_cnvs(max_threads=self.max_threads, project_dir=self.project_dir, cfg=self.cfg, case_dict=self.case_dict),
		variant_calling.annotate_pindel(project_dir=self.project_dir, max_threads=self.max_threads, cfg=self.cfg, case_dict=self.case_dict)]

	def output(self):
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'all_samples', 'mut_mat.tsv')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'all_samples', 'cnv_mat.tsv')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'all_samples', 'mut_counts.tsv'))]
		# return self.input()

	def run(self):
		for output in self.output():
			pipeline_utils.confirm_path(output.path)
		mafs = [input_file.path for input_file in self.input()[0][1:]]
		cnvs = [input_files[-1].path for input_files in self.input()[1][1:]]
		pindel = [input_file.path for input_file in self.input()[2]]
		misc_utils.create_mut_mats(mafs=mafs, cnvs=cnvs, pindel=pindel, mut_mat_file=self.output()[0].path, cnv_mat_file=self.output()[1].path, mut_counts_file=self.output()[2].path)

class cases(luigi.Task):
	# generated parameters
	sample_dict = luigi.DictParameter()
	project_dir = luigi.Parameter()
	sample_threads = luigi.IntParameter()

	# cfg parameters
	fasta_file = luigi.Parameter()
	library_bed = luigi.Parameter()
	gatk4_location = luigi.Parameter()
	gatk3_location = luigi.Parameter()
	known_vcf = luigi.Parameter()
	germline_resource = luigi.Parameter()
	picard_location = luigi.Parameter()
	vardict_location = luigi.Parameter()
	mills = luigi.Parameter()
	kg = luigi.Parameter()
	#omni = luigi.Parameter()
	#hapmap = luigi.Parameter()
	library_prep = luigi.Parameter()
	platform = luigi.Parameter()
	base_name = luigi.Parameter()
	samtools_location = luigi.Parameter()
	bowtie_build_location = luigi.Parameter()
	bowtie_location = luigi.Parameter()
	fastqc_location = luigi.Parameter()
	trim_location = luigi.Parameter()
	insert_size = luigi.Parameter()
	freebayes_location = luigi.Parameter()
	vcffilter_location = luigi.Parameter()
	cnvkit_location = luigi.Parameter()
	refFlat = luigi.Parameter()
	cnvkit_seg_method = luigi.Parameter()
	cnvkit_genemetrics_threshold = luigi.Parameter()
	cnvkit_genemetrics_minprobes = luigi.Parameter()
	pindel_min_reads = luigi.IntParameter()
	pindel_min_qual = luigi.IntParameter()
	pindel_max_inv_length = luigi.IntParameter()
	genmap = luigi.Parameter()
	exons_bed = luigi.Parameter()

	def requires(self):
		cfg = {
			'fasta_file': self.fasta_file,
			'library_bed': self.library_bed,
			'gatk4_location': self.gatk4_location,
			'gatk3_location': self.gatk3_location,
			'known_vcf': self.known_vcf,
			'germline_resource': self.germline_resource,
			'picard_location': self.picard_location,
			'vardict_location': self.vardict_location,
			'mills': self.mills,
			'kg': self.kg,
			#'omni': self.omni,
			#'hapmap': self.hapmap,
			'library_prep': self.library_prep,
			'platform': self.platform,
			'base_name': self.base_name,
			'samtools_location': self.samtools_location,
			'bowtie_build_location': self.bowtie_build_location,
			'bowtie_location': self.bowtie_location,
			'fastqc_location': self.fastqc_location,
			'trim_location': self.trim_location,
			'insert_size': self.insert_size,
			'freebayes_location': self.freebayes_location,
			'vcffilter_location': self.vcffilter_location,
			'cnvkit_location': self.cnvkit_location,
			'refFlat': self.refFlat,
			'cnvkit_seg_method': self.cnvkit_seg_method,
			'cnvkit_genemetrics_threshold': self.cnvkit_genemetrics_threshold,
			'cnvkit_genemetrics_minprobes': self.cnvkit_genemetrics_minprobes,
			'pindel_min_reads': self.pindel_min_reads,
			'pindel_min_qual': self.pindel_min_qual,
			'pindel_max_inv_length': self.pindel_max_inv_length,
			'genmap': self.genmap,
			'exons_bed': self.exons_bed,
			'tmp_dir': os.path.join(self.project_dir, 'tmp')
		}
		pipeline_utils.confirm_path(cfg['tmp_dir'])
		# return [aggregate_variants(case=case, tumor=self.sample_dict[case]['T'], matched_n=self.sample_dict[case]['N'], project_dir=self.project_dir, max_threads=self.sample_threads, case_dict=self.sample_dict) for case in self.sample_dict]

		return [create_mut_mats(max_threads=self.sample_threads, project_dir=self.project_dir, cfg=cfg, case_dict=self.sample_dict)] + \
		[germline.filter_germline(project_dir=self.project_dir, max_threads=self.sample_threads, cfg=cfg, case_dict=self.sample_dict)] + \
		[variant_calling.msisensor(max_threads=self.sample_threads, project_dir=self.project_dir, case=case, tumor=self.sample_dict[case]['T'], matched_n=self.sample_dict[case]['N'], cfg=cfg, vcf_path=os.path.join(self.project_dir, 'output', case, 'variants')) for case in self.sample_dict]
		# [variant_calling.annotate_pindel(project_dir=self.project_dir, max_threads=self.sample_threads, cfg=cfg, case_dict=self.sample_dict)]
		# [msi(case=case, tumor=self.sample_dict[case]['T'], matched_n=self.sample_dict[case]['N'], project_dir=self.project_dir, max_threads=self.sample_threads, case_dict=self.sample_dict, cfg=cfg, vcf_path=os.path.join(self.project_dir, 'output', case, 'variants')) for case in self.sample_dict]

		# [variant_analysis.vep(case=case, tumor=self.sample_dict[case]['T'], matched_n=self.sample_dict[case]['N'], project_dir=self.project_dir, max_threads=self.sample_threads, case_dict=self.sample_dict, cfg=cfg) for case in self.sample_dict] + \
		
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
		# return luigi.LocalTarget(os.path.join(self.project_dir, 'output', 'all_samples.vcf.gz'))
		return self.input()

	# def run(self):
	# 	input_vcfs = []
	# 	for variant_caller in self.input():
	# 		if not isinstance(variant_caller, list):
	# 			vcf_path = variant_caller.path
	# 			if vcf_path[-3:] != '.gz':
	# 				# with gzip.open(vcf_path, 'rb') as f:
	# 				# 	with open(vcf_path.split('.gz')[0], 'wb') as new_f:
	# 				# 		new_f.write(f.read())
	# 				# vcf_path = vcf_path.split('.gz')[0]
	# 				cmd = ['bgzip', vcf_path]
	# 				pipeline_utils.command_call(cmd, [vcf_path + '.gz'], sleep_time=0.05)
	# 				vcf_path = vcf_path + '.gz'
	# 			if not os.path.exists(vcf_path + '.tbi'):
	# 				cmd = ['tabix', vcf_path]
	# 				pipeline_utils.command_call(cmd, [vcf_path], sleep_time=0.05)
	# 			input_vcfs.append(vcf_path)
	# 	cmd = ['bcftools merge %s | bgzip -c > %s' % (' '.join(input_vcfs), self.output().path)]
	# 	pipeline_utils.command_call(cmd, [self.output()], sleep_time=1)
	# 			variant_analysis.vep(max_threads=self.max_threads, project_dir=self.project_dir, case=self.case, tumor=self.tumor, matched_n=self.matched_n, vcf_path=vcf_path, cfg=self.cfg)
