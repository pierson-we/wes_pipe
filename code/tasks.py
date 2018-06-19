#!/usr/bin/env python3
import subprocess
import luigi
import luigi.interface
import os
import sys
import time

# ***
# GATK Resource Bundle Files
# https://software.broadinstitute.org/gatk/download/bundle
# ***

def confirm_path(file):
	if not os.path.exists(os.path.dirname(file)):
		try:
			os.makedirs(os.path.dirname(file))
		except OSError as exc: # Guard against race condition
			if exc.errno != errno.EEXIST:
				raise

def command_call(cmd, threads_needed=1, sleep_time=1):
	global thread_count, max_threads
	while max_threads - thread_count < max_threads:
		time.sleep(sleep_time)
	thread_count += threads_needed
	subprocess.call(cmd, stdout=subprocess.PIPE)
	thread_count -= threads_needed

class genome_index(luigi.Task):
	max_threads = luigi.Parameter()
	fasta_file = luigi.Parameter()
	# threads = luigi.Parameter()
	base_name = luigi.Parameter()

	fasta_dir = os.path.join('/', *luigi.Parameter().task_value('genome_index', 'fasta_file').split('/')[:-1])

	def output(self):
		return luigi.LocalTarget(os.path.join(self.fasta_dir, 'index', self.base_name + '.1.bt2'))
	
	def run(self):
		cwd = os.getcwd()
		os.chdir(self.fasta_dir)
		if not os.path.exists('./index'):
			os.mkdir('./index')
		os.chdir('./index')
		subprocess.call(['bowtie2-build', '--threads=%s' % self.max_threads, self.fasta_file, self.base_name], stdout=subprocess.PIPE)

		os.chdir(cwd)

class samtools_index(luigi.Task):
	max_threads = luigi.Parameter()
	fasta_file = luigi.Parameter()

	def output(self):
		return luigi.LocalTarget(self.fasta_file + '.fai')
	
	def run(self):
		subprocess.call(['samtools', 'faidx', '--nthreads=%s' % self.max_threads, self.fasta_file], stdout=subprocess.PIPE)

		# samtools faidx Homo_sapiens_assembly18.fasta

class picard_index(luigi.Task):
	fasta_file = luigi.Parameter()
	picard_location = luigi.Parameter()

	def output(self):
		return luigi.LocalTarget(self.fasta_file + '.dict')
	
	def run(self):
		cmd = ['java', '-jar', self.picard_location, 'CreateSequenceDictionary', 'R=%s' % self.fasta_file, 'O=%s' % self.output().path]
		subprocess.call(cmd, stdout=subprocess.PIPE)

# FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# FASTX-Toolkit: http://hannonlab.cshl.edu/fastx_toolkit/commandline.html
# class quality_control(luigi.Task):
# 	fastq_file = luigi.Parameter()
# 	fastqc_file = luigi.Parameter()

# 	def output(self):
# 		return luigi.LocalTarget(self.fastqc_file)

# 	def run(self):
# 		confirm_path(self.fastqc_file)


class bowtie(luigi.Task):
	max_threads = luigi.Parameter()
	fastq_file = luigi.Parameter()
	# sam_file = luigi.Parameter()
	threads = luigi.Parameter()
	fasta_file = luigi.Parameter()
	base_name = luigi.Parameter()
	sample = luigi.Parameter()
	# fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	# fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])

	def requires(self):
		yield genome_index(max_threads=self.max_threads) #threads=self.threads, base_name=self.base_name, fasta_path=self.fasta_path)
		yield samtools_index(max_threads=self.max_threads)
		yield picard_index()

	def output(self):
		print(self.fasta_file.split('/')[:-1])
		return luigi.LocalTarget(os.path.join(*['/', self.fasta_file.split('/')[:-1], 'index', self.sample + '_raw.sam']))

	def run(self):
		# try:
		cwd = os.getcwd()
		os.chdir(os.path.join(*['/', self.fasta_file.split('/')[:-1], 'index']))

		subprocess.call(['bowtie2', '-x', self.base_name, '--threads=%s' % self.threads, '-U', self.fastq_file, '-S', self.output().path], stdout=subprocess.PIPE)

		os.chdir(cwd)
		# except KeyboardInterrupt:
		# 	sys.exit()

# class move_file(luigi.Task):
# 	from_file = luigi.Parameter()
# 	to_file = luigi.Parameter()
# 	required = luigi.Parameter()

# 	def requires(self):
# 		return self.required
	
# 	def output(self):
# 		return luigi.LocalTarget(self.to_file)

# 	def run(self):
# 		os.rename(self.from_file, self.to_file)

# class convert_bam(luigi.Task):
# 	project_dir = luigi.Parameter()
# 	fastq_file = luigi.Parameter()

# 	# case_dir = luigi.Parameter()
# 	sample = luigi.Parameter()
# 	# fastq_file = luigi.Parameter()

# 	# sam_file = luigi.Parameter() # os.path.join(luigi.Parameter().task_value('convert_bam', 'case_dir'), 'alignment', luigi.Parameter().task_value('convert_bam', 'sample'), '.sam')
# 	# bam_file = luigi.Parameter() # os.path.join(luigi.Parameter().task_value('convert_bam', 'case_dir'), 'alignment', luigi.Parameter().task_value('convert_bam', 'sample'), '.bam')

# 	def requires(self):
# 		return bowtie(fastq_file=self.fastq_file, sample=self.sample)

# 	def output(self):
# 		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.sample + '_raw.bam'))

# 	def run(self):
# 		try:
# 			confirm_path(self.output().path)
# 			p = subprocess.Popen("samtools view -Sb %s > %s" % (self.input().path, self.output().path), shell=True, stdout=subprocess.PIPE)
# 			p.wait()
# 			#os.remove(self.input().path)
# 		except KeyboardInterrupt:
# 			try:
# 				os.remove(self.output().path)
# 			except:
# 				pass

class add_read_groups(luigi.Task):
	max_threads = luigi.Parameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	picard_location = luigi.Parameter()
	library_prep = luigi.Parameter()
	platform = luigi.Parameter()

	# input_bam = luigi.Parameter()
	# output_bam = luigi.Parameter()
	sample = luigi.Parameter()

	def requires(self):
		return bowtie(fastq_file=self.fastq_file, sample=self.sample, max_threads=self.max_threads)

	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.sample + '_preprocessed.bam'))

	def run(self):
		confirm_path(self.output().path)
		cmd = ['java', '-jar', self.picard_location, 'AddOrReplaceReadGroups', 'I=%s' % self.input().path, 'O=%s' % self.output().path, 'SORT_ORDER=coordinate', 'RGID=1', 'RGLB=%s' % self.library_prep, 'RGPL=%s' % self.platform, 'RGPU=barcode', 'RGSM=%s' % self.sample]
		subprocess.call(cmd, stdout=subprocess.PIPE)

# # https://broadinstitute.github.io/picard/command-line-overview.html
# class sort_bam(luigi.Task):
# 	project_dir = luigi.Parameter()
# 	fastq_file = luigi.Parameter()
# 	picard_location = luigi.Parameter()

# 	# input_bam = luigi.Parameter()
# 	# output_bam = luigi.Parameter()
# 	sample = luigi.Parameter()

# 	def requires(self):
# 		return add_read_groups(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample)

# 	def output(self):
# 		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.sample + '_sorted.bam'))

# 	def run(self):
# 		confirm_path(self.output().path)
# 		cmd = ['java', '-jar', self.picard_location, 'SortSam', 'I=%s' % self.input().path, 'O=%s' % self.output().path, 'SORT_ORDER=coordinate']
# 		subprocess.call(cmd, stdout=subprocess.PIPE)

# https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
class mark_duplicates(luigi.Task):
	max_threads = luigi.Parameter()
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
		return [luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.sample + '_marked_dups.bam')), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.sample + '_marked_dups_metrics.txt'))]

	def run(self):
		confirm_path(self.output()[0].path)
		cmd = ['java', '-jar', self.picard_location, 'MarkDuplicatesWithMateCigar', 'I=%s' % self.input().path, 'O=%s' % self.output()[0].path, 'M=%s' % self.output()[1].path]
		subprocess.call(cmd, stdout=subprocess.PIPE)

class index_bam(luigi.Task):
	max_threads = luigi.Parameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()

	sample = luigi.Parameter()

	# input_bam = luigi.Parameter()
	# recalibration_report = luigi.Parameter()
	# output_bam = luigi.Parameter()

	def requires(self):
		return mark_duplicates(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads)

	def output(self):

		return [self.input()[0], luigi.LocalTarget(self.input()[0].path + '.bai')]

	def run(self):
		confirm_path(self.output()[1].path)
		cmd = ['samtools', 'index', '-b', '--nthreads=%s' % self.max_threads, self.input()[0].path]
		subprocess.call(cmd, stdout=subprocess.PIPE)

# ~20 mins w/2 cores
class realigner_target(luigi.Task):
	max_threads = luigi.Parameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	gatk3_location = luigi.Parameter()

	sample = luigi.Parameter()
	fasta_file = luigi.Parameter()
	known_vcf = luigi.Parameter()

	# input_bam = luigi.Parameter()
	# output_file = luigi.Parameter()

	def requires(self):
		return index_bam(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads)

	def output(self):
		return [self.input()[0], luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.sample + '_forIndelRealigner.intervals'))]

	def run(self):
		confirm_path(self.output()[1].path)
		cmd = ['java', '-jar', self.gatk3_location, '-nt', self.max_threads, '-T', 'RealignerTargetCreator', '-R', self.fasta_file, '-I', self.input()[0].path, '--known', self.known_vcf, '-o', self.output()[1].path]
		subprocess.call(cmd, stdout=subprocess.PIPE)

# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php
# ~60 mins (no multiprocessing support)
class indel_realignment(luigi.Task):
	max_threads = luigi.Parameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	gatk3_location = luigi.Parameter()

	fasta_file = luigi.Parameter()
	known_vcf = luigi.Parameter()
	sample = luigi.Parameter()
	
	# input_bam = luigi.Parameter()
	# output_bam = luigi.Parameter()

	def requires(self):
		return realigner_target(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads)

	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.sample + '_realigned.bam'))

	def run(self):
		confirm_path(self.output().path)
		cmd = ['java', '-jar', self.gatk3_location, '-T', 'IndelRealigner', '-R', self.fasta_file, '-I', self.input()[0].path, '-known', self.known_vcf, '-targetIntervals', self.input()[1].path, '-o', self.output().path]
		subprocess.call(cmd, stdout=subprocess.PIPE)

# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
class bqsr(luigi.Task):
	max_threads = luigi.Parameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	gatk3_location = luigi.Parameter()

	fasta_file = luigi.Parameter()
	known_vcf = luigi.Parameter()
	sample = luigi.Parameter()
	
	# input_bam = luigi.Parameter()
	# output_table = luigi.Parameter()

	def requires(self):
		return indel_realignment(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads)

	def output(self):
		return [self.input(), luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.sample + '_recalibrated.table'))]

	def run(self):
		confirm_path(self.output()[1].path)
		cmd = ['java', '-jar', self.gatk3_location, '-nt', self.max_threads, '-T', 'BaseRecalibrator', '-R', self.fasta_file, '-I', self.input().path, '-knownSites', self.known_vcf, '-o',  self.output()[1].path]
		subprocess.call(cmd, stdout=subprocess.PIPE)

class recalibrated_bam(luigi.Task):
	max_threads = luigi.Parameter()
	project_dir = luigi.Parameter()
	fastq_file = luigi.Parameter()
	gatk3_location = luigi.Parameter()

	fasta_file = luigi.Parameter()
	sample = luigi.Parameter()

	# input_bam = luigi.Parameter()
	# recalibration_report = luigi.Parameter()
	# output_bam = luigi.Parameter()

	def requires(self):
		return bqsr(fastq_file=self.fastq_file, project_dir=self.project_dir, sample=self.sample, max_threads=self.max_threads)

	def output(self):
		return luigi.LocalTarget(os.path.join(self.project_dir, 'output', self.sample[:-2], self.sample + '_recalibrated.bam'))

	def run(self):
		confirm_path(self.output().path)
		cmd = ['java', '-jar', self.gatk3_location, '-nt', self.max_threads, '-T', 'PrintReads', '-R', self.fasta_file, '-I', self.input()[0].path, '-BQSR', self.input()[1].path, '-o',  self.output().path]
		subprocess.call(cmd, stdout=subprocess.PIPE)

class run_variant_caller(luigi.Task):
	max_threads = luigi.Parameter()
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
	fasta_dir = os.path.join('/', *luigi.Parameter().task_value('bowtie', 'fasta_file').split('/')[:-1])
	
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
	max_threads = luigi.Parameter()
	project_dir = luigi.Parameter()
	case = luigi.Parameter()
	tumor = luigi.Parameter()
	matched_n = luigi.Parameter()

	
	def requires(self):
		case_dir = os.path.join(self.project_dir, self.case)
		vcf_path = os.path.join(case_dir, 'variants')

		kwargs = {
		'project_dir': self.project_dir,
		'case_dir': case_dir,
		'vcf_path': vcf_path,
		'case': self.case,
		'tumor': self.tumor,
		'matched_n': self.matched_n,
		'max_threads': self.max_threads
		}
		yield run_variant_caller(caller='MuTect', **kwargs)
		yield run_variant_caller(caller='VarDict', **kwargs)
		yield run_variant_caller(caller='FreeBayes', **kwargs)
		yield run_variant_caller(caller='VarScan', **kwargs)

		yield run_variant_caller(caller='Scalpel', **kwargs)

		yield run_variant_caller(caller='CNVkit', **kwargs)

		yield run_variant_caller(caller='LUMPY', **kwargs)
		yield run_variant_caller(caller='DELLY', **kwargs)
		yield run_variant_caller(caller='WHAM', **kwargs)

	def output(self):
		return self.input()

class cases(luigi.Task):
	# sample_csv = luigi.Parameter()
	project_dir = '/Users/wep/Documents/Research/Rare_Tumors/pipeline/'
	max_threads = luigi.Parameter()

	def requires(self):
		# sample_df = pd.read_csv(sample_csv, header=True, index_col='sample_id')
		# sample_dict = {}

		# for sample in sample_df.index.tolist():
		# 	case = sample_df.iloc[sample]['case']
		# 	if case not in sample_dict:
		# 		sample_dict[case] = {'T':'', 'N':''}
		# 	sample_type = sample_df.iloc[sample]['type']
		# 	sample_dict[case][sample_type] = sample_df.iloc[sample]['file']

		sample_dict = {'ERR031838_1': {'T': '/Users/wep/Documents/Research/Rare_Tumors/pipeline/test_data/ERR031838_1.fastq.gz', 'N': ''}}
		for case in sample_dict:
			tumor = sample_dict[case]['T']
			matched_n = sample_dict[case]['N']
			yield aggregate_variants(case=case, tumor=tumor, matched_n=matched_n, project_dir=self.project_dir, max_threads=self.max_threads)

if __name__ == '__main__':
	# sample_csv = sys.argv[1]
	# sample_df = pd.read_csv(sample_csv, header=True, index_col='sample_id')
	# sample_dict = {}

	# for sample in sample_df.index.tolist():
	# 	case = sample_df.iloc[sample]['case']
	# 	if case not in sample_dict:
	# 		sample_dict[case] = {'T':'', 'N':''}
	# 	sample_type = sample_df.iloc[sample]['type']
	# 	sample_dict[case][sample_type] = sample_df.iloc[sample]['file']

	# sample_dict = {'ERR031838_1': {'T': '/Users/wep/Documents/Research/Rare_Tumors/pipeline/test_data/ERR031838_1.fastq.gz', 'N': ''}}
	# for case in sample_dict:
	# 	tumor = sample_dict[case]['T']
	# 	matched_n = sample_dict[case]['N']
	# 	luigi.build([variant_calls(case=case, tumor=tumor, matched_n=matched_n)], workers=1, local_scheduler=False)
	luigi.build([cases(max_threads=2)], workers=1, local_scheduler=False)
	# luigi.build([bowtie(fastq_path=fastq_path, sam_path=sam_path, threads=threads, fasta_path=fasta_path), convert_bam(sam_path=sam_path, bam_path=bam_path)], workers=1, local_scheduler=False)
