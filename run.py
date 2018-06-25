#!/usr/bin/env python3
import os
import sys
import argparse
from code import pipeline_utils, global_vars

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='wes pipeline parser')
	parser.add_argument('-O', action='store', dest='project_dir', default=os.getcwd(), help='Directory in which the program will create an "output" directory containing all output files.')
	parser.add_argument('-I', action='store', dest='sample_dir', help='Directory containing the input fastq files. Directory contents must be structured as follows:\n[sample_dir]\n\t[sample_name]\n\t\ttumor\n\t\t\t[tumor].fastq.gz\n\t\tnormal\n\t\t\t[normal].fastq.gz')
	parser.add_argument('--threads', action='store', dest='max_threads', default=4, type=int, help='The maximum number of threads available for use. The program will attempt to distribute available threads amongst samples as efficiently and equally as possible.')
	parser.add_argument('--threads_per_sample', action='store', dest='threads_per_sample', default=0, type=int, help='The number of threads to be allowed per sample. If not specified, the program will divide the threads evenly among samples. Specifying it allows you to prioritize single sample completion over equal prcessing of samples simultaneously.')
	parser.add_argument('--workers', action='store', dest='workers', default=1, type=int, help='The number of workers that should be used by the pipeline scheduler (Luigi) to schedule jobs - it will not necessarily determine the number of threads utilized. Recommended: set equal to the number of sequencing files being processed.')
	parser.add_argument('-m', action='store_false', dest='mutect', default=True, help='This flag suppresses analysis with Mutect2')
	parser.add_argument('-d', action='store_false', dest='vardict', default=True, help='This flag suppresses analysis with VarDict')
	parser.add_argument('-f', action='store_false', dest='freebayes', default=True, help='This flag suppresses analysis with FreeBayes')
	parser.add_argument('-v', action='store_false', dest='varscan', default=True, help='This flag suppresses analysis with VarScan')
	parser.add_argument('-c', action='store_false', dest='cnvkit', default=True, help='This flag suppresses analysis with CNVKit')
	parser.add_argument('-s', action='store_false', dest='scalpel', default=True, help='This flag suppresses analysis with Scalpel')
	parser.add_argument('-l', action='store_false', dest='lumpy', default=True, help='This flag suppresses analysis with LUMPY')
	parser.add_argument('-D', action='store_false', dest='delly', default=True, help='This flag suppresses analysis with DELLY')
	parser.add_argument('-w', action='store_false', dest='wham', default=True, help='This flag suppresses analysis with WHAM')
	parser.add_argument('--local_scheduler', action='store_true', dest='local_scheduler', default=False, help='This flag will use the local luigi scheduler as opposed to the luigi server.')

	args = parser.parse_args()

	if not os.path.exists(os.path.join(os.getcwd(), 'code')):
		os.chdir('/'.join(sys.argv[0].split('/')[:-1]))
		sys.path.append(os.getcwd())
	if not os.path.exists(os.path.join(os.getcwd(), 'code')):
		raise ValueError('you must run script from "wes_pipe" directory, as relative references are used throughout the analysis.')
	
	sys.path.append('./code')
	sys.path.append('./packages')
	import luigi
	import bam_processing
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
	
	# project_dir = sys.argv[1]
	# sample_dir = sys.argv[2]
	
	luigi.build([bam_processing.cases(max_threads=args.max_threads, project_dir=args.project_dir, sample_dir=args.sample_dir, threads_per_sample=args.threads_per_sample)], workers=args.workers, local_scheduler=args.local_scheduler)
		
	# luigi.build([bowtie(fastq_path=fastq_path, sam_path=sam_path, threads=threads, fasta_path=fasta_path), convert_bam(sam_path=sam_path, bam_path=bam_path)], workers=1, local_scheduler=False)
