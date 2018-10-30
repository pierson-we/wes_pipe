#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import time
import luigi
# from code import pipeline_utils, global_vars

def run_pipeline(args):
	import variant_analysis
	import global_vars
	import pipeline_utils
	timestamp = str(int(time.time()))

	global_vars.global_max_threads = args.max_threads
	global_vars.thread_file = os.path.join(os.getcwd(), 'thread_count_temp_%s.txt' % timestamp)
	print(global_vars.thread_file)
	pipeline_utils.init_thread_file(global_vars.thread_file)
	global_vars.working_files = os.path.join(os.getcwd(), 'working_files_%s.pkl' % timestamp)
	# pipeline_utils.init_working_files(global_vars.working_files)
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
	# with open('size_of_vars.txt', 'w') as f:
	# 	f.write('vars().values() is %s' % sys.getsizeof(vars().values()))
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
	luigi.build([variant_analysis.cases(sample_dict=sample_dict, project_dir=args.project_dir, sample_threads=sample_threads)], workers=args.workers, local_scheduler=args.local_scheduler, worker_scheduler_factory=worker_scheduler_factory) # , workers=args.workers #, scheduler_port=int(args.port)) # workers=sample_threads

		# [(max_threads=args.max_threads, project_dir=args.project_dir, sample_dir=args.sample_dir, threads_per_sample=args.threads_per_sample, timestamp=timestamp)], workers=args.workers, local_scheduler=args.local_scheduler)
		# yield aggregate_variants(case=case, tumor=tumor, matched_n=matched_n, project_dir=self.project_dir, max_threads=sample_threads, case_dict=sample_dict)

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
	parser.add_argument('--port', action='store', dest='port', default='8082', help='If using the central luigi scheduler, use this parameter to specify a custom port for the luigi server to operate on (defaults to 8082)')

	args = parser.parse_args()

	# make sure we're in the correct working directory so relative references work. If not, change to the correct directory
	if not os.path.exists(os.path.join(os.getcwd(), 'code')):
		os.chdir('/'.join(sys.argv[0].split('/')[:-1]))
		sys.path.append(os.getcwd())
	if not os.path.exists(os.path.join(os.getcwd(), 'code')):
		raise ValueError('you must run script from "wes_pipe" directory, as relative references are used throughout the analysis.')
	
	sys.path.append('./code')
	sys.path.append('./packages')
	#sys.path.append('./packages/lib/python3.6/site-packages')
	# import luigi
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
	# if not args.local_scheduler:
		# print('about to start luigi server...')
		# subprocess.call('python3 /home/wpierson/luigi/src/luigi/bin/luigid --background --pidfile /home/wpierson/projects/wes_pipe/luigi_pidfile.txt --logdir /home/wpierson/projects/wes_pipe --state-path /home/wpierson/projects/wes_pipe/luigi_statepath.txt --port %s &' % args.port)
		# print('Starting luigi server...\n\n')
		# sys.stdout.flush()
		# time.sleep(2)

	run_pipeline(args)
	# luigi.build([bam_processing.cases(max_threads=args.max_threads, project_dir=args.project_dir, sample_dir=args.sample_dir, threads_per_sample=args.threads_per_sample, timestamp=timestamp)], workers=args.workers, local_scheduler=args.local_scheduler)
		
	# luigi.build([bowtie(fastq_path=fastq_path, sam_path=sam_path, threads=threads, fasta_path=fasta_path), convert_bam(sam_path=sam_path, bam_path=bam_path)], workers=1, local_scheduler=False)
