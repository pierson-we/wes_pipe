#!/usr/bin/env python3
import fcntl
import os
import sys
import time
import subprocess
import random
import global_vars
import bam_processing
import variant_calling
import pickle

# global_max_threads = 0
# thread_file = ''
# working_files = {}
# cwd = ''

def confirm_path(file):
	wait_time = random.uniform(0,1)
	time.sleep(wait_time)
	if not os.path.exists(os.path.dirname(file)):
		try:
			os.makedirs(os.path.dirname(file))
		except OSError as exc: # Guard against race condition
			if exc.errno != errno.EEXIST:
				raise

def command_call(cmd, outputs, cwd=os.getcwd(), threads_needed=1, sleep_time=1):
	# try:
	wait_time = random.uniform(0,3)
	# print(wait_time)
	time.sleep(wait_time)
	# print('max threads: %s' % global_max_threads)
	# print('thread count: %s' % get_thread_count(thread_file))
	# # while global_max_threads - thread_count < threads_needed:
		
	# # 	time.sleep(sleep_time)
	# # thread_count += threads_needed
	# # print('new thread count: %s' % thread_count)
	# global working_files

	cmd = [str(x) for x in cmd]

	# print(global_vars.thread_file)
	# print(global_vars.global_max_threads)
	# print(global_vars.working_files)
	sys.stdout.flush()
	while not add_thread_count(global_vars.thread_file, threads_needed):
		# print('waiting for godot...')
		time.sleep(sleep_time)
	
	# add_working_files(global_vars.working_files, outputs)
	
	print('\n\n' + ' '.join(cmd) + '\n\n')
	
	subprocess.call(' '.join(cmd), shell=True)
	
	# rm_working_files(global_vars.working_files, outputs)

	while not sub_thread_count(global_vars.thread_file, threads_needed):
		# print('waiting for godot...')
		time.sleep(sleep_time)
	# thread_count -= threads_needed
	# except Exception as e:
	# 	os.chdir(cwd)
	# 	for output in outputs:
	# 		output.remove()
	# 		# print('file_removed')
	# 	# time.sleep(random.uniform(0,0.3))
	# 	os.remove(thread_file)
	# 	raise e


def init_thread_file(thread_file):
	with open(thread_file, 'w') as f:
		f.write('0')

def get_thread_count(thread_file):
	while True:
		try:
			with open(thread_file, 'r') as f:
				fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
				threads = int(f.readlines()[-1])
				fcntl.flock(f, fcntl.LOCK_UN)
				return threads
		except IOError as e:
			time.sleep(0.05)

def add_thread_count(thread_file, threads=1, sleep_time=0.05):
	while True:
		try:
			with open(thread_file, 'r+') as f:
				fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
				thread_count = int(f.readlines()[-1])
				if global_vars.global_max_threads - thread_count >= threads:
					f.write('\n' + str(thread_count + threads))
					fcntl.flock(f, fcntl.LOCK_UN)
					return True
				else:
					fcntl.flock(f, fcntl.LOCK_UN)
					return False
		except IOError as e:
			time.sleep(sleep_time)

def sub_thread_count(thread_file, threads=1, sleep_time=0.05):
	while True:
		try:
			with open(thread_file, 'r+') as f:
				fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
				thread_count = int(f.readlines()[-1])
				new_threads = max(0, thread_count - threads)
				f.write('\n' + str(new_threads))
				fcntl.flock(f, fcntl.LOCK_UN)
				return True
		except IOError as e:
			time.sleep(sleep_time)

def init_working_files(file_log):
	with open(file_log, 'wb') as f:
		pickle.dump({}, f)

def add_working_files(file_log, outputs, sleep_time=0.05):
	read_file = False
	write_file = False
	while not read_file:
		wait_time = random.uniform(0,1)
		time.sleep(wait_time)
		try:
			with open(file_log, 'rb') as f:
				fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
				working_files = pickle.load(f)
				fcntl.flock(f, fcntl.LOCK_UN)
				read_file = True
		except IOError as e:
			time.sleep(sleep_time)

	while not write_file:
		try: 
			with open(file_log, 'wb') as f:
				fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
				for output in outputs:
					working_files[output.path] = ''
				pickle.dump(working_files, f)
				fcntl.flock(f, fcntl.LOCK_UN)
				write_file = True
				print(working_files)
				return
		except IOError as e:
			time.sleep(sleep_time)

def rm_working_files(file_log, outputs, sleep_time=0.05):
	read_file = False
	write_file = False
	while not read_file:
		wait_time = random.uniform(0,3)
		# print(wait_time)
		time.sleep(wait_time)
		try:
			with open(file_log, 'rb') as f:
				fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
				working_files = pickle.load(f)
				fcntl.flock(f, fcntl.LOCK_UN)
				read_file = True
		except IOError as e:
			time.sleep(sleep_time)

	while not write_file:
		try: 
			with open(file_log, 'wb') as f:
				fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
				for output in outputs:
					try:
						del working_files[output.path]
					except:
						pass
				pickle.dump(working_files, f)
				fcntl.flock(f, fcntl.LOCK_UN)
				write_file = True
				return
		except IOError as e:
			time.sleep(sleep_time)
# def error_handling(exception):
# 	global working_files
# 	print('Current working files at time of interruption:')
# 	print(pipeline_utils.working_files)
# 	print(cwd)
# 	os.chdir(cwd)
# 	for file in pipeline_utils.working_files:
# 		os.remove(file)
# 	raise exception