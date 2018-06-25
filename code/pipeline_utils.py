#!/usr/bin/env python3
import fcntl
import os
import time
import subprocess
import random

global_max_threads = 0
thread_file = ''

def confirm_path(file):
	if not os.path.exists(os.path.dirname(file)):
		try:
			os.makedirs(os.path.dirname(file))
		except OSError as exc: # Guard against race condition
			if exc.errno != errno.EEXIST:
				raise

def command_call(cmd, outputs, cwd=os.getcwd(), threads_needed=1, sleep_time=1):
	# try:
	# wait_time = random.uniform(0,0.3)
	# print(wait_time)
	# time.sleep(wait_time)
	# print('max threads: %s' % global_max_threads)
	# print('thread count: %s' % get_thread_count(thread_file))
	# # while global_max_threads - thread_count < threads_needed:
		
	# # 	time.sleep(sleep_time)
	# # thread_count += threads_needed
	# # print('new thread count: %s' % thread_count)
	# print('\n\n' + ' '.join(cmd) + '\n\n')
	
	while not add_thread_count(thread_file, threads_needed):
		print('waiting for godot...')
		time.sleep(sleep_time)
	subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	while not sub_thread_count(thread_file, threads_needed):
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
				if global_max_threads - thread_count >= threads:
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
