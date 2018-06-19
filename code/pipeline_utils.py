import os
import time
import subprocess

def confirm_path(file):
	if not os.path.exists(os.path.dirname(file)):
		try:
			os.makedirs(os.path.dirname(file))
		except OSError as exc: # Guard against race condition
			if exc.errno != errno.EEXIST:
				raise

def command_call(cmd, threads_needed=1, sleep_time=1):
	# global thread_count, global_max_threads
	# while global_max_threads - thread_count < threads_needed:
	# 	time.sleep(sleep_time)
	# thread_count += threads_needed
	# print('\n\n\n\n' + ' '.join(cmd) + '\n\n\n\n')
	subprocess.call(cmd, stdout=subprocess.PIPE)
	# thread_count -= threads_needed
