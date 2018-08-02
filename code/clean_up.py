import sys
import os
import shutil
import subprocess

def retrieve_filenames(directory):
	f = []
	for (dirpath, dirnames, filenames) in os.walk(directory):
		for filename in filenames:
			f.append(os.path.join(dirpath,filename))
		# break
	return f

def confirm_path(file):
	if not os.path.exists(os.path.dirname(file)):
		try:
			os.makedirs(os.path.dirname(file))
		except OSError as exc: # Guard against race condition
			if exc.errno != errno.EEXIST:
				raise

if __name__ == '__main__':
	new_dir = sys.argv[1]
	old_dir = sys.argv[2]
	# print(new_dir)
	# print(old_dir)
	new_files = retrieve_filenames(new_dir)
	old_files = retrieve_filenames(old_dir)
	# old_files_relative = []
	# for file in old_files:
	# 	relative_path = file.split(old_dir)[-1]
	# 	if relative_path.startswith('/'):
	# 			relative_path = relative_path[1:]
	# 	old_files_relative.append(relative_path)
	old_files = [file.split(old_dir)[-1] for file in old_files]

	# print(new_files)
	# print(old_files)
	# exit()
	for file in new_files:
		if file.split(new_dir)[-1] not in old_files:
			relative_path = file.split(new_dir)[-1]
			if relative_path.startswith('/'):
				relative_path = relative_path[1:]
			
			confirm_path(os.path.join(old_dir, relative_path))
			print(os.path.join(old_dir, relative_path))

			subprocess.call('mv %s %s' % (file, os.path.join(old_dir, relative_path))
			# os.rename(file, os.path.join(old_dir, relative_path))

	shutil.rmtree(new_dir)
