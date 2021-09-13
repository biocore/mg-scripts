import os
from os.path import join, exists, dirname, isdir, isfile, basename
from os import makedirs, listdir
import shutil
from sequence_processing_pipeline.Job import Job
from time import sleep
from zipfile import ZipFile


#source ~/.seqvars
#module load fastp_0.20.1 samtools_1.12 minimap2_2.18


class TrimFile():
	def __init__(self, trim_file, current_directory):
		self.path = trim_file
		self.list = []
		# current_directory should be something like:
		# /pscratch/seq_test/210518_A00953_0305_AHCJT7DSX2/Data/Fastq/THDMI_US_10317
		if exists(self.path):
			with open(self.path, 'r') as f:
				l = f.readlines()
				l = [x.strip() for x in l]
				self.list += l
		else:
			# don't create a trim_file just yet. keep it in memory.
			# with open(self.path, 'w') as f:
			for item in listdir(current_directory):
				if isfile(join(current_directory, item)):
					if item.endswith('.fastq.gz'):
						if '_R1_' in item:
							self.list.append(item)
						# f.write('%s\n' % item)


class Encapsulation(Job):
	def __init__(self, dir, file_dir, output_dir, base_seq_dir, project, final_output_dir, trim_file, current_directory):
		super().__init__()
		self.trim_file_obj = TrimFile(trim_file, current_directory)
		self.output_dir = output_dir

		# this should be written out to $final_output/${project}/empty_file_list.txt
		# by either this method or the caller().
		self.empty_file_list_txt = []

		# we have to create a lot of new directories and the conditions
		# involved for creating each one can be involved. It's better to
		# simply create all of the required directories at the beginning and
		# walk the tree at the end of processing to remove any empty trees.
		path_list = []
		self.project = project
		self.final_output_dir = final_output_dir

		dir = join(dir, file_dir)
		self.dir = dir
		path_list.append(join(dir, project, 'html'))
		path_list.append(join(dir, project, 'json'))

		final_output = join(output_dir, base_seq_dir)
		path_list.append(final_output)
		path_list.append(join(final_output, 'trim_logs'))
		path_list.append(join(final_output, 'zero_files'))

		self.fastp_qc_output = join(dir, project, 'fastp_qc')

		path_list.append(self.fastp_qc_output)
		path_list.append(join(self.fastp_qc_output, 'fastp_logs'))
		path_list.append(join(self.fastp_qc_output, 'fastp_fastqc'))

		path_list.append(join(final_output_dir, project, 'amplicon'))

		final_output_project = join(final_output, project)
		path_list.append(final_output_project)
		path_list.append(join(final_output_project, 'filtered_sequences'))
		path_list.append(join(final_output_project, 'trimmed_sequences'))
		path_list.append(join(final_output_project, 'html'))
		path_list.append(join(final_output_project, 'json'))
		path_list.append(join(final_output_project, 'trim_logs'))
		path_list.append(join(final_output_project, 'zero_files'))
		path_list.append(join(final_output_project, 'trimmed_sequences'))

		path_list.append(self.fastp_qc_output)
		path_list.append(join(self.fastp_qc_output, 'fastp_logs'))
		path_list.append(join(self.fastp_qc_output, 'fastp_fastqc'))
		path_list.append(join(self.fastp_qc_output, 'trimmed_sequences'))
		path_list.append(join(self.fastp_qc_output, 'processed_qc'))

		# this appears to be in addition to the identically named directories
		# created above.
		final_output = join(dir, project, 'filtered_sequences')
		path_list.append(final_output)
		path_list.append(join(final_output, 'trim_logs'))
		path_list.append(join(final_output, 'zero_files'))

		for path in path_list:
			makedirs(path, exist_ok=True)

	def zip_files(self, file_paths, zip_path):
		with ZipFile(zip_path, 'w') as zip:
			for file in file_paths:
				zip.write(file)

	def copy_directory_contents(self, source_directory, destination_directory, ends_with=None, contains=None, starts_with=None, delete_after=False):
		files = listdir(source_directory)
		# convert from bytes output of listdir to string
		files = [str(x) for x in files]
		for some_file in files:
			sentinel = False
			if starts_with:
				if some_file.startswith(starts_with):
					sentinel = True
			elif ends_with:
				if some_file.endswith(ends_with):
					sentinel = True
			elif contains:
				if contains in some_file:
					sentinel = True

			if sentinel:
				if delete_after:
					shutil.move(join(source_directory, some_file), destination_directory)
				else:
					shutil.copyfile(join(source_directory, some_file), destination_directory)

	def possible_amplicon(self, adapter_a, adapter_A, nprocs, chemistry, final_output, qiita_proj, copy_file):
		for file in self.trim_file_obj.list:
			parent_dir = dirname(file)
			tmp = parent_dir.split('/')
			# if tmp begins with a leading empty string '', this will remove it.
			tmp = [x for x in tmp if x]
			project_dir = tmp[0]
			filename1 = basename(file)
			filename1_short = filename1.replace('.fastq.gz', '')
			filename2 = filename1.replace('_R1_00', '_R2_00')
			filename2_short = filename2.replace('.fastq.gz', '')
			filename_index1 = filename1.replace('_R1_', '_I1_')
			filename_index2 = filename2.replace('_R1_', '_I1_')

			cmd_list = ['/full/path/to/fastp']

			if adapter_a == None or adapter_A == None:
				pass
			else:
				cmd_list.append('--adapter_sequence')
				cmd_list.append(adapter_a)
				cmd_list.append('--adapter_sequence_r2')
				cmd_list.append(adapter_A)

			partial_path = join(final_output, self.project)
			f1 = join(partial_path, 'trimmed_sequences', filename1_short + '.fastp.fastq.gz')
			f2 = join(partial_path, 'trimmed_sequences', filename2_short + '.fastp.fastq.gz')

			cmd_list += ['-l', '100',
					'-i', filename1,
					'-I', filename2,
					'-w', nprocs,
					'-j', join(partial_path, 'json', filename1_short + '.json'),
					'-h', join(partial_path, 'html', filename1_short + '.html'),
					'-o', f1,
					'-O', f2,
					'-R', filename1_short + '_report'
					]

			out, err, rc = self._system_call(cmd_list)

			size1 = os.stat(f1).st_size
			size2 = os.stat(f2).st_size

			### amplicon/16s data. include index files
			if chemistry == 'Amplicon':
				dst = join(partial_path, 'trimmed_sequences')
				# TODO: We're going to need to replace self.dir with
				#  the directory that contains these files.
				self.copy_directory_contents(self.dir, dst, contains=filename_index1)
				self.copy_directory_contents(self.dir, dst, contains=filename_index2)

			if size1 <= 500 or size2 <= 500:
				src = join(partial_path, 'trimmed_sequences')
				dst = join(partial_path, 'zero_files')
				self.copy_directory_contents(self.dir, dst, contains=filename1_short, delete_after=True)
				self.copy_directory_contents(self.dir, dst, contains=filename2_short, delete_after=True)

				self.empty_file_list_txt.append(join(partial_path, 'trimmed_sequences', filename1_short + '.fastp.fastq.gz'))
				self.empty_file_list_txt.append(size1)
				self.empty_file_list_txt.append(join(partial_path, 'trimmed_sequences', filename2_short + '.fastp.fastq.gz'))
				self.empty_file_list_txt.append(size2)

		for file_transfer in self.trim_file_obj.list:
			transfer_1 = 1  # $(echo "$file_transfer" | cut -f1 -d".")
			transfer_2 = 1  # $(echo "$transfer_1" | sed -e 's/_R1_00/_R2_00/g')
			transfer_3 = 1  # $(echo "$transfer_1" | sed -e 's/_R1_00/_I1_00/g')
			transfer_4 = 1  # $(echo "$transfer_1" | sed -e 's/_R1_00/_I2_00/g')
			if chemistry == 'Amplicon':
				pass
				#rsync -avp --progress ${final_output}/${project}/*.fastp.fastq.gz ${final_output_dir}/${base_seq_dir}/${project}
			else:
				pass
				#rsync -avp --progress ${final_output}/${project}/${transfer_1}.fastp.fastq.gz ${final_output}/${project}/${transfer_2}.fastp.fastq.gz ${final_output_dir}/${base_seq_dir}/${project}
			if qiita_proj == 'NA':
				if chemistry == 'Amplicon' and copy_file == 'TRUE':
					pass
					#sudo -u qiita rsync -avp --progress --chown 5500:5500 ${final_output}/${project}/trimmed_sequences/*_R[12]*.fastq.gz ${final_output}/${project}/trimmed_sequences/*_I1*.fastq.gz /qmounts/qiita_data/uploads/${qiita_proj}/

				elif copy_file == 'TRUE':
					pass
					#sudo -u qiita rsync -avp --chown 5500:5500 ${final_output}/${project}/trimmed_sequences/${transfer_1}.fastp.fastq.gz ${final_output}/${project}/trimmed_sequences/${transfer_2}.fastp.fastq.gz /qmounts/qiita_data/uploads/${qiita_proj}/

			sleep(1)
			# rsync -avp --progress ${final_output}/${project}/ ${final_output_dir}/${base_seq_dir}

	def possible_amplicon2(self, adapter_a, adapter_A, final_output, copy_file, qiita_proj):
		partial_path = join(final_output, self.project)
		for file in self.trim_file_obj.list:
			parent_dir = dirname(file)
			tmp = parent_dir.split('/')
			# if tmp begins with a leading empty string '', this will remove it.
			tmp = [x for x in tmp if x]
			project_dir = tmp[0]
			filename1 = basename(file)
			filename1_short = filename1.replace('.fastq.gz', '')
			filename2 = filename1.replace('_R1_00', '_R2_00')
			filename2_short = filename2.replace('.fastq.gz', '')
			if adapter_a == 'NA' or adapter_A == 'NA':
				pass
			# $fastp -l 100 -i $filename1 -I $filename2 -w $NPROCS --stdout -j ${dir}/${project}/json/${filename1_short}.json -h ${dir}/${project}/html/${filename1_short}.html | $minimap2 -ax sr -t $NPROCS $db - -a | $samtools fastq -@ $NPROCS -f 12 -F 256 -1 $final_output/${project}/filtered_sequences/${filename1_short}.trimmed.fastq.gz -2 $final_output/${project}/filtered_sequences/${filename2_short}.trimmed.fastq.gz
			else:
				pass
			# $fastp --adapter_sequence ${adapter_a} --adapter_sequence_r2 ${adapter_A} -l 100 -i $filename1 -I $filename2 -w $NPROCS --stdout -j ${dir}/${project}/json/${filename1_short}.json -h ${dir}/${project}/html/${filename1_short}.html | $minimap2 -ax sr -t $NPROCS $db - -a | $samtools fastq -@ $NPROCS -f 12 -F 256 -1 $final_output/${project}/filtered_sequences/${filename1_short}.trimmed.fastq.gz -2 $final_output/${project}/filtered_sequences/${filename2_short}.trimmed.fastq.gz

			f1 = join(partial_path, 'filtered_sequences', filename1_short + '.trimmed.fastq.gz')
			f2 = join(partial_path, 'filtered_sequences', filename2_short + '.trimmed.fastq.gz')

			size1 = os.stat(f1).st_size
			size2 = os.stat(f2).st_size

			if size1 <= 500 or size2 <= 500:
				src = join(partial_path, 'trimmed_sequences')
				dst = join(partial_path, 'zero_files')

				# mv $final_output/${project}/filtered_sequences/${filename1_short}* $final_output/${project}/filtered_sequences/${filename2_short}* ${final_output}/${project}/zero_files
				self.copy_directory_contents(self.dir, dst, contains=filename1_short, delete_after=True)

				self.empty_file_list_txt.append(
					join(partial_path, 'filtered_sequences$', filename1_short + '.trimmed.fastq.gz'))
				self.empty_file_list_txt.append(size1)
				self.empty_file_list_txt.append(
					join(partial_path, 'filtered_sequences$', filename2_short + '.trimmed.fastq.gz'))
				self.empty_file_list_txt.append(size2)

		### always executes
		if copy_file:
			for file_transfer in self.trim_file_obj.list:
				# transfer_1=$(echo "$file_transfer" | cut -f1 -d".")
				# transfer_2=$(echo "$transfer_1" | sed -e 's/_R1_00/_R2_00/g')
				# rsync -avp --progress ${final_output}/${project}/filtered_sequences/${transfer_1}.trimmed.fastq.gz ${final_output}/${project}/filtered_sequences/${transfer_2}.trimmed.fastq.gz ${final_output_dir}/${base_seq_dir}/${project}/filtered_sequences/
				if qiita_proj != 'NA' and copy_file == 'TRUE':
					pass

	# sudo -u qiita rsync -avp --chown 5500:5500 ${final_output}/${project}/filtered_sequences/${transfer_1}.trimmed.fastq.gz ${final_output}/${project}/filtered_sequences/${transfer_2}.trimmed.fastq.gz /qmounts/qiita_data/uploads/${qiita_proj}/

	# rsync -avp --progress ${final_output} ${final_output_dir}/

	def main(self, rundir, filter_db, a_trim, h_filter, chemistry, qiita_proj, copy_file, adapter_a, adapter_A, nprocs, final_output):
		if not isinstance(a_trim, bool):
			raise ValueError("a_trim must be a boolean value.")

		if not isinstance(h_filter, bool):
			raise ValueError("h_filter must be a boolean value.")

		if exists(join(rundir, 'run_config.h')):
			pass
			# TODO: source the file - we'll need to find an alternative method
			#  to take in the settings

		filter_db="/databases/bowtie/Human/Human"

		# move sample sheets into place
		# pull everything from original location to top level output
		self.copy_directory_contents(self.dir, self.output_dir, '.csv')

		### amplicon MiSeq sequencing
		if not a_trim and not h_filter:
			### just copy data to qiita study.  16s data.
			if chemistry == 'Amplicon':
				src = join(self.dir, self.project)
				dst = join(self.final_output_dir, self.project, 'amplicon')
				self.copy_directory_contents(src, dst, '.fastq.gz')
			if qiita_proj:
				if qiita_proj != 'NA' and copy_file != 'FALSE':
					src = join(self.dir, self.project)
					dst = join('qmounts', 'qiita_data', 'uploads', qiita_proj)
					# these probably need to change ownership to qiita and group to the proper group as well.
					self.copy_directory_contents(src, dst, '.fastq.gz')

		if a_trim and not h_filter:
			self.possible_amplicon(adapter_a, adapter_A, nprocs, chemistry, final_output, qiita_proj, copy_file)
		elif a_trim and h_filter:
			self.possible_amplicon2(adapter_a, adapter_A, final_output, copy_file, qiita_proj)

		if not h_filter:
			#if the files failed to transfer then exit with -1
			pass
		else:
			exit(0)

		do_bowtie = True

		if do_bowtie:
			for some_file in self.trim_file_obj.list:
				final_output = join(self.dir, self.project, 'filtered_sequences')
				parent_dir = dirname(some_file)
				tmp = parent_dir.split('/')
				# if tmp begins with a leading empty string '', this will remove it.
				tmp = [x for x in tmp if x]
				project_dir = tmp[0]
				filename1 = basename(some_file)
				filename1_short = filename1 + '.fastq.gz'
				filename2 =  filename1.replace('_R1_', '_R2_')
				filename2_short = filename2 + '.fastq.gz'

				from subprocess import Popen, PIPE

				p1 = join(self.fastp_qc_output, filename1)
				p2 = join(self.fastp_qc_output, filename2)

				bowtie_log = join(final_output, 'trim_logs', filename1_short + '.log')
				with open(bowtie_log, 'a') as f:
					bowtie = Popen(['bowtime','-p', nprocs, '-x', filter_db, '--very-sensitive', '-1', p1, '-2', p2], stdout=PIPE, stderr=f)
					samtools1 = Popen(['samtools', 'view', '-f', '12', '-F', '256'], stdout=PIPE, stderr=PIPE, stdin=bowtie.stdout)
					samtools2 = Popen(['samtools', 'sort', '-@', '16', '-n'], stdout=PIPE, stderr=PIPE, stdin=samtools1.stdout)
					samtools3 = Popen(['samtools', 'view', '-bS'], stdout=PIPE, stderr=PIPE, stdin=samtools2.stdout)
					bedtools = Popen(['bedtools', 'bamtofastq', '-i', '-', '-fq', join(final_output, filename1_short + '.trimmed.fastq'), '-fq2',join(final_output, filename1_short + '.trimmed.fastq')], stdout=f, stderr=PIPE, stdin=samtools3.stdout)

					stdout, stderr = bedtools.communicate()
					return_code = bedtools.returncode

					file1 = join(final_output, filename1_short + '.trimmed.fastq')
					self.zip_files(file1, file1 + '.zip')
					file2 = join(final_output, filename2_short + '.trimmed.fastq')
					self.zip_files(file2, file2 + '.zip')

					# TODO: Find an alternate method to check and see if the two
					#  files gzipped above should be added to zero_files.txt and
				    #  the gzips deleted.
					for some_file in self.trim_file_obj.list:
						tmp = some_file.split('_')
						# extract 'X00178928_S1_L003' from
						# 'X00178928_S1_L003_R1_001.fastq.gz' ...
						tmp = '_'.join(tmp[0:3])

						# this code may be vestigial. the contents of both
						# conditionals were commented out.
						if qiita_proj != 'NA':
							# echo copying files ${final_output}/${file}_R1_*.trimmed.fastq.gz ${final_output}/${file}_R2_*.trimmed.fastq.gz
							pass
						else:
							# echo syncing file ${final_output}/${file}_R1_*.trimmed.fastq.gz ${final_output}/${file}_R2_*.trimmed.fastq.gz
							pass
