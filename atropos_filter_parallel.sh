#!/bin/bash
#set -x

dir=$1
trim_file=$2
NPROCS=$3

source /home/jede9131/miniconda3/bin/activate test_env_2

module load bedtools_2.26.0 samtools_1.3.1 bowtie2_bowtie2-2.2.3

bowtie=$(which bowtie2)
samtools=$(which samtools)
bedtools=$(which bedtools)
suffix=R*.fastq*
NPROCS=16 

filter_db="/databases/bowtie/Human_phiX174/Human_phix174"
tar="/bin/tar"

#atropos_param="-a ATCTCGTATGCCGTCTTCTGCTTG -A GTGTAGATCTCGGTGGTCGCCGTATCATT -q 15 --minimum-length 100 --pair-filter any"
atropos_param="-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT -q 15 --minimum-length 100 --pair-filter any"

atropos_qc_output=$dir/atropos_qc

pushd $dir

  if [[ ! -d $atropos_qc_output ]]; then
    mkdir -p $atropos_qc_output
    mkdir ${atropos_qc_output}/atropos_logs
		mkdir ${atropos_qc_output}/atropos_fastqc
	fi

#for file in `find $dir -maxdepth 1 -type f -name "*.fastq.gz" | grep R1`;
#for file in `cat ${atropos_qc_output}/${trim_file}`;
for file in `cat ${dir}/${trim_file} | grep "_R1_"`; 
do
  parent_dir=$(dirname $file)
  project_dir=$(echo $parent_dir | cut -f 2 -d"/")
      ### actual filename minus the preceding path
  filename1=$(basename "$file") ### .fastq.gz)
      ### strip fastq.gz for better file naming regarding output
  filename1_short=$(basename "$filename1" .fastq.gz)
      ### replace R1 with R2 to be able to pass both forward and reverse reads to bowtie
  filename2=$(echo "$filename1" | sed -e 's/_R1_/_R2_/g')
  filename2_short=$(basename "$filename2" .fastq.gz)

  atropos --threads 16 ${atropos_param} --report-file ${atropos_qc_output}/atropos_logs/${filename1_short}.log --report-formats txt -o ${atropos_qc_output}/${filename1_short}.fastq -p ${atropos_qc_output}/${filename2_short}.fastq -pe1 $filename1 -pe2 $filename2

done

#for file in `find ${atropos_qc_output} -maxdepth 1 -type f -name "*.fastq" | grep _R1_`;
#for file in `cat ${atropos_qc_output}/${trim_file}`;
for file in `cat ${dir}/${trim_file}`; 
do
	#final_output=./filtered_sequences
	final_output=$dir/filtered_sequences
	if [[ ! -d ${final_output} ]]; then
		mkdir ${final_output}
		mkdir ${final_output}/trim_logs
	fi

	if [[ ! -d ${atropos_qc_output}/processed_qc ]]; then
		mkdir ${atropos_qc_output}/processed_qc
	fi

	parent_dir=$(dirname $file)
  project_dir=$(echo $parent_dir | cut -f 2 -d"/")
      ### actual filename minus the preceding path
  filename1=$(basename "$file" .gz) ### .fastq.gz)
      ### strip fastq.gz for better file naming regarding output
  filename1_short=$(basename "$filename1" .fastq)
	filename1_trim=$(echo "filename1_short" | sed -e 's/atropos//g')
      ### replace R1 with R2 to be able to pass both forward and reverse reads to bowtie
  filename2=$(echo "$filename1" | sed -e 's/_R1_/_R2_/g')
  filename2_short=$(basename "$filename2" .fastq)
	filename2_trim=$(echo "filename2_short" | sed -e 's/atropos//g')

        $bowtie -p 16 -x ${filter_db} --very-sensitive -1 ${atropos_qc_output}/$filename1 -2 ${atropos_qc_output}/$filename2 | $samtools view -f 12 -F 256 | $samtools sort -@ 16 -n | $samtools view -bS | $bedtools bamtofastq -i - -fq $final_output/${filename1_short}.trimmed.fastq -fq2 $final_output/${filename2_short}.trimmed.fastq &> $final_output/trim_logs/${filename1_short}.log

	gzip -f ${final_output}/${filename1_short}.trimmed.fastq
	gzip -f ${final_output}/${filename2_short}.trimmed.fastq

	mv ${atropos_qc_output}/${filename1} ${atropos_qc_output}/processed_qc
  mv ${atropos_qc_output}/${filename2} ${atropos_qc_output}/processed_qc

done

for file in `cat $trim_file`; do
	file=$(basename $file .gz)
	echo "stripped filename == " $file
	#rm ${atropos_qc_output}/processed_qc/${file}

done

### option for removing index files 

#	if [[ ! -d ${final_output}/index_files ]]; then
#		mkdir ${final_output}/index_files}
#	fi
#	mv $final_output/*I[12]*.filtered.fastq.gz

popd
