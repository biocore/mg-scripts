#!/bin/bash
set -x

### 263 -- no job script found
### 265 -- data is already processing
### 266 -- bcl already converted
### 267 -- labadmin run directory does not exist at /opt/samplesheets/

source ~/.seqvars
source ~/.bashrc
source /etc/profile.d/pbs-maui.sh

seqpath="$1"

echo seqpath=$seqpath
curl=$(which curl)

fastq_output=$output_dir
### mac does not support -printf unless install findutils
### in which case you need to use gfind
function platform() {
  if [[ $(uname) == 'Darwin' ]]; then
    FIND=$(which gfind)
    echo $FIND
  elif [[ $(uname) == 'Linux' ]]; then
    FIND=$(which find)
    echo $FIND
  fi
}
myfind=$(platform myvar)

function prep_data_loc() {
### check for RTAComplete and samplesheet
### technically shouldn't get here unless RTAComplete exists but...

	if [[ -e $seqdir/RTAComplete.txt ]] && [[ ! -e $seqdir/*.csv ]]; then
		:
	else		
#		echo "RTAComplete present but no sample sheet for $seqdir" | /bin/mailx -s "imcomplete information for $seqdir" \
		#	-r ${recipient_name} <<EOF
		return -1
	fi

  if [[ -e $seqdir/alockfile ]]; then
	   echo "data already processing?"
	 return 265
  elif [[ -e $seqdir/processed ]]; then
	   echo "bcl conversion complete?"
     return 266
  elif [[ ! -e $seqdir/alockfile ]] || [[ ! -e $seqdir/processed ]]; then
    touch $seqdir/alockfile
    if [[ $? == 0 ]]; then
      pushd $seqdir/
      if [[ ! -d samplesheets ]]; then
				:
      fi

		  if [[ ! -d ${seqdir}/Data/Fastq ]]; then
      mkdir -p $seqdir/Data/Fastq
      echo making output directory $seqdir/Data/Fastq
    fi	
      fi
      popd
      return 0
    else
      echo failed to acquire lockfile
      return 1
    fi
  elif [[ -e $seqdir/processed ]]; then
    return 1
    :
  fi
}

function process_data() {
echo cwd == $(pwd)
for csvfile in $(ls ${seqdir}/*.csv | grep -v "sav.csv")
do
  labadmin_run=$(awk '/Description/{getline;print;}' $csvfile | cut -f3 -d",")
  if [[ $(ls /opt/samplesheets/$labadmin_run/ | wc -l) > "1" ]]; then
    multi_convert="TRUE"
  elif [[ $(ls /opt/samplesheets/$labadmin_run/ | wc -l) < "1" ]]; then
    echo "where are my samplesheets?"
    multi_convert="FALSE"
  fi

	mail_list=$(awk '/# Contact/{getline;print;}' $csvfile | cut -d',' -f 2-10 | tr -s ,)
	if [[ -z ${mail_list} ]]; then
		mail_some_people=$maillist
	else
		mail_some_people=${mail_list}
	fi


  rn1=$(expr $(awk '/Reads/{getline;print;}' $csvfile | cut -f1 -d",")) ### - 1)
  rn2=$(expr $(awk '/Reads/{getline;getline;print}' $csvfile | cut -f1 -d",")) ### - 1)
  echo read1 $rn1
  echo read2 $rn2
  if [[ -z $rn2 ]]; then
    ### single lane
    direction="1"
  else
    direction="2"
  fi
  ###rn1=$(expr $rn1 - 1 )
  job_read_val="Y$rn1"

### get qiita project identifier
  ###qiita_proj=$(awk '/Description/{getline;print;}' $csvfile | cut -f3 -d",")

### get comma count.  8 == 16s, 12bp
### comma count > 8 == metagenomic
  column_count=$(awk '/Sample_ID/{getline;print;}' $csvfile | tr -cd , | wc -c)
  echo column count $column_count
  if [[ $column_count -eq "8" ]]; then
    if [[ $(awk '/NNNNNNNNNNNN/{getline;print;}' $csvfile | cut -f6 -d",") -eq "NNNNNNNNNNNN" ]]; then
      echo "strip line from file with false barcode"
    fi
    cutval="7"
  elif [[ $column_count > "8" ]]; then
    cutval=""
  fi

  index_value_6=$(awk '/Sample_ID/{getline;print;}' $csvfile | cut -f6 -d",")
	index_value_8=$(awk '/Sample_ID/{getline;print;}' $csvfile | cut -f8 -d",")
  index_size=${#index_value}
  index_6_size=${#index_value_6}

  if [[ -z $index_value_6 ]]; then
    job_index_val="I12"
  elif [[ $index_value_6 -eq "NNNNNNNNNNNN" ]]; then
		#if [[ $index_value_8 -eq "NNNNNNNN" ]]; then
		#	sed -i s/NNNNNNNN//g $csvfile
		#fi
    sed -i.bak s/NNNNNNNNNNNN//g $csvfile
    job_index_val="I12"
  elif [[ $index_6_size -eq "8" ]]; then
    job_index_val="I8"
  elif [[ $index_6_size -eq "12" ]]; then
    job_index_val="I12"
  fi

  ### if index 7 == NNNNNNNNNNNN then must be metagenomic with multiple rows
  ### because of additional first column
	### should have label in sample sheet
  if [[ $index_value_7 -eq "NNNNNNNNNNNN" ]]; then
    sed -i.bak s/NNNNNNNNNNNN/g $csvfile
    job_index_val="I12"
  fi

  if [[ $direction -eq "1" ]]; then
    base_mask="--use-bases-mask Y$rn1,$job_index_val"
  elif [[ $direction -eq "2" ]]; then
    base_mask="--use-bases-mask Y$rn1,$job_index_val,$job_index_val,Y$rn1"
  fi

  ### get experiment name
  exp_name=$(awk -F',' '/Experiment Name/{print $2}' $csvfile)
	echo experiment name == $exp_name

	root_sequence=$(basename $dirpath)

	if [ -z $exp_name ]; then
		exp_name=$(basename $dirpath)
		echo experiment name 1 == $exp_name
	else
		exp_name=$exp_name
		echo experiment name 2 == $exp_name
	fi

  email_list=null
  job_o_out=localhost://${seqdir}
  job_e_out=localhost://${seqdir}

  bcl_template=${job_template}

  ### submit job
  if [[ ! -e $bcl_template ]]; then
    echo "can't locate job submit script" $bcl_template
    return 263
  else
		fastq_output="${output_dir}$(basename $dirpath)"
    #fastq_output="${output_dir}$(basename $dirpath)/Data/Fastq/ "
    ### this is trying to be created twice
    if [[ ! -d ${fastq_output} ]]; then
      mkdir -p $fastq_output
			echo creating output directory $fastq_output
    fi

    if [[ $? == 0 ]] || [[ -d $fastq_output ]]; then
			### changed output_dir=$fastq_output to $seqdir to write untrimmed files to private location
      jobnum=$(qsub -V -v seqdir="$seqdir",outputdir="$fastq_output",csvfile="$csvfile",job_o_out="$job_o_out",job_e_out="$job_e_out",base_mask="$base_mask" -N $exp_name ${bcl_template} )

			
			touch $fastq_output/$jobnum.txt
			echo start date == `date` > $fastq_output/$jobnum.txt
			echo experiment name == $exp_name >> $fastq_output/$jobnum.txt
			echo job_id == $jobnum >> $fastq_output/$jobnum.txt
			echo submit args == qsub -v seqdir="$seqdir",outputdir="$fastq_output",csvfile="$csvfile",job_o_out="$job_o_out",job_e_out="$job_e_out" -N $exp_name ${bcl_template} >> $fastq_output/$jobnum.txt 
			### removed for now.  
			### this should check if assay is metagenomics
			### if not skip
			### if so, filter sequences and move to final output
			#filter_job_num=$(qsub -V -v outputdir="$fastq_output" -W depend=afterok:$jobnum atropos_filter_project.sh)

      while [[ $(qstat $jobnum) ]]; do
	      if [[ $? -eq "153" ]]; then
	          break
	      else
        sleep 60
	      fi
      done

			echo end date == `date` >> $fastq_output/$jobnum.txt
					
      ### in case run spans a day change
      lock_file_age=$((($(date +%s) - $(date +%s -r $seqdir/alockfile)) / 86400 ))
	### shouldn't need this now that we are dumb parallel
			new_age=$(expr $lock_file_age + 1)
      if [[ $(tracejob -n${new_age} $jobnum | grep "Exit_status" | awk '{print $4}' | cut -f 2 -d"=") == 0 ]]; then
				echo exit status == 0 for $jobnum >> $fastq_output/$jobnum.txt

/bin/mailx -s "processing complete for project $exp_name " \
					-a $csvfile -r jdereus@ucsd.edu $mail_some_people <<EOF	

Job $jobnum complete on `date`
Run processed for experiment $exp_name
Data located on barnacle at $fastq_output
EOF

			/bin/cp -f $csvfile $fastq_output
			else 
				echo exit status == $($tracejob -n${new_age} $jobnum | grep "Exit_status" | awk '{print $4}' | cut -f 2 -d"=") >> $fastq_output/$jobnum.txt	
      fi
      :
    fi
  fi
  :
done

touch $seqdir/processed && sleep 2
rm $seqdir/alockfile
}

function human_filter() {
	pushd $output_dir

	for dir in `ls -d */ | egrep -v 'Reports|Stats'`; 
		do
			sh ~/atropos_filter.sh ${output_dir}/${dir}
		done

	popd
}

function concat_fastq() {
	pushd $output_dir
	if [[ ! -d ${output_dir}/concat_files ]]; then
		mkdir ${output_dir}/concat_files}
	fi

		for lane in L001; do
			for direction in R1 R2; do
				for file in `ls *.fastq.gz | grep $lane`; do
					unset filename1
					unset filename2
					filename1=$(basename $file)
					filename2=$(echo "$filename1" | sed -e 's/L001/L002/g')
					filename1_short=$(basename "$filename1" .fastq.gz)

					cat ${filename1} ${filename2} > concat_files/${filename1_short}.concat.fastq.gz
				done
			done
		done
	popd
}



FIND=$(platform FIND)
  ### check filesystem for new sequence directory
  ### if exists, check for RTAComplete.txt files
pushd $seqpath/
path_count=(echo ${#raw_sequencing_loc[@]})
declare -a newdirs=$(echo "("; $FIND $seqpath/ -maxdepth 2 ! -path $seqpath/ -type d -mtime -1 ; echo ")")
count=${#newdirs[@]}
popd
labadmin_run=""


for dirpath in "${newdirs[@]}"
  do
    echo dirpath==$dirpath
		echo myfind==$myfind
		echo dirpath==$dirpath
    seqloc=$($myfind $dirpath -maxdepth 1 ! -path $dirpath -type f -name "RTAComplete.txt")

		#if [[ -e $dirpath/RTAComplete.txt ]]; then
    seqdir=$dirpath    
		echo seqdir==$seqdir		
    if [[ $seqloc ]]; then
      prep_data_loc "$seqdir" "$labadmin_run"

      if [[ $? == 0 ]]; then
				
        process_data "$seqdir" "$dirpath" "$curl"
      else
        :
      fi
	#			human_filter "$seqdir" "$output_dir"
    fi

		### end checking for RTAComplete.txt
		#fi
  done
