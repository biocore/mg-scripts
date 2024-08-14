#!/bin/bash -l
#SBATCH -J cleanup             # cleanup
#SBATCH --time 24:00:00  # 24:00:00
#SBATCH --mem 8G        # 8G
#SBATCH -N 1           # 1
#SBATCH -c 1       # 1
#SBATCH -p qiita           # qiita

# for now comment these out as qiita is responsible for notifying users.
###SBATCH --mail-user=qiita.help@gmail.com
###SBATCH --mail-type=BEGIN,FAIL

# for now these can be left hard-coded.
#SBATCH --output %x-%A.out
#SBATCH --error %x-%A.err

if [[ -z "${OUTPUT}" ]]; then
    echo "OUTPUT is not specified"
    exit 1
fi

# remove unused large outputs
rm -rf ${OUTPUT}/biosample_format ${OUTPUT}/1_demult ${OUTPUT}/Full