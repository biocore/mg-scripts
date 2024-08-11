#!/bin/bash -l
#SBATCH -J cleanup
#SBATCH --time 24:00:00
#SBATCH --mem 8gb
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mail-user=qiita.help@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --output %x-%A.out
#SBATCH --error %x-%A.err
#SBATCH -p qiita

if [[ -z "${OUTPUT}" ]]; then
    echo "OUTPUT is not specified"
    exit 1
fi

# remove unused large outputs
rm -fr ${OUTPUT}/biosample_format ${OUTPUT}/1_demult ${OUTPUT}/Full
