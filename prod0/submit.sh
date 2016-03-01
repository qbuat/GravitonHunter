#!/bin/bash

### Merge stdout et stderr in a single file
#$ -j y

### use sps
#$ -l sps=1

### use xrootd
#$ -l xrootd=1


run_dir=/afs/in2p3.fr/home/q/qbuat/mygroup/GravitonHunter_new_geo/prod1
log_dir=${run_dir}/log/

ISO=${1}
INPUT_FILE=${2}
OUTPUT_FILE=${3}

echo "Hello World!"

echo "My working directory is"
pwd

echo "On the host"
hostname

echo "Athena is setup ?"
echo $TestArea

echo "library path is"
echo $LD_LIBRARY_PATH

cd ${run_dir}
echo "My working directory is now"
pwd

echo "RUN !"
echo ${ISO}
echo ${INPUT_FILE}
echo ${OUTPUT_FILE}
echo './create_skimmed_files_x "${ISO}" 1 "" "${INPUT_FILE}" "${OUTPUT_FILE}"'
./create_skimmed_files_x ${ISO} 1 "" ${INPUT_FILE} ${OUTPUT_FILE}
