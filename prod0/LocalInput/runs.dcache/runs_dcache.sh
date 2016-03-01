#!/bin/bash


for inds in `grep data runstofix.txt`
#for inds in `grep data onerun.txt`
#for inds in `grep data runs2012_d3pd.txt`
do 
    runs=`dq2-ls -f -p -L IN2P3-CC_PHYS-HIGGS $inds | awk '/ccsrm.in2p3.fr/ { print $0 }' | sort`
    runnumber=`dq2-ls $inds | awk -F"data12_8TeV.00" '{print $2}'| awk -F".physics_Egamma" '{print $1}'`
    echo $runnumber

    files_list=""
    let Nfiles=0
    for run in $runs 
    do 
	echo $run
	if [ $Nfiles -eq 0 ]; then
	    files_list="$run"
	else
	    files_list="${files_list},$run"
	fi
	  let "Nfiles += 1" 
    done
    echo $files_list > input_data_run${runnumber}_new.txt
    sed -i 's/srm:\/\/ccsrm.in2p3.fr/dcap:\/\/ccdcapatlas.in2p3.fr:22125/g' input_data_run${runnumber}_new.txt
    echo "**** DIFF OLD NEW ***"
    diff input_data_run${runnumber}.txt input_data_run${runnumber}_new.txt
done



# dcap_path="dcap://ccdcapatlas.in2p3.fr:22125/pnfs/in2p3.fr/data/atlas/atlasgroupdisk/phys-higgs/data12_8TeV/NTUP_PHOTON/"
# for inds in `grep data runs2012_d3pd.txt`
# for inds in `grep data onerun.txt`
# do 

#     dsdir=`dq2-ls -r $inds | grep tid | awk -F":" '{print $1}' ` 
#     tag=`dq2-ls $inds | awk -F"NTUP_PHOTON." '{print $2}'| awk -F"/" '{print $1}'`
#     run=`dq2-ls $inds | awk -F"data12_8TeV.00" '{print $2}'| awk -F".physics_Egamma" '{print $1}'`
#     echo "BASIC INFO ON THE RUN"
#     echo $inds
#     echo $dsdir
#     echo $tag
#     echo $run

#     files_list=""
#     let Nfiles=0
#  ##---> loop over the rootfiles in the run
#     for file in `dq2-ls -f "*root*" $inds | grep root | awk -F" " '{print $3}'`
#     do
# 	echo $dsdir"/"$file
# 	file_withpath=$dcap_path$tag"/"$dsdir"/"$file
# 	if [ $Nfiles -eq 0 ]; then
# 	    files_list="$file_withpath"
# 	else
# 	    files_list="${files_list},$file_withpath"
# 	fi
# 	  let "Nfiles += 1" 
#     done
#     echo $files_list > input_data_run${run}.txt
# done



