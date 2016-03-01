#!/bin/bash

#---------------------------------------------------------------------------------------
##################   RUN CONFIGURATION
#---------------------------------------------------------------------------------------
debug_mode=0 #--> Set to 1 to print the command that will be run without running it

data_debug=0
data=0                  #-->Set to 1 to create the skimmed data files
mc_rs=0                 #-->Set to 1 to create the skimmed RS files
mc_gg=1                 #-->Set to 1 to create the skimmed GG files
mc_gj=0                 #-->Set to 1 to create the skimmed GJ files
iso_type="TOPO_PTDEP"   #-->Choice of the isolation criteria ("TOPO","TOPO_PTDEP","CONE")
#prod_version="prod_12_07_13" #-->NAME OF THE OUTPUT DIRECTORY (MUST BE CREATED BEFORE RUNNING)
prod_version="prod_07_10_13" #-->NAME OF THE OUTPUT DIRECTORY (MUST BE CREATED BEFORE RUNNING)
#outputpath="/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis8TeV/" #--> Path to find the output directory
#outputpath="OutputDir/" #--> Path to find the output directory
outputpath="/sps/atlas/j/jbrown/OUTDS/GravitonAnalysis8TeV/" #--> Path to find the output directory

outputdir="${outputpath}${prod_version}/"
inputdir="LocalInput/"
datafiledir="${inputdir}runs/"
debugfiledir="${inputdir}debug/"


#---------------------------------------------------------------------------------------
##################   DATA - DATA - DATA
#---------------------------------------------------------------------------------------
if [ "$data" == "1" ] ; then 
    let Nrun=0
    for file in `ls $datafiledir* | grep input_data_run | awk -F"runs/" '{print $2}'`
    do
	array[$Nrun]=$file
	let "Nrun += 1"
    done
    echo $Nrun
    # for i in `seq 1 ${Nrun}`
    for i in `seq 1 ${Nrun}`
    do 
	input=${array[i-1]}
	runnumber=`echo $input | awk -F"_" '{print $3}' | awk -F"." '{print $1}'`
	output="datafile_${iso_type}_${runnumber}.root"
	
	if [ "$debug_mode" == "1" ] ; then 
	    echo "./create_skimmed_files_x ${iso_type} 1 '' $datafiledir$input $outputdir$output"
	else
	    ./create_skimmed_files_x $iso_type 1 "" $datafiledir$input $outputdir$output
	fi
    done
    echo "DATA SKIMMED FILES HAVE BEEN CREATED"
fi



#---------------------------------------------------------------------------------------
##################   DATA DEBUG STREAM - DATA DEBUG STREAM - DATA DEBUG STREAM
#---------------------------------------------------------------------------------------
if [ "$data_debug" == "1" ] ; then 
    let Nrun=0
    for file in `ls $debugfiledir* | grep input_debug_run | awk -F"debug/" '{print $2}'`
    do
	array[$Nrun]=$file
	let "Nrun += 1"
    done
    echo $Nrun
    # for i in `seq 1 ${Nrun}`
    for i in `seq 1 ${Nrun}`
    do 
	input=${array[i-1]}
	runnumber=`echo $input | awk -F"_" '{print $3}' | awk -F"." '{print $1}'`
	output="debugfile_${iso_type}_${runnumber}.root"
	
	if [ "$debug_mode" == "1" ] ; then 
	    echo "./create_skimmed_files_x ${iso_type} 1 '' $debugfiledir$input $outputdir$output"
	else
	    ./create_skimmed_files_x $iso_type 1 "" $debugfiledir$input $outputdir$output
	fi
    done
    echo "DATA DEBUG STREAM SKIMMED FILES HAVE BEEN CREATED"
fi







#---------------------------------------------------------------------------------------
##################   PYTHIA RS FILES - PYTHIA RS FILES
#---------------------------------------------------------------------------------------
if [ "$mc_rs" == "1" ] ; then 
    for file in `ls $inputdir* | grep input_mc12_Ggg | awk -F"${inputdir}" '{print $2}'`
    do
	output_label=`echo $file | awk -F"input" '{print $2}' | awk -F".txt" '{print $1}'`
	output="sig${output_label}_${iso_type}.root"
	if [ "$debug_mode" == "1" ] ; then 
	    echo "./create_skimmed_files_x $iso_type 0 'pythia_rs' $inputdir$file $outputdir$output"
	else
	    ./create_skimmed_files_x $iso_type 0 "pythia_rs" $inputdir$file $outputdir$output
	fi
    done
    echo "PYTHIA RS SKIMMED FILES HAVE BEEN CREATED"
fi

#---------------------------------------------------------------------------------------
##################   PYTHIA GAMMAGAMMA FILES - PYTHIA GAMMAGAMMA FILES
#---------------------------------------------------------------------------------------
if [ "$mc_gg" == "1" ] ; then 
    for file in `ls $inputdir* | grep input_mc12_gamgam | grep all | awk -F"${inputdir}" '{print $2}'`
    do
	output_label=`echo $file | awk -F"input" '{print $2}' | awk -F".txt" '{print $1}'`
	output="bkg${output_label}_${iso_type}.root"
	if [ "$debug_mode" == "1" ] ; then 
	    echo "./create_skimmed_files_x $iso_type 0 'pythia_gg' $inputdir$file $outputdir$output"
	else
	    ./create_skimmed_files_x $iso_type 0 "pythia_gg" $inputdir$file $outputdir$output
	fi
    echo "PYTHIA GAMMAGAMMA SKIMMED FILES HAVE BEEN CREATED"
    done
fi
#---------------------------------------------------------------------------------------
##################   PYTHIA GAMMAJET(+1FRAG) FILES - PYTHIA GAMMAJET(+1FRAG) FILES
#---------------------------------------------------------------------------------------
if [ "$mc_gj" == "1" ] ; then 
    #for file in `ls $inputdir* | grep input_mc12_gamjet | grep all | awk -F"${inputdir}" '{print $2}'`
    for file in `ls $inputdir* | grep input_mc12_gamjet | grep -v all | awk -F"${inputdir}" '{print $2}'`
    do
	output_label=`echo $file | awk -F"input" '{print $2}' | awk -F".txt" '{print $1}'`
	output="bkg${output_label}_${iso_type}.root"
	if [ "$debug_mode" == "1" ] ; then 
	    echo "./create_skimmed_files_x $iso_type 0 'pythia_gj' $inputdir$file $outputdir$output"
	else
	    #echo "./create_skimmed_files_x $iso_type 0 'pythia_gj' $inputdir$file $outputdir$output"
	    ./create_skimmed_files_x $iso_type 0 "pythia_gj" $inputdir$file $outputdir$output
	fi
    done
    echo "PYTHIA GAMMAJET SKIMMED FILES HAVE BEEN CREATED"
fi

