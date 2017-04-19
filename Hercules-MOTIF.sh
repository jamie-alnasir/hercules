#!/bin/bash
#//==============================================================================
#// Hercules - PySpark framework for transcriptomic analysis
#// Centre for Systems and Synthetic Biology, Dept of Computer Science
#// Royal Holloway, University of London
#// By Jamie Alnasir, 07/2015 
#// Copyright (c) 2015 Jamie J. Alnasir, All Rights Reserved
#//==============================================================================
#// Version: Bash boot-strap
#//==============================================================================

#// Executes and Manages Hercules Distributed PySpark program motif counting
#// MapReduce steps (MOTIF_POSITION MAP, MOTIF_POSITION_VECTOR_MAP,
#// MOTIF_POSITION_REDUCE).


#//------------------------------------------------------------------------------
#// YARN "yarn-client" specific parameters
#// ** Crucial for distribution and optimisation **
_YARN_NUM_EXECUTORS_=5
_YARN_NUM_EXECUTOR_CORES_=20
#//------------------------------------------------------------------------------

_SEARCH_MOTIF_=$1;
#_JOB_PREFIX_="mus-mus-wholebrain-repl1-hex-"; # Trail with - for readability
_JOB_PREFIX_="calibration2-"; 
_SUB_DIR_="4mer";

#//------------------------------------------------------------------------------

# Input, cleaned data, single text file from GTF_SAM_MAP
_HDFS_MAP_OUTPUT_CLEAN_="hdfs://bigdata/user/jamie/herc-txt-out-calibration2-ALL.txt";

_HDFS_MAP_OUTPUT_MOTIF_="hdfs://bigdata/user/jamie/$_SUB_DIR_/herc-txt-out-motif-$_JOB_PREFIX_$_SEARCH_MOTIF_";
_HDFS_MAP_OUTPUT_MOTIF2_="hdfs://bigdata/user/jamie/$_SUB_DIR_/herc-txt-out-motif2-$_JOB_PREFIX_$_SEARCH_MOTIF_";
_HDFS_MAP_OUTPUT_MOTIF2_CLEAN_="hdfs://bigdata/user/jamie/$_SUB_DIR_/herc-txt-out-motif2-clean-$_JOB_PREFIX_$_SEARCH_MOTIF_.txt";
_HDFS_REDUCE_OUTPUT_MOTIF_="hdfs://bigdata/user/jamie/$_SUB_DIR_/herc-txt-out-motifcounts-$_JOB_PREFIX_$_SEARCH_MOTIF_";
_TEMP_DIR_="/home/local/jamie/Spark/_working/$_SUB_DIR_-$_JOB_PREFIX_/$_SEARCH_MOTIF_";

#//------------------------------------------------------------------------------


# Do GTF Stuff here (already done manually for now)
# GTF and SAM input files specified in Hercules script.

echo "WELCOME TO HERCULES"

#read -r -p "Run job? Previous HDFS data at $_HDFS_MAP_OUTPUT_ will be first removed...Are you sure? [y/N] " response
#if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
#then

	echo "initiating system bootstrap..."
	
	echo "removing existing data..."

	# Remove current intermediate.data
	# rm $_TEMP_DIR_/intermediate.data

	# Remove existing GTF_SAM_MAP data
	hadoop fs -rm -r $_HDFS_MAP_OUTPUT_MOTIF_

	# Motif Position Map
	echo "excuting the PySpark job [Operation: MOTIF_MAP] on the cluster..."        
	time spark-submit --master yarn-client --num-executors $_YARN_NUM_EXECUTORS_ --executor-cores $_YARN_NUM_EXECUTOR_CORES_ Hercules.py MOTIF_MAP $_SEARCH_MOTIF_ $_HDFS_MAP_OUTPUT_CLEAN_ $_HDFS_MAP_OUTPUT_MOTIF_

	# Motif Position Vector Map
	echo "excuting the PySpark job [Operation: VECTOR_MAP] on the cluster..."       
	time spark-submit --master yarn-client --num-executors $_YARN_NUM_EXECUTORS_ --executor-cores $_YARN_NUM_EXECUTOR_CORES_ Hercules.py VECTOR_MAP $_SEARCH_MOTIF_ $_HDFS_MAP_OUTPUT_MOTIF_  $_HDFS_MAP_OUTPUT_MOTIF2_
	

	echo "pulling off intermediate MOTIF_MAP vector map-step results from HDFS into $_TEMP_DIR_"
	time hadoop fs -getmerge $_HDFS_MAP_OUTPUT_MOTIF2_ $_TEMP_DIR_/intermediate-vectormap.data

	rm $_TEMP_DIR_/intermediate-clean.data
	# Remove any previous hadoop crc files which sometimes result in crc error on upload
    rm $_TEMP_DIR_/intermediate-vectormap-clean.data.crc
	time ./Hercules-clean2.sh $_TEMP_DIR_/intermediate-vectormap.data $_TEMP_DIR_/intermediate-vectormap-clean.data

	# Re-upload cleaned intermediate data to hdfs
	echo "Uploading cleaned intermediate data to hdfs: $_HDFS_MAP_OUTPUT_MOTIF2_CLEAN_"
	time hadoop fs -put $_TEMP_DIR_/intermediate-vectormap-clean.data $_HDFS_MAP_OUTPUT_MOTIF2_CLEAN_

	# Remove current intermediate-vectormarp-clean.data
        rm $_TEMP_DIR_/intermediate-vectormap.data
	
	# Remove existing MOTIF_REDUCE data
	hadoop fs -rm -r $_HDFS_REDUCE_OUTPUT_MOTIF_

	# Motif Map count ReduceByKey step
	echo "excuting the PySpark job [Operation: MOTIF_REDUCE] on the cluster..."	
	time spark-submit --master yarn-client --num-executors $_YARN_NUM_EXECUTORS_ --executor-cores $_YARN_NUM_EXECUTOR_CORES_ Hercules.py MOTIF_REDUCE $_SEARCH_MOTIF_ $_HDFS_MAP_OUTPUT_MOTIF2_CLEAN_ $_HDFS_REDUCE_OUTPUT_MOTIF_
	
	echo "pulling off final MOTIF_REDUCE results from HDFS into $_TEMP_DIR_"
	time hadoop fs -getmerge $_HDFS_REDUCE_OUTPUT_MOTIF_ $_TEMP_DIR_/final-tmp.data
	
	
	# Clean-up final data
	time ./Hercules-clean2.sh $_TEMP_DIR_/final-tmp.data $_TEMP_DIR_/final-tmp2.data
	tr -d "L" < $_TEMP_DIR_/final-tmp2.data >$_TEMP_DIR_/final-tmp3.data
	time sort $_TEMP_DIR_/final-tmp3.data >$_TEMP_DIR_/final-tmp4.data
	time ./Hercules-clean3.sh $_TEMP_DIR_/final-tmp4.data $_TEMP_DIR_/herc-final-$_SEARCH_MOTIF_.csv

	# Remove intermediate data because we maybe running many jobs...
	rm $_TEMP_DIR_/intermediate-clean.data
	rm $_TEMP_DIR_/intermediate-vectormap-clean.data

#else
#    echo "Nothing to do here. Bye";
#	exit;
#fi
