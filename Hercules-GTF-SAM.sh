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

#// Executes and Manages Hercules Distributed PySpark program GTF-SAM map step.
#// The GTF-MAP step  buckets ALL sequence alignment reads (from input
#// SAM file) into Exon (or CDS) feature regions (determined from input GTF
#// annotation file), this is upstream of the motif counting MapReduce steps
#// (MOTIF_POSITION MAP, MOTIF_POSITION_VECTOR_MAP, MOTIF_POSITION_REDUCE).

#//------------------------------------------------------------------------------
#// YARN "yarn-client" specific parameters
#// ** Crucial for distribution and optimisation **
_YARN_NUM_EXECUTORS_=5
_YARN_NUM_EXECUTOR_CORES_=20
#//------------------------------------------------------------------------------

_SEARCH_MOTIF_="ALL";
_JOB_PREFIX_="calibration2-"; # Trail with - for readability


#//------------------------------------------------------------------------------

_HDFS_SAM_READS_="hdfs://bigdata/user/jamie/calibration2.sam";


_HDFS_MAP_OUTPUT_="hdfs://bigdata/user/jamie/herc-txt-out-$_JOB_PREFIX_$_SEARCH_MOTIF_";
_HDFS_MAP_OUTPUT_CLEAN_="hdfs://bigdata/user/jamie/herc-txt-out-clean-$_JOB_PREFIX_$_SEARCH_MOTIF_.txt";
_HDFS_MAP_OUTPUT_MOTIF_="hdfs://bigdata/user/jamie/herc-txt-out-motif-$_JOB_PREFIX_$_SEARCH_MOTIF_";
_HDFS_MAP_OUTPUT_MOTIF2_="hdfs://bigdata/user/jamie/herc-txt-out-motif2-$_JOB_PREFIX_$_SEARCH_MOTIF_";
_HDFS_MAP_OUTPUT_MOTIF2_CLEAN_="hdfs://bigdata/user/jamie/herc-txt-out-motif2-clean-$_JOB_PREFIX_$_SEARCH_MOTIF_.txt";
_HDFS_REDUCE_OUTPUT_MOTIF_="hdfs://bigdata/user/jamie/herc-txt-out-motifcounts-$_JOB_PREFIX_$_SEARCH_MOTIF_";

_HDFS_MAP_OUTPUT_GC_="hdfs://bigdata/user/jamie/herc-txt-out-GC-$_JOB_PREFIX_$_SEARCH_MOTIF_";
_TEMP_DIR_="/home/local/jamie/Spark/_working/$_JOB_PREFIX_$_SEARCH_MOTIF_";

#//------------------------------------------------------------------------------


# Do GTF Stuff here (already done manually for now)
# GTF and SAM input files specified in Hercules script.

echo "WELCOME TO HERCULES"

#read -r -p "Run job? Previous HDFS data at $_HDFS_MAP_OUTPUT_ will be first removed...Are you sure? [y/N] " response
#if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
#then

	echo "initiating system bootstrap..."
	
	echo "removing existing data..."

	# Remove existing GTF_SAM_MAP data
	hadoop fs -rm -r $_HDFS_MAP_OUTPUT_
	
	# Ensure temp folder exists
	mkdir $_TEMP_DIR_

	# Execute the Map Job to assign key, value pairs to Aligned Reads (SAM) by using
	# GTF annotation to find the discrete regions of interest (Genes/Coding segments)
	echo "excuting the PySpark job [Operation: GTF_SAM_MAP] on the cluster..."	
	time spark-submit --master yarn-client --num-executors $_YARN_NUM_EXECUTORS_ --executor-cores $_YARN_NUM_EXECUTOR_CORES_ Hercules.py GTF_SAM_MAP $_SEARCH_MOTIF_ $_HDFS_SAM_READS_ $_HDFS_MAP_OUTPUT_

	rm $_TEMP_DIR_/intermediate.data
	
	echo "pulling off intermediate GTF_SAM map-step results from HDFS into $_TEMP_DIR_"
	time hadoop fs -getmerge $_HDFS_MAP_OUTPUT_ $_TEMP_DIR_/intermediate.data

	rm $_TEMP_DIR_/intermediate-clean.data
	# Remove any previous hadoop crc files which sometimes result in crc error on upload
    rm $_TEMP_DIR_/intermediate-clean.data.crc
	time ./Hercules-clean.sh $_TEMP_DIR_/intermediate.data $_TEMP_DIR_/intermediate-clean.data.tmp

	# Remove un-mapped (-1) reads for optimisation
	echo "Removing un-mapped (-1) reads for optimisation..."
	time ./Hercules-clean-GTF-SAM.sh intermediate-clean.data.tmp intermediate-clean.data
	rm intermediate-clean.data.tmp
	
	# Re-upload cleaned intermediate data to hdfs
	echo "Uploading cleaned intermediate data to hdfs: $_HDFS_MAP_OUTPUT_CLEAN_"
	time hadoop fs -put $_TEMP_DIR_/intermediate-clean.data $_HDFS_MAP_OUTPUT_CLEAN_

	# Remove current intermediate.data
	rm $_TEMP_DIR_/intermediate.data

	# Execute the GC-Map step to compute the median read GC content per exon as one lookup file
	# indexed by ExonID for the species
	echo "excuting the PySpark job [Operation: GC_MAP] on the cluster..."	
	time spark-submit --master yarn-client --num-executors $_YARN_NUM_EXECUTORS_ --executor-cores $_YARN_NUM_EXECUTOR_CORES_ Hercules.py GC_MAP $_SEARCH_MOTIF_ $_HDFS_MAP_OUTPUT_CLEAN_ $_HDFS_MAP_OUTPUT_GC_
	
	time hadoop fs -getmerge $_HDFS_MAP_OUTPUT_GC_ $_TEMP_DIR_/GC-content.data
	time ./Hercules-clean2.sh $_TEMP_DIR_/GC-content.data $_TEMP_DIR_/GC-content.csv

	# Remove temporary data
	rm $_TEMP_DIR_/GC-content.data

	
#else
#    echo "Nothing to do here. Bye";
#	exit;
#fi
