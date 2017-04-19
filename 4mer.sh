#!/bin/bash

# Format $_SUB_DIR_-$_JOB_PREFIX_ (i.e. 4mer-mut-r2-)
#_LOG_DIR_="/home/local/jamie/Spark/_working/4mer-mus-mus-wholebrain-repl1-hex-/logs/";
_LOG_DIR_="/home/local/jamie/Spark/_working/4mer-calibration2-/logs/";


mkdir -p $_LOG_DIR_

while read fourmer; do
  echo "Executing Spark Hercules job for motif 4-mer $fourmer";
  _LOG_FILE_=$_LOG_DIR_$fourmer.log;

  echo "logging output to $_LOG_FILE_";

  ./Hercules-MOTIF.sh $fourmer >$_LOG_FILE_;

  echo "job ($fourmer) completed or terminated."

done <4mer.txt
