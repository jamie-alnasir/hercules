
HERCULES SYSTEM FOR ANALYSING K-MER MOTIFS IN TRANSCRIPTOMICS DATA

By Jamie Alnasir (PhD candidate)
Hugh Shanahan (PhD Supervisor)
Royal Holloway, University of London


PREREQUISITES

Apache Spark cluster configured to run PySpark jobs
GTF genome annotation file ** MUST be sorted on chromosome then feature start columns in ascending order
SAM aligned reads file

1) Transcriptomic feature i.e. CDS, EXON must be filtered from GTF so that the GTF contains
only those features to be used in the analysis. This can be achieved using "grep":

grep exon input.gtf >new.gtf

2) SAM reads file must not contain any commas ',' (sometimes found in the quality score).
Removal of this specific character can be achieved using the "tr" command:

tr -d ',' < input.sam >output.sam

3) ** GTF file must be sorted ** on chromosome and feature start columns in ascending order to enable
partitioning of the reads which uses a binary search optimisation. This can be achieved as follows:

cat input.gtf | sort -k1 | sort -k3n,3 $1 > output.gtf

CODE FILES AND THEIR PURPOSE

Phase I - for k-mer motif counting on Spark cluster:
Hercules.py				The PySpark code to partition reads by GTF feature and count motifs
Hercules-GTF-SAM.sh		The Bash script to run Hercules.py on Spark to partition reads
Hercules-MOTIF.sh		The Bash script to run Hercules.py on Spark count k-mer for a given Motif (parameter)
4mer.sh					The MAIN Bash script to run Hercules-MOTIF.sh once for each motif to count ALL k-mers

Phase II - post-cluster analysis of k-mer correlations:
Hercules-graphit.py				The plain Python code (requirements: Numpy, Scipy and Mattplotlib) to compute
                          k-mer correlations at various distances, by exon, and motif.

STEPS TO EXECUTE THE ANALYSIS

A. Set all parameters (mainly paths to input data such as _LFS_GTF_FILE_FILTERED_) at top of each code file.
B. Set number of Executors and Cores appropriate to cluster in the .sh scripts, these spawn Hercules PySpark jobs
on the Spark cluster


Running the job itself:

Phase I - On the Spark cluster:

1. Execute Hercules-GTF-SAM.sh to partition ALL the reads in input GTF-SAM file by features in GTF annotation file.
- This runs a single PySpark job to partition the reads.

2. Execute 4mer.sh to count ALL k-mers (example is AAAA through to GGGG, i.e. 4^4 = 256 permutations)
- This runs the Hercules-MOTIF.sh PySpark job 256 times on the cluster, once for each motif which is passed to
the job from 4mer.sh.

3. Final output will be on HDFS in the _HDFS_REDUCE_OUTPUT_MOTIF_ folder (specified in the Hercules-MOTIF.sh file),
but is copied automatically from HDFS to local fille system by the Getmerge command to the _TEMP_FOLDER_
(specified in the Hercules-MOTIF.sh file).

- The 256 k-mer result files are eached placed in a subfolder of the motif, i.e. AAAA and labeled i.e. herc-final-AAAA.csv
- A single csv file called GC-content.csv will be written to _TEMP_FOLDER_ which contains the mean GC content for
reads partitioned to each GTF feature.

Phase II -- Offline from cluster:

Run Hercules-graphit.py (this should work on data from the cluster stored on local file system in _TEMP_FOLDER_). It has
various parameters for different correlations calculations and motif plotting operations. See the code.
