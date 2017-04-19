#//==============================================================================
#//    _    _                     _           
#//   | |  | |                   | |          
#//   | |__| | ___ _ __ ___ _   _| | ___  ___ 
#//   |  __  |/ _ \ '__/ __| | | | |/ _ \/ __|
#//   | |  | |  __/ | | (__| |_| | |  __/\__ \
#//   |_|  |_|\___|_|  \___|\__,_|_|\___||___/
#//
#//==============================================================================
#// Hercules - PySpark framework for transcriptomic analysis
#// Centre for Systems and Synthetic Biology, Dept of Computer Science
#// Royal Holloway, University of London
#// By Jamie Al-Nasir, 07/2015 
#// Copyright (c) 2015 Jamie J. Alnasir, All Rights Reserved
#//==============================================================================
#// Version: Python (PySpark) edition
#//==============================================================================
#//
#// Change-log:
#//
#// 17/05/2016 - Added fork for MOTIF-MAP - read-start vs exon-start motif positions
#// 			 both modes produce the same format output, downstream procs. identical
#// 18/05/2016 - Added GC-MAP step to compute median GC content for reads in each exon


"""PDB-Spark.py"""
import sys;
import os;
import os.path;
import re;
from pyspark import SparkContext
from multiprocessing import Pool


#//------------------------------------------------------------------------------
# GTF Annotation file path (on local, shared filesystem)

# REGULAR JOB FILES

_LFS_GTF_FILE_FILTERED_  = "/home/local/jamie/hadoop-temp/flybase2006-exons.gtf";
#_LFS_GTF_FILE_FILTERED_  = "/home/local/jamie/hadoop-temp/mus-mus-wholebrain-r1.gtf";

# CALIBRATION JOB FILES
#_LFS_GTF_FILE_FILTERED_  = "/home/local/jamie/hadoop-temp/calibration.gtf";


# SAM reads file path
# NB: Changed to be an input parameter, no longer needs to be specified herein...
#_HDFS_SAM_READS_ = "hdfs://bigdata/user/jamie/Drosophila/sam/sam_reads.sam";
# NB: Changed to be an input parameter
#_SEARCH_MOTIF_
#//------------------------------------------------------------------------------

# Set to exclude preceeding nucleotides from analysis (Hansen's paper on Random Hexamer Primers)
_PRIMER_EXCLUDE_=True;
_READ_RELATIVE_PRIMER_EXCLUDE_=10;

# Nucleotide bases
bases = "ATGC";

# preload the GTF annotation
fr = open(_LFS_GTF_FILE_FILTERED_, 'r');	
gtf_lines = fr.readlines();
lstGTF_feat = [];
lstGTF_chr = [];

def getCol(lst, col):
	return [row[col] for row in lst];
	
def getChrList(chr):
	for chromes in lstGTF_chr:
		if (chromes[0] == chr):
			return lstGTF_feat[chromes[1]:chromes[2]];
			
def getMiddle(lst):
	if lst is None:
		return -1;
	return len(lst) // 2;
	

def readlinesFileSection(aFile, start, end):
# Read only lines of aFile specified by 0-based start and end (inclusive)
    lstR = [];
    fp = open(aFile);
    for i, line in enumerate(fp):
        if (i < start):
            continue;
        elif (i >= start) and (i <= end):
            lstR.append(line);
        elif (i > end):
            break;
    fp.close()
    return lstR;


def CigarToTupleList(aCigarStr):
# Parse a CIGAR string into a list of tuples [(length,operation),(length,operation),(length,operation)]
    return re.findall(r'(\d+)(\w)', aCigarStr);

def getFirstCigarAlignedMatchLen(aCigarStr):
# Return the Length of the first Alignment match 
# Kludge -- there can be multiple alignment matches separated by insertions, deletions or skips. We only use
# the first portion of the read that is an alignment match, regardless of whether it is a sequence match or not.
	lstCIG = CigarToTupleList(aCigarStr);
	if (lstCIG is None): return -1;
	for item in lstCIG:
		if (item[1] == 'M'):
			return long(item[0]);
	return -1;

def isCigarSingleAlignMatch(aCigarStr):
# Returns True if the CIGAR string is a single, contiguous aligned match
	lstCIG = CigarToTupleList(aCigarStr);
	if (lstCIG is None): return False;
	if (len(lstCIG) <> 1): return False;
	return (lstCIG[0][1] == 'M');

def _isPosWithinRange(feat_start, feat_end, pos):
	return (long(pos) >= long(feat_start)) and (long(pos) <= long(feat_end));

def _isWithinRange(feat_start, feat_end, ReadPos_start, ReadPos_end):
	# Old functionality
	# Return feature key if read start and end fall totally within the feature range
	#return long(ReadPos_start) >= long(feat_start) and long(ReadPos_end) <= long(feat_end);
	# New functionality (03/05/2016)
	# Return feature key if read falls within the feature range OR straddles feature range
	bLeftStraddle  = ( long(ReadPos_start) < long(feat_start) and _isPosWithinRange(feat_start, feat_end, ReadPos_end) );
	bRightStraddle = ( _isPosWithinRange(feat_start, feat_end, ReadPos_start) and long(ReadPos_end) > long(feat_end) );
	bStraddling    = (bLeftStraddle or bRightStraddle);
	#if bStraddling:
	#	print "straddling: ", ReadPos_start, ReadPos_end;
	bWithin 	   = long(ReadPos_start) >= long(feat_start) and long(ReadPos_end) <= long(feat_end);
	return bStraddling or bWithin;
	

def _gtf_lookup(lst, ReadPos_start, ReadPos_end):
# Recursive GTF lookup using binary search (sequential halving of the search list)
	m = getMiddle(lst);
	if (m == -1):
		return -1;
	if _isWithinRange(lst[m][1], lst[m][2], ReadPos_start, ReadPos_end):
		return lst[m][3]; 	
	if (m == 0):
		return -1;
	if (int(ReadPos_start) < lst[m][1]):
		return _gtf_lookup(lst[:m], ReadPos_start, ReadPos_end);
	if (int(ReadPos_start) > lst[m][1]):
		return _gtf_lookup(lst[m:], ReadPos_start, ReadPos_end);		
	

def gtf_lookup(chr, ReadPos_start, ReadPos_end):
# Function to look up the feature range within which a read lies using the input
# read and the GTF annotation using a binary search technique

	key = _gtf_lookup(getChrList(chr), ReadPos_start, ReadPos_end);
	if (key is None):
		return -1;
	return key;


def BaseCounts(aSeq):
# Count frequencies of bases in given sequence string
	lstBaseCounts = [];
	for base in bases:
		lstBaseCounts.append([base, len(filter(lambda x: (x == base), aSeq))]);
	return lstBaseCounts;


def GC_Content(aSeq):
# Return GC content of given sequence string
	lstB = BaseCounts(aSeq);
	dictB = dict(lstB);
	bases_all = int( dictB['A'] + dictB['T'] + dictB['G'] + dictB['C'] );
	bases_GC = int( dictB['G'] + dictB['C'] );
	gc = float(bases_GC) / float(bases_all) * 100;
	return gc;


def motifs(aSeqStr, aMotif):
# Return an array of motif string 0-based index positions within aSeqStr
	x = 0;
	lstResult = [];
	while (x <> -1):	
		if (x == 0):
			x = aSeqStr.find(aMotif, x);
			if (x <> -1):						
				lstResult.append(x);
				x = x + 1;							
				continue;
			else:
				return None;				
		else:	
			x = aSeqStr.find(aMotif, x + len(aMotif));						
		if (x <> -1):
			lstResult.append(x);			
	return lstResult;


	
def gtf_sam_mapper(sam_read_line):
	# GTF/SAM READ LOOKUP Map step
	# Sorts sequence alignment (sam read fragments) (v)
	# into GTF gene/coding regions (k)
	sam_data = sam_read_line.split('\t');
	chr  = sam_data[2];
	rpos = sam_data[3];
	rlen = getFirstCigarAlignedMatchLen(sam_data[5]);
	endpos = long(rpos) + long(rlen);
	#print 'looking up read @ pos ' + rpos + '-' + endpos;
	mapkey = gtf_lookup(chr, rpos, endpos);
#	if (mapkey is None):
#		_debug('None,' + sam_read_line);

	#print 'key-gen: ' + mapkey;
	yield (mapkey, sam_read_line);


def _motif_mapper(aReadLine, bReadRelative):
	# MOTIF Map step
	# Itemise Motif occurrence positions (offset from
	# start of read-range)
	# Set bReadRelative true for motif positions relative to read-start
	# otherwise positions are relative to exon-start
	lstMotifCorrect = None;
	print 'READLINE=' + aReadLine;
	key, value = aReadLine.strip('\'').split(',');
	key = key.strip('(');
	value = value.rstrip(')\n');
	if (long(key) <> -1):
		key_start = long(key[:10]);
		key_end = long(key[10:]);
	else:
		key_start = 0;
		key_end = 0;
	sam_data = value.split('\t');
	rpos = long(sam_data[3]);
	rlen = getFirstCigarAlignedMatchLen(sam_data[5]);
	endpos = rpos + rlen;		
	offset = rpos - key_start;
	keylen =  key_end - key_start;
	#print "SEQ=" + sam_data[9];
	lstMotif = motifs(sam_data[9], _SEARCH_MOTIF_);

	if (bReadRelative):
	# No offset correction as position is relative to read-start
		lstMotifCorrect = lstMotif;
	else:
		# Correct for read-offset so position is relative to exon-start
		if lstMotif is not None:					
			lstMotifOffset = map(lambda x: x+offset, lstMotif);
			lstMotifCorrect = filter(lambda x: (x > 0 and x <= keylen), lstMotifOffset);

	if lstMotifCorrect is None:
		lstMotifCorrect = [];

	#print "motifs=" + str(lstMotifCorrect);
	return key + ", " + str(lstMotifCorrect);


def motif_mapper_exon(aReadLine):	
# Perform exon-start relative motif map step
	return _motif_mapper(aReadLine, False);

def motif_mapper_read(aReadLine):	
# Perform read-start relative motif map step
	return _motif_mapper(aReadLine, True);

	
def gc_mapper(aReadLine):
	# GC_MAP Map step
	# Return the GC content of the given read
	key, value = aReadLine.strip('\'').split(',');
	key = key.strip('(');
	value = value.rstrip(')\n');
	sam_data = value.split('\t');
	gc = GC_Content(sam_data[9]);
	return (str(key), gc);

	
def vector_mapper(aVectorReadLine):
        # Motif Position map step
        # Input 1 value, return 1 or more value, must be called via flatMap
        #print "LINE: " + aVectorReadLine;
        key, value = aVectorReadLine.split(',', 1);
        key = key.strip();
        value = value.strip().lstrip('[').rstrip(']');
        lstVector = list(value.split(','));
        lstVector = [v.strip().replace('L','') for v in lstVector];
        lstVector = filter(lambda x: x.isdigit(), lstVector);
        for motif_pos in lstVector:
                #print key + "_" + motif_pos.strip() + ", 1";
                yield (key + "_" + motif_pos.strip().zfill(6).replace('L', ''), 1);
                #yield ("-1_0", 1);


def motif_reducer(v1, v2):
	return v1+v2;
	

print "WELCOME TO HERCULES!";	
	
# Spark Context
# argv[2] is the searcg motif, add to job context name
sc = SparkContext("yarn-client", sys.argv[2] + " Hercules");

# Load data for Binary Search method
cur_chr = "";
for x in range(0, len(gtf_lines)):
	chr, source, feature, feat_start, feat_end, score, strand, frame, attribute = gtf_lines[x].strip().split('\t');
	gtf_key = feat_start.zfill(10) + feat_end.zfill(10);
	lstGTF_feat.append([chr, int(feat_start), int(feat_end), gtf_key]);

	# Optimisation
	# Build list of chromosome ranges
	if (chr <> cur_chr):
			cur_chr = chr;
			lstGTF_chr.append([chr, x, -1]);
			if (x > 0):
				lstGTF_chr[len(lstGTF_chr) - 2][2] = x - 1;
lstGTF_chr[len(lstGTF_chr) - 1][2] = len(gtf_lines) - 1;

# Optimisation
# Pre-sort features within chromosomes
lstGTF_featTmp = [];
for chromes in lstGTF_chr:
	lstTmp = lstGTF_feat[chromes[1]:chromes[2]+1];
	lstGTF_featTmp.extend( sorted(lstTmp, key=lambda x: (x[0], x[1])) );
lstGTF_feat = [];
lstGTF_feat.extend(lstGTF_featTmp);
del lstGTF_featTmp;


def _safeCigarProc(read_str):
# Exception Wrapper
# Discard reads with non-predictable tab spacing
	try:
		return isCigarSingleAlignMatch(read_str.split('\t')[5]);
	except:
		return False;

def _safeFind(read_str, aMotif):
# wrapper
# Locates the motif position from within the SeqStr tab-delimited section of the read
        try:
                return read_str.split('\t')[9].find(aMotif);
        except:
                return -1;


if (len(sys.argv) > 4):
	_SEARCH_MOTIF_ = sys.argv[2]; # Pick up Search Motif

	if (sys.argv[1] == "GTF_SAM_MAP"):
		print "performing operation: GTF_SAM_MAP";
		_HDFS_SAM_READS_=sys.argv[3];
		_HDFS_OUT_=sys.argv[4];
	
		sam_reads = sc.textFile(_HDFS_SAM_READS_);
		
		# Filter for only single-continguous alignment matches
		sam_reads = sam_reads.filter(lambda read: _safeCigarProc(read) );
		
		if (_SEARCH_MOTIF_ == "ALL"):
			# No filtering of reads
			hercules_output = sam_reads.flatMap(gtf_sam_mapper);
		else:
			# Filter out only reads containing _SEARCH_MOTIF_
			# NB: Crude, whole-read filter, Should re-implement to only use the SEQ str, OK for now
			filtered_reads = sam_reads.filter(lambda read: _SEARCH_MOTIF_ in read);
			hercules_output = filtered_reads.flatMap(gtf_sam_mapper);

		hercules_output.saveAsTextFile(_HDFS_OUT_);
		
	elif (sys.argv[1] == "MOTIF_MAP"):		
		print "performing operation: MOTIF_MAP";
		_HDFS_MAPPED_READS_=sys.argv[3];
		_HDFS_OUT_=sys.argv[4];
	
		mapped_reads = sc.textFile(_HDFS_MAPPED_READS_);
		
		# filter motifs occuring within the primer region
	        if (_PRIMER_EXCLUDE_):
			filtered_reads = mapped_reads.filter(lambda read: _safeFind(read, _SEARCH_MOTIF_) >  _READ_RELATIVE_PRIMER_EXCLUDE_);
		else:
			# Use filtering for non-filtered GTF-MAPPED READS...
	                filtered_reads = mapped_reads.filter(lambda read: _SEARCH_MOTIF_ in read);

		hercules_output = filtered_reads.map(motif_mapper_exon);
		# else
		# hercules_output = mapped_reads.map(motif_mapper);

		hercules_output.saveAsTextFile(_HDFS_OUT_);

	elif (sys.argv[1] == "MOTIF_MAP_READ_RELATIVE"):		
		print "performing operation: MOTIF_MAP_READ_RELATIVE";
		_HDFS_MAPPED_READS_=sys.argv[3];
		_HDFS_OUT_=sys.argv[4];
	
		mapped_reads = sc.textFile(_HDFS_MAPPED_READS_);
		
		# Use filtering for non-filtered GTF-MAPPED READS...
		filtered_reads = mapped_reads.filter(lambda read: _SEARCH_MOTIF_ in read);		
		hercules_output = filtered_reads.map(motif_mapper_read);
		# else
		# hercules_output = mapped_reads.map(motif_mapper);

		hercules_output.saveAsTextFile(_HDFS_OUT_);	
		
	elif (sys.argv[1] == "VECTOR_MAP"):		
		print "performing operation: VECTOR_MAP";
		_HDFS_MOTIF_VECTORS_=sys.argv[3];
		_HDFS_OUT_=sys.argv[4];
	
		motif_vectors = sc.textFile(_HDFS_MOTIF_VECTORS_);
		hercules_output = motif_vectors.flatMap(vector_mapper);		

		hercules_output.saveAsTextFile(_HDFS_OUT_);

	elif (sys.argv[1] == "MOTIF_REDUCE"):
		print "performing operation: MOTIF_REDUCE";
		_HDFS_MOTIF_MAP2_=sys.argv[3];
		_HDFS_OUT_=sys.argv[4];
	
		motif_map2 = sc.textFile(_HDFS_MOTIF_MAP2_);
		
		# Clean data and package into tuples in the form of (k,1),(k,1)...for stock "add" reducer
		motif_map3 = motif_map2.map( lambda x: (x.split(',')[0], long(x.split(',')[1])) );
		
		# Perform ReduceByKey
		hercules_output = motif_map3.reduceByKey( lambda x,y: ( x+y )  );

		hercules_output.saveAsTextFile(_HDFS_OUT_);
		
	elif (sys.argv[1] == "GC_MAP"):		
		print "performing operation: GC_MAP";
		_HDFS_MAPPED_READS_=sys.argv[3];
		_HDFS_OUT_=sys.argv[4];
	
		mapped_reads = sc.textFile(_HDFS_MAPPED_READS_);
		
		if (_SEARCH_MOTIF_ == "ALL"):
			# Compute median read GC content per exon for ALL reads in the exon
			hercules_output = mapped_reads.map(gc_mapper);
		else:
			# Only compute median read GC content per exon for reads containing the motif
			filtered_reads = mapped_reads.filter(lambda read: _SEARCH_MOTIF_ in read);
			hercules_output = filtered_reads.map(motif_mapper_exon);

		# Compute median by combineByKey method
		sumCount = hercules_output.combineByKey(lambda value: (value, 1),
						 lambda x, value: (x[0] + value, x[1] + 1),
						 lambda x, y: (x[0] + y[0], x[1] + y[1]))

		gc_averageByKey = sumCount.map(lambda (label, (value_sum, count)): (label, value_sum / count));
		gc_averageByKey.saveAsTextFile(_HDFS_OUT_);

	else:
		print "Invalid mode specified";
		print "Usage sparksubmit ... Hercules.py <GTF_SAM_MAP|MOTIF_MAP> SEQMOTIF> <hdfs_in_dir> <hdfs_out_dir>";
		
else:
	print "Usage sparksubmit ... Hercules.py <GTF_SAM_MAP|MOTIF_MAP> SEQMOTIF> <hdfs_in_dir> <hdfs_out_dir>"
	
