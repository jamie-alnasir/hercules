#!/usr/bin/python
#//==============================================================================
#// Hercules - PySpark framework for transcriptomic analysis
#// Centre for Systems and Synthetic Biology, Dept of Computer Science
#// Royal Holloway, University of London
#// By Jamie Alnasir, 07/2015 
#// Copyright (c) 2015 Jamie J. Alnasir, All Rights Reserved
#//==============================================================================
#// Version: Python edition
#//==============================================================================
#// Change-log:
#// 22/03/2017 - Fixed bug in the way in which SCATTERMATRIX plot GC concentrations
#// 			 are computed. The spacing correlations are computed for each exon
#// 			 group rather than by motif group. This fixed the problem with
#// 			 spread of Motif GC as a function of mean exon GC content of the
#// 			 reads (proxy for mean read GC content). ** This change was wrong!
#//
#// 24/03/2017 - Changed the above fix which was incorrect. We can't group by exon
#// 			 as there would be insufficient data, causing poor correlation.
#// 			 We therefore compute for Motif/spacing and then bin the reads
#// 			 by filtering for exon-GC.
#// NB:			 If _GC_FILTERING_ = True, *MUST* define following variables as
#//				 writable global variables: _GC_FILTER_MIN_, _GC_FILTER_MAX_
#//				 in the method that sets them for filtering, otherwise filter will
#//				 not be applied. -- They are set as globally writable in main method.
#//
#// 27/03/2017 - Removed max cap for COUNTS > 1000.

 
#// Module to process Hercules Spark MapReduce output data (exon motif positions and
#// counts) and filter/format for plotting
import sys;
import os.path;
import copy;
import numpy;
import pprint;
import simplejson;
import time;
import random;
from itertools import groupby;
 
# Stats libraries
import numpy as np
from scipy import stats;
from scipy.special import stdtr;
 
# Plotting libraries
import matplotlib;
matplotlib.use('Agg');
import matplotlib.pyplot as plt;
 
bases = "ATGC";
#_CAP_COUNT_=1000; # REMOVED THIS FEATURE
_NAN_VALUE_=float('nan');#8;
 
#_JOB_FRIENDLY_TITLE_ = "[Wild-type D. melanogaster]";
#_JOB_FRIENDLY_TITLE_ = "[Mutant-r2 type D. melanogaster]";
#_JOB_FRIENDLY_TITLE_ = "[Mus-musculus Wholebrain-r1]";
#_JOB_FRIENDLY_TITLE_ = "[Mus-musculus Adrenal-r1]";
#_JOB_FRIENDLY_TITLE_ = "[Synthetic dataset]";
 
#_SPECIES_LABEL_ = "Mutant-r2-type"; # For latex labels, capitalize first letter
#_SPECIES_LABEL_ = "Wild-type"; # For latex labels, capitalize first letter
 _SPECIES_LABEL_ = "Mus-musculusAdrenalR1"; # For latex labels, capitalize first letter ** ONE WORD ONLY OR CamelCase
#_SPECIES_LABEL_ = "Mutant-r2-type"; # For latex labels, capitalize first letter
#_SPECIES_LABEL_ = "Mus-musculus Wholebrain-r1"; # For latex labels, capitalize first letter


# N-mer mode, use 4 for Fourmer, 5 for Fivemer etc...
_NMER_MODE_=4;
_NMER_JOB_PATH_ = "/home/local/mxba001/PhD/_working/4mer-wild-/";
#_NMER_JOB_PATH_ = "/home/local/mxba001/PhD/_working/4mer-mutant-r2-/";
#_NMER_JOB_PATH_ = "/home/local/mxba001/PhD/_working/4mer-mus-mus-wholebrain-r1-/";
#_NMER_JOB_PATH_ = "/home/local/mxba001/PhD/_working/4mer-mus-mus-r1-/";
#_NMER_JOB_PATH_ = "/home/local/mxba001/PhD/_working/4mer-mus-mus-r2-/";
#_NMER_JOB_PATH_ = "/home/local/mxba001/PhD/_working/4mer-calibration2-/";
 
# Path to csv file of GC content for all exons of the organism
#_GC_FILE_PATH_ = "/home/local/mxba001/PhD/_working/4mer-wild-/GC-content.csv";
#_GC_FILE_PATH_ = "/home/local/mxba001/PhD/_working/4mer-mutant-r2-/GC-content.csv";
#_GC_FILE_PATH_ = "/home/local/mxba001/PhD/_working/4mer-mus-mus-wholebrain-r1-/GC-content.csv";
_GC_FILE_PATH_ = "/home/local/mxba001/PhD/_working/4mer-mus-mus-r1-/GC-content.csv";
#_GC_FILE_PATH_ = "/home/local/mxba001/PhD/_working/4mer-calibration-/GC-content.csv";
 
# Set to true to read in Read-Relative motif position file [DEPRECATED]
_READ_RELATIVE_=False;
  
# For exon reads median GC content lookup
dictGC = {};
 
#lstExonGC_0_25 = [];
#lstExonGC_25_50 = [];
#lstExonGC_50_75 = [];
#lstExonGC_75_100 = [];
 
# Enable for Filtering reads/counts based on exon % GC content
_GC_FILTERING_ = False;
_GC_FILTER_MIN_ = 60;
_GC_FILTER_MAX_ = 70;
 
# Set minimum number of reads (pairs of counts) to compute the correlation, typically about 10.
_MIN_CORRELATION_COUNTS_ = 10; # NORMALLY 10
 
# Plotting constants
_PLOT_SCATTER_PT_SIZE_ = 4;
_PLOT_OUT_PATH_ = "/home/local/mxba001/Dropbox/PySpark/Results/Drosophila/wild-type/plots/"; # add trailing /
#_PLOT_OUT_PATH_ = "/home/local/mxba001/Dropbox/PySpark/Results/Drosophila/mutant-r2-type/plots/"; # add trailing /
#_PLOT_OUT_PATH_ = "/home/local/mxba001/Dropbox/PySpark/Results/Mus-musculus/wholebrain-repl1/plots/"; # add trailing /
#_PLOT_OUT_PATH_ = "/home/local/mxba001/Dropbox/PySpark/Results/Mus-musculus/adrenal-repl1/plots/"; # add trailing /
#_PLOT_OUT_PATH_ = "/home/local/mxba001/Dropbox/PySpark/Results/Mus-musculus/adrenal-repl2/plots/"; # add trailing /
#_PLOT_OUT_PATH_ = "/home/local/mxba001/Dropbox/PySpark/Results/Drosophila/Calibration/plots/"; # add trailing /


# Correlations data
_CORREL_FILE_ = "./Correls-wild.csv";
#_CORREL_FILE_ = "./Correls-mutant-r2.csv";
#_CORREL_FILE_ = "./Correls-calibration.csv";
#_CORREL_FILE_ = "./Correls-mus-mus-wholebrain-r1.csv";#
#_CORREL_FILE_ = "./Correls-mus-mus-wholebrain-RenLab-r1.csv";
#_CORREL_FILE_ = "./Correls-mus-mus-r1.csv";
#_CORREL_FILE_ = "./Correls-mus-mus-r2.csv";

# Caching of PEARSON-TABLE1 results for certain jobs
_PEARSON_TABLE_CACHE_ = False;
#_PEARSON_TABLE_CACHE_FILE_ = "/home/local/mxba001/Dropbox/PySpark/wild.data";
#_PEARSON_TABLE_CACHE_FILE_ = "/home/local/mxba001/Dropbox/PySpark/mutant-r2.data";
#_PEARSON_TABLE_CACHE_FILE_ = "/home/local/mxba001/Dropbox/PySpark/Motif-calibration.data";
#_PEARSON_TABLE_CACHE_FILE_ = "/home/local/mxba001/Dropbox/PySpark/wholebrain-r1-.data";
lstPearsonTblCache = [];

# Latex Generation
_LATEX_MODE_ = True;
lstLatex = [];



def SaveData(aOutFile, data):
# Save as SimpleJSON dump
	with open(aOutFile, 'w') as f:
		simplejson.dump(data, f);
 
def LoadData(aInFile):
# Load from SimpleJSON dump
	with open(aInFile) as f:
		data = simplejson.load(f);
		return data;
		
		
def svLinesToArray(data, svChar):
# Convert array of separated value string lines into list of lists
	lstR = [];
	for x in data:
		lstR.append(x.strip().split(svChar));
		# Individually trim the separated values
		for y in range(0, len(lstR[-1])):
			lstR[-1][y] = lstR[-1][y].strip();
	return lstR;

def loadText(aInFile):
# Load array of strings from text file
	with open(aInFile) as f:
		data = f.readlines();
		return data;
		
def saveText(aOutFile, data):
# Load array of strings from text file
	with open(aOutFile, "w") as f:
		for i in data:
			f.write(i + "\n");	

def saveNumArray(aOutFile, data):
# Load array of strings from text file
	with open(aOutFile, "w") as f:
		for i in data:
			f.write(str(i) + "\n");	

			
def LoadGC(gc_file):
	global dictGC;
	fr = open(gc_file, 'r');
	lstGC = fr.read().splitlines();
	for i in lstGC:
		exon, gc_str = i.split(',');
		dictGC[exon] = float(gc_str);
		 
def GC_Lookup(exonIDStr):
	return dictGC[exonIDStr];
 
 
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
 
def getCol(lst, col):
	return [row[col] for row in lst];
	
def sliceLargerBySmaller(lstA, lstB):
# return lstA, lstB sliced to match the length of the smaller list
	if len(lstA) > len(lstB):
		lstA = lstA[0: len(lstB)];
	if len(lstB) > len(lstA):
		lstB = lstB[0: len(lstA)];
	return lstA, lstB;
	 
def prettyFloat(afloat):
	return "{0:.4f}".format(afloat);
 
def chunks(l, n):
	n = max(1, n);
	return [l[i:i + n] for i in range(0, len(l), n)];
	
def removeNaNs(lstCorrels):
	return filter(lambda x: x == x, lstCorrels); # NB non-nans are equal to themselves
	
def cleanR(r):
	if r != r:
		return _NAN_VALUE_;
	else:
		return r;
	
def median(lstV):
# Python 2.7 has no median function, use numpy's
	v = np.array(lstV);
	return np.median(v);
	
def Q1(lstV):
# Use numpy to calculate Q1 (q25)	
	q75, q25 = np.percentile(lstV, [75 ,25], interpolation = 'higher');
	#return 6; # INVESTIGATION 21/03/2017 -- TEST;
	return q25;
	
 
def calcFirstQuartile(sorted_motifs_file):
# Calculate the First Quartile (value between lowest value and median) for ALL of the counts EXCEPT
# those which were not associated with an exon (have Exon key of -1 from the Hercules Map-Reduce Job).
	global dictMotifs;
	global FIRST_QUARTILE;
	f = open(sorted_motifs_file, 'r');
	motifs_raw = f.read().splitlines();
 
	#print sorted_motifs_file;
	 
	lstCounts = [];
 
	for motif_line in motifs_raw:
		exon, pos, count = motif_line.split(',');
		try:
			# Only interested in
			if (exon <> "-1"):				
				lstCounts.append(int(count));
		except:
			(1 == 1);
			#print "error: " + motif_line;

 
	if (len(lstCounts) == 1):
		print sorted_motifs_file;
	 	 
	#c_avg = float( sum(lstCounts) ) / float( len(lstCounts) );
	# Calc 1st quartile
	c_fq = Q1(lstCounts);
	return c_fq;
 
def removeNoise(exon_motif_reads):
	if (exon_motif_reads == []):
		return [];
 
	global _CAP_COUNT_;
	global FIRST_QUARTILE;
	lstWorking = [];
	counts = [];
	#print "FIRST QUARTILE = ", FIRST_QUARTILE;
	# itertools group iterator, pull off array
	for item in exon_motif_reads:
		counts.append(item[2]);
		lstWorking.append([item[0], item[1], item[2]]);
 
	# Calc 1st quartile (OLD CALCULATION, NOT CORRECT ACCORDING TO Q1 DEFINITION!)
	#c_avg = sum(counts)/len(counts);
	#c_min = min(counts);
	#c_fq  = float(c_avg - c_min)/2;
	
	# USE FIRST_QUARTILE GLOBAL CALCULATED WITH NUMPY
	#return filter(lambda x: (x[2] >= FIRST_QUARTILE and x[2] <= _CAP_COUNT_), lstWorking); # REMOVED CAP
	return filter(lambda x: (x[2] >= FIRST_QUARTILE), lstWorking);

	
def meanGC(lstExons):
# Return the mean GC content of all the exons supplied

	lstR = [];
	for e in lstExons:
		lstR.append(GC_Lookup(e));

	return np.mean(lstR);

# Code from StatsModels library -- source code to avoid dependencies
def _ecdf(x):
    '''no frills empirical cdf used in fdrcorrection;
    ''';
    nobs = len(x);
    return np.arange(1,nobs+1)/float(nobs);

def fdrcorrection(pvals, alpha=0.05, method='indep'):
    pvals = np.asarray(pvals);
    pvals_sortind = np.argsort(pvals);
    pvals_sorted = pvals[pvals_sortind];
    sortrevind = pvals_sortind.argsort();
    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = _ecdf(pvals_sorted);
    elif method in ['n', 'negcorr']:
        cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))  # corrected this
        ecdffactor = _ecdf(pvals_sorted) / cm;
    else:
        raise ValueError('only indep and necorr implemented');
    reject = pvals_sorted <= ecdffactor*alpha;
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0]);
        reject[:rejectmax] = True;
    pvals_corrected_raw = pvals_sorted / ecdffactor;
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1];
    pvals_corrected[pvals_corrected > 1] = 1;
    return reject[sortrevind], pvals_corrected[sortrevind]; 

	
def _safefdrcorrection(pvals):
# The purpose of this function is to remove Nans from the source list for the calculations
# perform the FDR correction histogram calc then re-add the nans at the same index so as
# not to affect the FDR calculation.

	filtered = [];
	lstNan = [];      # store position of NaNs
	for i in range(0, len(pvals)):
		if (pvals[i] == pvals[i]): # test for non NaN
			filtered.append(pvals[i]);
		else:
			lstNan.append(i); # position of NaN
	result = fdrcorrection(filtered);

	if len(lstNan) > 0:
		#print pvals;
		#print "NAN=", lstNan;
		for i in lstNan:
			# make new tuple because they're immutable
			r1 = np.insert(result[1], i, float('NaN')); # re-add the NaN
			r2 = np.insert(result[0], i, False); # re-add the Corresponding Bool
			result = [r2,r1];

	return result;
	

def filterByGC(exon_motif_reads, minGC, maxGC):
# Filter exon counts by their GC content
	lstWorking = [];
	for item in exon_motif_reads:
		lstWorking.append([item[0], item[1], item[2]]);
 
	#print len(lstWorking), minGC, maxGC;
		 
	# return only exon counts for exons with median GC content within minGC and maxGC
	lstR =  filter(lambda x: (GC_Lookup(x[0]) >= minGC and GC_Lookup(x[0]) < maxGC), lstWorking);
	#print len(lstR);
	return lstR;

	 
def filterPosition(exon_motif_reads, position_bp, read_tolerance_bp):
	# Filter exons with motifs occurring at given position +/- tolerance
	lstWorking = exon_motif_reads;
	return filter(lambda x: (x[1] >= (position_bp - read_tolerance_bp) and x[1] <= (position_bp + read_tolerance_bp)), lstWorking);
 
 
def filter_and_combinePositions(exon_motif_reads, position1_bp, position2_bp, read_tolerance_bp):
	# Filter exons with motifs occurring at given position1 AND position2 +/- tolerance
	# Value copy exon lists for two separate filtering operations
	lstExons1 = list(exon_motif_reads);
	lstExons2 = list(lstExons1);
	 
	lstPos1 = filterPosition(lstExons1, position1_bp, read_tolerance_bp);
	lstPos2 = filterPosition(lstExons2, position2_bp, read_tolerance_bp);
	 
	if len(lstPos1)>0 and len(lstPos2)>0:
		lstResult = [];
		lstResult.append(lstPos1[0] + lstPos2[0]);
		return lstResult;
	else:
		return [];
	 
	 
def pairsWithinSpace(exon_motif_reads, space_bp, tolerance_bp):
	# Filter exons with motifs occurring within a given spacing (+/- tolerance)
	# Value copy exon lists for two separate filtering operations
	lstExons1 = list(exon_motif_reads);
	lstExons2 = list(lstExons1);
 
	lstResult = [];
	for i in range (0, len(lstExons1)):
		for j in range (i+1, len(lstExons2)):
			#print lstExons1[i], lstExons2[j];
			if (abs(lstExons2[j][1] - lstExons1[i][1]) >= space_bp-tolerance_bp) and (abs(lstExons2[j][1] - lstExons1[i][1]) <= space_bp+tolerance_bp):
					lstResult.append(lstExons1[i] + lstExons2[j]);
	 
	return lstResult;	   
 
 
def selectMotifSignals(sorted_motifs_file):
	global dictMotifs;  
	f = open(sorted_motifs_file, 'r');
	motifs_raw = f.read().splitlines();
 
	lstMotifs = [];
 
	for motif_line in motifs_raw:
		exon, pos, count = motif_line.split(',');
		try:
			lstMotifs.append([exon, int(pos), int(count)]);
		except:
			(1 == 1);
			#print "error: " + motif_line;

	for key, group in groupby(lstMotifs, lambda x: x[0]):
		if (key == "-1"):
			continue;
 
		#print key;
 
		#filtered_group = filterByGC(filtered_group, 100, 150);
		filtered_group = removeNoise(group);
		 
		#for count in filtered_group:
		#   print count;
		counts = getCol(filtered_group,2);
		if len(counts) > 1:
			chunk = chunks(counts, 2)[0];
			print chunk;
 
 
def selectMotifSignalsPosition(sorted_motifs_file, position_bp, tolerance_bp):
# Return motif occurrences on the same exon at a given position (+/- tolerance)
	global dictMotifs;  
	f = open(sorted_motifs_file, 'r');
	motifs_raw = f.read().splitlines();
 
	lstMotifs = [];
 
	for motif_line in motifs_raw:
		exon, pos, count = motif_line.split(',');
		try:
			lstMotifs.append([exon, int(pos), int(count)]);
		except:
			(1 == 1);
			#print "error: " + motif_line;
 
	for key, group in groupby(lstMotifs, lambda x: x[0]):
		if (key == "-1"):
			continue;
 
		#filtered_group = filterByGC(group, 100, 150);
		filtered_group = removeNoise(group);
		filtered_group_pos = filterPosition(filtered_group, position_bp, tolerance_bp);
 
		#print len(filtered_group);
		#print len(filtered_group_pos);
	 
		for item in filtered_group_pos:
			print item[0] + ", " + str(item[1]) + ", " + str(item[2]);		  
 
 
def selectMotifSignalsPositions(sorted_motifs_file, position1_bp, position2_bp, tolerance_bp):
# Return motif occurrences on the same exon at a given two positions (+/- tolerance)
	global dictMotifs;  
	f = open(sorted_motifs_file, 'r');
	motifs_raw = f.read().splitlines();
 
	lstMotifs = [];
 
	for motif_line in motifs_raw:
		exon, pos, count = motif_line.split(',');
		try:
			lstMotifs.append([exon, int(pos), int(count)]);
		except:
			(1 == 1);
			#print "error: " + motif_line;
			 
	for key, group in groupby(lstMotifs, lambda x: x[0]):
		if (key == "-1"):
			continue;
 
		#filtered_group = filterByGC(group, 100, 150);
		filtered_group = removeNoise(group);
		filtered_group_pos = filter_and_combinePositions(filtered_group, position1_bp, position2_bp, tolerance_bp);
 
		#print len(filtered_group);
		#print len(filtered_group_pos);
	 
		for item in filtered_group_pos:
			print item[0] + ", " + str(item[1]) + ", " + str(item[2]) + ", " + item[3] + ", " + str(item[4]) + ", " + str(item[5]);
 

def _motifsBySpace(sorted_motifs_file, space_bp, tolerance_bp):
# Return motif occurrences on the same exon within a given spacing
	global dictMotifs;  
	f = open(sorted_motifs_file, 'r');
	motifs_raw = f.read().splitlines();
 
	lstMotifs = [];
	lstResults = [];
 
	for motif_line in motifs_raw:
		exon, pos, count = motif_line.split(',');
		try:
			lstMotifs.append([exon, int(pos), int(count)]);
		except:
			(1 == 1);
			#print "error: " + motif_line;
	 
	# Filter reads/counts based on exon GC content
	if _GC_FILTERING_:
		#print "GC-FILTERING OF counts for the motif by EXON GC:", _GC_FILTER_MIN_, _GC_FILTER_MAX_;
		lstMotifs = filterByGC(lstMotifs, _GC_FILTER_MIN_, _GC_FILTER_MAX_);
 

	for key, group in groupby(lstMotifs, lambda x: x[0]):
		if (key == "-1"):
			continue;
 
		filtered_group = removeNoise(group);
		filtered_group2 = pairsWithinSpace(filtered_group, space_bp, tolerance_bp);
 
	 
		for item in filtered_group2:
			#print item[0] + ", " + str(item[1]) + ", " + str(item[2]) + ", " + item[3] + ", " + str(item[4]) + ", " + str(item[5]);
			lstResults.append(item);
			 
	return lstResults;
			 
			 
def motifsBySpace(sorted_motifs_file, space_bp, tolerance_bp):
# Return motif occurrences on the same exon within a given spacing
	 
	lstResults = _motifsBySpace(sorted_motifs_file, space_bp, tolerance_bp);
	 
	for item in lstResults:
		print item[0] + ", " + str(item[1]) + ", " + str(item[2]) + ", " + item[3] + ", " + str(item[4]) + ", " + str(item[5]);


def _safeCalcPearson(lstExonCounts, lstCounts1, lstCounts2):
# Calculate Pearson Co-efficient

	if (len(lstExonCounts) >= _MIN_CORRELATION_COUNTS_):
		try:
			Rvalue = numpy.corrcoef(lstCounts1, lstCounts2)[0, 1];		
		except:
			Rvalue = 0;
			print "NUMPY ERR:", len(lstCounts1);
	else:
		Rvalue = 0;
		
	return Rvalue;
		


		
def rawCorrelationsData(motif, motifFile, bPrint):
# Return array of simple tuples for motif for Mean ExonGC, MotifGC, Spacing, R-value
# There will be 4 tuples in the result, one for each spacing.

	global _GC_FILTERING_;
	global _GC_FILTER_MIN_, _GC_FILTER_MAX_; # to be writable to set filter range

	# We'll be cycling the filter range for the computations
	_GC_FILTERING_=True;

	# GC content bins
	lstREAD_GC_bins  = [[30,40], [40,50], [50,60], [60,70]];

	lstRawCorrelData = [];


	dictRval = {10:0, 50:0, 100:0, 200:0}
	for spacing in [10,50,100,200]:
	
		if spacing in [10,50]:
			tol = 2;
		else:
			tol = 4;

		# filter by mean GC content of the exon reads
		for gc_range in lstREAD_GC_bins:

			_GC_FILTER_MIN_ = gc_range[0];
			_GC_FILTER_MAX_ = gc_range[1];

			bin_index_exon_gc = gc_range[0];
			#print "Filtering for Median Exon GC content range: ", gc_range[0], '-', gc_range[1];
			MeanExonGC = (float(gc_range[1] - gc_range[0]) / 2) + gc_range[0];

			# n-spaced n-mers by exon
			#print "Mean GC=", MeanExonGC, " Motif GC=", GC_Content(motif), " spacing=", spacing;
			lstMotifs = _motifsBySpace(motifFile, spacing, tol); # NB filtered by Read/Exon GC
	
			lstCounts1 = getCol(lstMotifs, 2);
			lstCounts2 = getCol(lstMotifs, 5);

			print "correls = ", len(lstCounts1);
			
			if len(lstCounts1) >= _MIN_CORRELATION_COUNTS_: # we need at least this many pairs to compute correlation
				
				# Test 29/03/2017 -- use real average exon GC instead of middle of bin
				#meanExonGC = meanGC(getCol(lstMotifs, 0))
				#print "meanGC = ", meanExonGC;

				
				Rvalue = _safeCalcPearson(lstCounts1, lstCounts1, lstCounts2);
				print "Mean GC=", MeanExonGC, "Motif GC=", GC_Content(motif), "spacing=", spacing, "r2=", Rvalue;
				lstRawCorrelData.append([MeanExonGC, GC_Content(motif), spacing, Rvalue]);
				#lstRawCorrelData.append([meanExonGC, GC_Content(motif), spacing, Rvalue]);

	
	if bPrint:
		print lstRawCorrelData;
	
	return lstRawCorrelData;


def rawMotifData(motif, motifFile, bPrint):
# Return array of simple tuples for motif for Mean ExonGC, MotifGC

	lstRawMotifData = [];

	f = open(motifFile, 'r');
	motifs_raw = f.read().splitlines();
 
	for motif_line in motifs_raw:
		exon, pos, count = motif_line.split(',');
		try:
			lstRawMotifData.append([dictGC[exon], GC_Content(motif)]);
		except:
			(1 == 1);
	
	return lstRawMotifData;

	
def getCorrels(aCorrelsFile, exonGC, motifGC, spacing):
# Load data from raw correlations file
    lines = loadText(aCorrelsFile);
    lstData = [];
    for l in lines:
        _exonGC, _motifGC, _spacing, _r = l.split(',');
        _exonGC = float(_exonGC);
        _motifGC = float(_motifGC);
        _spacing = float(_spacing);
        _r = float(_r);
        lstData.append([_exonGC, _motifGC, _spacing, _r]);
    return filter(lambda x: (x[0] == exonGC) and (x[1] == motifGC) and (x[2] == spacing), lstData);
	
		
def _calcPearson(motif, motifFile, bPrint):
# Compute Pearson correlation co-efficients for given fourmer
# return an array of motif and correlation co-efficients: [motif, Rvalue10, Rvalue50, Rvalue100, Rvalue200];
# _MIN_CORRELATION_COUNTS_ is used to assert a minimum number of reads to compute the correlation, typically about 10.
 
 
	if (_PEARSON_TABLE_CACHE_):
		lstResult = getPearsonCachedMotif(fourmer);
		if bPrint:
			print lstResult;
		return lstResult;

	global FIRST_QUARTILE;
 		
	dictRval = {10:0, 50:0, 100:0, 200:0}
	for spacing in [10,50,100,200]:
	
		if spacing in [10,50]:
			tol = 2;
		else:
			tol = 4;

		dictSp = {spacing: []}; # for holding multiple exons Pearson Correlations for each spacing.
								# NB: these will be averaged for each spacing.

		# n-spaced n-mers by exon
		lstMotifs = _motifsBySpace(motifFile, spacing, tol);
		lstCounts1 = getCol(lstMotifs, 2);
		lstCounts2 = getCol(lstMotifs, 5);
			
		
		if len(lstCounts1) >= _MIN_CORRELATION_COUNTS_: # we need at least this many pairs to compute correlation
			dictRval[spacing] = _safeCalcPearson(lstCounts1, lstCounts1, lstCounts2);
			
	
	if (bPrint):
		if _READ_RELATIVE_:
			print motif + ", " + "{0:.3f}".format(FIRST_QUARTILE) + ", " + str(dictRval[10]);
		else:
			print motif + ", " + "{0:.3f}".format(FIRST_QUARTILE) + ", " + str(dictRval[10]) + ", " + str(dictRval[50]) + ", " + str(dictRval[100])  + ", " + str(dictRval[200]);
 
	return [motif, dictRval[10], dictRval[50], dictRval[100], dictRval[200]];
 
 
def _calcPearsonCalib(motif, motifFile, bPrint):
# Compute Pearson Correlation co-efficient for given motif
	global FIRST_QUARTILE;
 
	# 10-spaced n-mers on the same exon
	lstExonCounts = _motifsBySpace(motifFile, 10, 2);								   
	lstCounts1 = getCol(lstExonCounts, 2);
	lstCounts2 = getCol(lstExonCounts, 5);
	# Calculate Pearson Co-efficient
	print lstCounts1;
	print lstCounts2;
	print numpy.corrcoef(lstCounts1, lstCounts2);
	Rvalue10 = numpy.corrcoef(lstCounts1, lstCounts2)[0, 1];
	print Rvalue10;
 
	# 50-spaced n-mers on the same exon
	lstExonCounts = _motifsBySpace(motifFile, 50, 2);								   
	lstCounts1 = getCol(lstExonCounts, 2);
	lstCounts2 = getCol(lstExonCounts, 5);
	# Calculate Pearson Co-efficient
	Rvalue50 = numpy.corrcoef(lstCounts1, lstCounts2)[0, 1];
								 
	if (bPrint):
		print motif + ", " + "{0:.3f}".format(FIRST_QUARTILE) + ", " + str(Rvalue10) + ", " + str(Rvalue50);
 
	return [motif, Rvalue10, Rvalue50];
	 
		 
def calcPearson(bPrint):
# Compute Pearson Correlation co-efficient for ALL fourmers in batch (given by folder below)
	global FIRST_QUARTILE;
	 
	if (_READ_RELATIVE_):
		read_rel_file = "-rr";
	else:
		read_rel_file = "";
 
	lstResult = [];
	#tmp = 0; # DEVELOPMENT ONLY

	if (_PEARSON_TABLE_CACHE_):
		return getPearsonTableCache();

	 
	if (_NMER_MODE_ == 4):
		for a in bases:
			for b in bases:
				for c in bases:
					for d in bases:
						fourmer=a+b+c+d;
						#tmp += 1;
						#if (tmp > 5): # DEVELOPMENT ONLY
						#   break;  
						job_prefix = "herc-final-"
						_NMER_JOB_FULL_PATH_ = _NMER_JOB_PATH_ + fourmer + "/" + job_prefix + fourmer + read_rel_file + ".csv";
						FIRST_QUARTILE = calcFirstQuartile(_NMER_JOB_FULL_PATH_);								   
						lstResult.append( _calcPearson(fourmer, _NMER_JOB_FULL_PATH_, bPrint) );
		return lstResult;
											 
	if (_NMER_MODE_ == 5):
		for a in bases:
			for b in bases:
				for c in bases:
					for d in bases:
						for e in bases:
							fivemer=a+b+c+d+e;
								#tmp += 1;
								#if (tmp > 5): # DEVELOPMENT ONLY
								#   break;  
							job_prefix = "herc-final-"
							_NMER_JOB_FULL_PATH_ = _NMER_JOB_PATH_ + fivemer + "/" + job_prefix + fivemer + read_rel_file + ".csv";
							FIRST_QUARTILE = calcFirstQuartile(_NMER_JOB_FULL_PATH_);								   
							lstResult.append( _calcPearson(fivemer, _NMER_JOB_FULL_PATH_, bPrint) );
		return lstResult;
 
def calcPearsonCalib(bPrint):
	global FIRST_QUARTILE;
	motif = "CC";   
	calibFile = "/home/local/mxba001/Dropbox/PySpark/Results/herc-final-CALIB50NORM-CC.csv";
	FIRST_QUARTILE = calcFirstQuartile(calibFile);
	return _calcPearsonCalib(motif, calibFile, bPrint);
 
 
def doPearsonCalib():
	# Compute Pearson correlation co-efficient for given
	print "motif, 1st_quartile, R(10), R(50), R(100), R(200)";
	calcPearsonCalib(True);
 
def doPearson():
	if _READ_RELATIVE_:
		print "motif, 1st_quartile, R(10)";
	else:
		print "motif, 1st_quartile, R(10), R(50), R(100), R(200)";
	calcPearson(True);
 
def RemoveAllNaNMotifs(lstCorells):
# Remove Motifs where all spacing correlations are NaN due to insufficient data

	# Columns for readability
	COL_MOTIF = 0;
	COL_R10   = 1;
	COL_R50   = 2;
	COL_R100  = 3;
	COL_R200  = 4;
	COL_VAL   = 1;

	lstResult = [];
	for c in lstCorells:
		# Column and NaN comparison:
		bSameColVals = (c[COL_R10] == c[COL_R50]) and (c[COL_R50] == c[COL_R100]) and (c[COL_R100] == c[COL_R200]);
		bNanVal = (c[COL_R10] == c[COL_R10]);
		if not ( bSameColVals and bNanVal ):
			lstResult.append(c);
			
	return lstResult;

 
def doPearson2(outlier_count):
	#print "motif, 1st_quartile, R(10), R(50), R(100), R(200)";
	 
	# Columns for readability
	COL_MOTIF = 0;
	COL_R10   = 1;
	COL_R50   = 2;
	COL_R100  = 3;
	COL_R200  = 4;
	COL_VAL   = 1;
	 
	# Calculate Pearson correlation co-efficients for ALL 256 (4^4) different
	# motif fourmers and their counts produced by the Spark pipeline
	lstPearsonTable = calcPearson(False);
	
	random.shuffle(lstPearsonTable);
	

	#lstPearsonTable = RemoveAllNaNMotifs(lstPearsonTable);
	
	# Obtain minimum and maximum range of Pearson correlation co-efficients
	# for the counts of each exon motif spacing, do this for each motif
	 
	# Create lists of the min and max ranges for each motif in the form of [ [motif, R], [motif, R], [motif, R] ]
	# for each exon motif spacing
	 
	lstPearsonTableNF = lstPearsonTable; # ALL Nan filtered
	 
	# 10
	lstPearsonTable = filter(lambda x: x[COL_R10] == x[COL_R10], lstPearsonTableNF);
	tmp = sorted(lstPearsonTable, key=lambda x: x[COL_R10], reverse=True)[:outlier_count];
	lstR10max = map(lambda a,b:[a[COL_MOTIF],b[COL_R10]], tmp, tmp);
	#print lstR10max;   
	lstPearsonTable = filter(lambda x: x[COL_R10] == x[COL_R10], lstPearsonTableNF);
	tmp = sorted(lstPearsonTable, key=lambda x: x[COL_R10], reverse=False)[:outlier_count];
	lstR10min = map(lambda a,b:[a[COL_MOTIF],b[COL_R10]], tmp, tmp);
	#print lstR10min;
	 
	# 50
	lstPearsonTable = filter(lambda x: x[COL_R50] == x[COL_R50], lstPearsonTableNF);
	tmp = sorted(lstPearsonTable, key=lambda x: x[COL_R50], reverse=True)[:outlier_count];
	lstR50max = map(lambda a,b:[a[COL_MOTIF],b[COL_R50]], tmp, tmp);
	#print lstR50min;
	lstPearsonTable = filter(lambda x: x[COL_R50] == x[COL_R50], lstPearsonTableNF);
	tmp = sorted(lstPearsonTable, key=lambda x: x[COL_R50], reverse=False)[:outlier_count];
	lstR50min = map(lambda a,b:[a[COL_MOTIF],b[COL_R50]], tmp, tmp);
	#print lstR50min;
	 
	# 100
	lstPearsonTable = filter(lambda x: x[COL_R100] == x[COL_R100], lstPearsonTableNF);
	tmp = sorted(lstPearsonTable, key=lambda x: x[COL_R100], reverse=True)[:outlier_count];
	lstR100max = map(lambda a,b:[a[COL_MOTIF],b[COL_R100]], tmp, tmp);
	#print lstR100max;
	lstPearsonTable = filter(lambda x: x[COL_R100] == x[COL_R100], lstPearsonTableNF);
	tmp = sorted(lstPearsonTable, key=lambda x: x[COL_R100], reverse=False)[:outlier_count];
	lstR100min = map(lambda a,b:[a[COL_MOTIF],b[COL_R100]], tmp, tmp);
	#print lstR100min;
	 
	# 200
	lstPearsonTable = filter(lambda x: x[COL_R200] == x[COL_R200], lstPearsonTableNF);
	tmp = sorted(lstPearsonTable, key=lambda x: x[COL_R200], reverse=True)[:outlier_count];
	lstR200max = map(lambda a,b:[a[COL_MOTIF],b[COL_R200]], tmp, tmp);
	#print lstR200max;
	lstPearsonTable = filter(lambda x: x[COL_R200] == x[COL_R200], lstPearsonTableNF);
	tmp = sorted(lstPearsonTable, key=lambda x: x[COL_R200], reverse=False)[:outlier_count];
	lstR200min = map(lambda a,b:[a[COL_MOTIF],b[COL_R200]], tmp, tmp);
	#print lstR200min;
	 
	print "Lowest " + str(outlier_count) + " Pearson-correlation outliers and their motifs"
	print "motif=R(10), motif=R(50), motif=R(100), motif=R(200)"
	for x in range(0,outlier_count):
		print (
		  lstR10min[x][COL_MOTIF]  + "=" + prettyFloat(lstR10min[x][COL_VAL]) + ", "
		+ lstR50min[x][COL_MOTIF]  + "=" + prettyFloat(lstR50min[x][COL_VAL]) + ", "
		+ lstR100min[x][COL_MOTIF] + "=" + prettyFloat(lstR100min[x][COL_VAL]) + ", "
		+ lstR200min[x][COL_MOTIF] + "=" + prettyFloat(lstR200min[x][COL_VAL])
		);
	 
	print "Highest " + str(outlier_count) + " Pearson-correlation outliers and their motifs"
	print "motif=R(10), motif=R(50), motif=R(100), motif=R(200)"
	for x in range(0,outlier_count):
		print (
		  lstR10max[x][COL_MOTIF]  + "=" + prettyFloat(lstR10max[x][COL_VAL]) + ", "
		+ lstR50max[x][COL_MOTIF]  + "=" + prettyFloat(lstR50max[x][COL_VAL]) + ", "
		+ lstR100max[x][COL_MOTIF] + "=" + prettyFloat(lstR100max[x][COL_VAL]) + ", "
		+ lstR200max[x][COL_MOTIF] + "=" + prettyFloat(lstR200max[x][COL_VAL])
		);


def loadPearsonTableCache():
	global lstPearsonTblCache;
	if (_PEARSON_TABLE_CACHE_):
		# If not already loaded...
		if (lstPearsonTblCache == []):
			if os.path.isfile(_PEARSON_TABLE_CACHE_FILE_):
				# Load intermediate data from cache
				print "Loading intermediate data from pearson table cache:", _PEARSON_TABLE_CACHE_FILE_;
				data = loadText(_PEARSON_TABLE_CACHE_FILE_);
				data2 = svLinesToArray(data, ',');
				for i in data2:
					lstPearsonTblCache.append([i[0], float(i[1]), float(i[2]), float(i[3]), float(i[4]), float(i[5])]);

	

def getPearsonTableCache():
	lstResult = [];
	if (_PEARSON_TABLE_CACHE_):
		loadPearsonTableCache();
		# Load intermediate data from cache
		for i in lstPearsonTblCache:
				lstResult.append([i[0], i[2], i[3], i[4], i[5]]);
	return lstResult;
	
def getPearsonCachedMotif(motif):
	loadPearsonTableCache();
	for i in lstPearsonTblCache:
		if (i[0] == motif):
			return [i[0], float(i[2]), float(i[3]), float(i[4]), float(i[5])];
			
def getPearsonCachedMotifQuartile(motif):
	#time.sleep(2);
	loadPearsonTableCache();
	for i in lstPearsonTblCache:
		if (i[0] == motif):
			return i[1];


# Plotting code [requires matplotlib module]
 
def PlotMotifsBySpace(sorted_motifs_file, space_bp, tolerance_bp, motif, bLog):
# Return motif occurrences on the same exon within a given spacing
 
	lstResults = _motifsBySpace(sorted_motifs_file, space_bp, tolerance_bp);
	fig = plt.figure();
	ax = fig.add_subplot(111);
	 
	lstCounts1 = getCol(lstResults, 2); # x
	lstCounts2 = getCol(lstResults, 5); # y
 
	print lstCounts1;
 
	#area = np.pi * (15 * np.random.rand(N))**2  # 0 to 15 point radiuses
	ax.scatter(lstCounts1, lstCounts2, s=_PLOT_SCATTER_PT_SIZE_, alpha=0.5)
 
	ax.set_title(str(space_bp) + ' spaced ' + motif + ' motif pairs +/-' + str(tolerance_bp) + 'bp');
	ax.set_xlabel('1st motif count');
	ax.set_ylabel('2nd motif count');
 
	# Hugh ex
	if bLog:
		ax.set_xscale("log", nonposx='clip');
		ax.set_yscale("log", nonposx='clip');
 
	plt.show()
 
 
def PlotGC(data, title, plot_filename, bSave):
# Plot spread of Pearson correlation co-efficients for varying GC content and spacings
 
	fig = plt.figure();
	ax = fig.add_subplot(111);
	 
	lstGC_cats = getCol(data, 1); # x
	lstCorells = getCol(data, 0); # y

	#area = np.pi * (15 * np.random.rand(N))**2  # 0 to 15 point radiuses
	ax.scatter(lstGC_cats, lstCorells, s=_PLOT_SCATTER_PT_SIZE_, alpha=0.5) 

	ax.set_title(title);
	ax.set_xlabel('correlation');   
	ax.set_ylabel('mean read GC content (exon)'); 
	major_ticks = np.arange(20, 80, 10);
	ax.set_yticklabels(["", "30-40%","40-50%","50-60%","60-70%", ""]);
	ax.set_yticks(major_ticks);

	# x-axis range, same for all plots for comparison purposes
	ax.set_xlim([-0.2,1.0]);
	 
	if (bSave):
		plt.savefig(_PLOT_OUT_PATH_ + plot_filename);

	#plt.show();
 

def PlotGCBox(data30_40, data40_50, data50_60, data60_70, title, plot_filename, bSave):
# Plot Box and Whisker of spread of Pearson correlation co-efficients for varying
# GC content and spacings

	fig = plt.figure();
	ax = fig.add_subplot(111);

	data = [data30_40, data40_50, data50_60, data60_70];

	ax.boxplot(data, 0, '.', 0);

	ax.set_title(title);
	ax.set_xlabel('correlation');   
	ax.set_ylabel('mean read GC content (exon)'); 
	ax.set_yticklabels(["30-40%","40-50%","50-60%","60-70%"]);

	# x-axis range, same for all plots for comparison purposes
	ax.set_xlim([-0.2,1.0]);
	 
	if (bSave):
		plt.savefig(_PLOT_OUT_PATH_ + plot_filename);

	#plt.show();

	

def saveRawCorrelationsData(lstRawCorrel):
# Save simeple table of variablesand correlations data as a csv file
	
	lines = [];	
	for i in lstRawCorrel:
		lineStr = ','.join(map(str, i));
		lines.append(lineStr);
		
	FileName = _CORREL_FILE_;
	FileName = FileName.replace(".txt", ".csv");	
	FileName = FileName.replace("Motif-GC-cube", "Correls");
	saveText(FileName, lines);	
	return FileName;

 
def Histogram(vals, plot_output_file, hist_title):
# Plot histogram of frequencies of given values
# NB: plot_output_file should be full path, ending in
# .png as that is the output format
 
	fig = plt.figure(figsize=(16,12)) # Size in inches
	ax = fig.add_subplot(111)
 
	numBins = len(set(vals)) / 8;
	#print numBins;
 
	ax.hist(vals, numBins,color='green',alpha=0.8)
 
	ax.set_title(hist_title);
	ax.set_xlabel('motif count');
	ax.set_ylabel('frequency');
 
	plt.savefig(plot_output_file);
	#plt.show()
	 
 
 
def Histogram4mer(fourmer_file, plot_output_path, fourmer):
# Plot histogram of frequencies of counts for given motif
	f = open(fourmer_file, 'r');
	motifs_raw = f.read().splitlines();
 
	lstMotifs = [];
 
	for motif_line in motifs_raw:
		exon, pos, count = motif_line.split(',');
		try:
			lstMotifs.append([exon, int(pos), int(count)]);
		except:
			(1 == 1);
			#print "error: " + motif_line;
 
	filtered_group = removeNoise(lstMotifs);
	lstCounts = getCol(filtered_group, 2);
	#print lstCounts;
	Histogram(lstCounts, plot_output_path, 'motif count distribution for  ' + fourmer);
	 
 
def PlotMotifCorrelation4mer(fourmer_file, plot_output_path, fourmer):
# Plot histogram of frequencies of counts for given motif
	f = open(fourmer_file, 'r');
	motifs_raw = f.read().splitlines();
 
	lstMotifs = [];
 
	for motif_line in motifs_raw:
		exon, pos, count = motif_line.split(',');
		try:
			lstMotifs.append([exon, int(pos), int(count)]);
		except:
			(1 == 1);
			#print "error: " + motif_line;
 
	filtered_group = removeNoise(lstMotifs);
	lstCounts = getCol(filtered_group, 2);
	#print lstCounts;
	Histogram(lstCounts, plot_output_path, 'motif count distribution for  ' + fourmer);
 
 
def PlotScatterMatrix(title, CorrelsFileName):
	
	# Load dataset
	names = ['exon-gc', 'motif-gc', 'spacing', 'r'];
	dataset = pandas.read_csv(CorrelsFileName, names=names);

	scatter_matrix(dataset);

	PlotFileName = CorrelsFileName.replace(".csv", ".png");
	plt.savefig(PlotFileName);


 
# Latex generating code
 
def latexDocumentHeader(aTitle):
	global lstLatex;
	lstLatex.append("\\documentclass[10pt,a4paper]{article}");
	lstLatex.append("\\usepackage[utf8]{inputenc}");
	lstLatex.append("\\usepackage{amsmath}");
	lstLatex.append("\\usepackage{amsfonts}");
	lstLatex.append("\\usepackage{amssymb}");
	lstLatex.append("\\usepackage[authoryear, round]{natbib}");
	lstLatex.append("\\usepackage{graphicx}");
	lstLatex.append("\\usepackage{subfigure}");
	lstLatex.append("\\usepackage{array}");
	lstLatex.append("\\usepackage{arydshln}");
	lstLatex.append("\\newcolumntype{?}{!{\\vrule width 1pt}}");
	lstLatex.append("\\author{Jamie Alnasir}");
	lstLatex.append("\\title{" + aTitle + "}");
	lstLatex.append("\\begin{document}");
	lstLatex.append("");
	lstLatex.append(latexBold(_JOB_FRIENDLY_TITLE_));
		
	
def latexDocumentFooter():	
	global lstLatex;
	lstLatex.append("\end{document}");

def latexStart():
	global lstLatex;
	latexDocumentHeader("Test document");
	
def latexPrint(aText):
	global lstLatex;
	lstLatex.append(aText);
	
def latexMotifSpacingTableStart(aSpacing):
	t = "\\begin{table}[h]\n\\scalebox{0.8}{\n\\begin{tabular}{?c?p{1.5cm}|p{1.5cm}?p{1.5cm}|p{1.5cm}?p{1.5cm}|p{1.5cm}?p{1.5cm}|p{1.5cm}?}\n\\hline";
	latexPrint(t);
	latexPrint("\multicolumn{9}{|c|}{ Motif spacing: " + latexBold(aSpacing) + "}\\\\");
	latexPrint("\\hline");
	
def latexMotifSpacingTable(dictExonGC):
	
	lstLabels = ["Exon GC\%", "-", "Motif GC\%", "\#Correlations", "p(t-test)", "p(Wilcoxon)", "-"];
	lstLabels = lstLabels + lstLabels + lstLabels + lstLabels;
	
	#latexPrint(str(dictExonGC));
	
	for l, gc3040, gc4050, gc5060, gc6070 in zip(lstLabels, dictExonGC[30], dictExonGC[40], dictExonGC[50], dictExonGC[60]):
		
		if (l == "-"):
			latexPrint("\\hline");
		else:
			latexPrint(l + " & " + gc3040 + " & " + gc4050 + " & " + gc5060 + " & " + gc6070 + " \\\\");
		
		if (l == "\#Correlations"):
			latexPrint("\hdashline[2pt/1.75pt]");
	

def latexFloat(n):
	if (n!=n):
		return "NaN";
	
	floatStr =  "{:.2E}".format(n);
	
	if floatStr.find("E-0") <> -1:	# search str must match replace in next line
		floatStr = floatStr.replace("E-0", "x10$^{-") + "}$";
	
	floatStr = floatStr.replace("E+00", "");
	return floatStr;
	
def latexBold(atxt):
	return "\\textbf{" + atxt + "}";
	
def latexMotifSpacingTableEnd(species, spacing):
	species = species.replace("[", "");
	species = species.replace("]", "");
	c = "T-test and Wilcoxon-test comparisons of Pearson correlations for motif-pairs at " + spacing + " spacing for varying motif GC and mean exon GC content in " + species + ". FDR corrected p-values in parenthesis. * suggests rejection of null hypothesis.";
	c2= "Pearson correlations comparisons for " + species + ", " + spacing + ".";
	t = """\\end{tabular}}
\caption[""" + c2 + """]{\label{tab:PearsonComp""" + _SPECIES_LABEL_ + spacing + """}""" + c +  """ }
\\end{table}\n""";

	latexPrint(t);	

def latexReplace(src, dst):
	global lstLatex;
	lstTmp = [];
	for l in lstLatex:
		l = l.replace(src, dst);
		lstTmp.append(l);
	lstLatex = lstTmp;


def latexEnd():
	global lstLatex;		
	latexDocumentFooter();
	saveText(_CORREL_FILE_.replace("./", "./Latex.tmp/") + ".tex", lstLatex);




 
if (_PEARSON_TABLE_CACHE_):
		loadPearsonTableCache();
		
if (len(sys.argv) > 2):
	global FIRST_QUARTILE;
	global _GC_FILTERING_;
	global _GC_FILTER_MIN_;
	global _GC_FILTER_MAX_;
	global _MOTIF_GC_CUBE_FILE_;
	 
	_INPUT_DATA_ = sys.argv[1]; # Pick up input data
	 
	# Load exon GC content lookup data
	LoadGC(_GC_FILE_PATH_);
	 
	FIRST_QUARTILE = calcFirstQuartile(_INPUT_DATA_);   
	#print "First quartile = ", FIRST_QUARTILE;
 
	if (sys.argv[2] == "FIRST-XY"):
		print "x, y";
		selectMotifSignals(_INPUT_DATA_);
	elif (sys.argv[2] == "POS"):
		print "exon, pos, count";
		selectMotifSignalsPosition(_INPUT_DATA_, long(sys.argv[3]), long(sys.argv[4]));
	elif (sys.argv[2] == "DUAL-POS"):
		print "exon, pos1, count1, exon, pos2, count2";
		selectMotifSignalsPositions(_INPUT_DATA_, long(sys.argv[3]), long(sys.argv[4]), long(sys.argv[5]));
	elif (sys.argv[2] == "SPACED-BY"):
		print "exon, pos1, count1, exon, pos2, count2";
		motifsBySpace(_INPUT_DATA_, long(sys.argv[3]), long(sys.argv[4]));
	elif (sys.argv[2] == "PLOT-SPACED-BY"):
		if len(sys.argv) > 6:
			bLog = (sys.argv[6].upper() == "LOG");
		else:
			bLog = False;
		PlotMotifsBySpace(_INPUT_DATA_, long(sys.argv[3]), long(sys.argv[4]), sys.argv[5], bLog);
	elif (sys.argv[2] == "PEARSON-SPACED-BY"):
		print "motif, 1st_quartile, R(10), R(50), R(100), R(200)";
		_calcPearson("CALC:", _INPUT_DATA_, True);
	else:
		print "Invalid mode specified";
		print "Usage Hercules-graphit.py <herc-final-data.csv> <FIRST-XY|POS|SPACED-BY> [*pos1 pos2 tolerance] *if POS mode";
		 
elif (len(sys.argv) == 2):

	# Load exon GC content lookup data
	LoadGC(_GC_FILE_PATH_);
	 
	print "-- " + _JOB_FRIENDLY_TITLE_ + " --";

	if (_READ_RELATIVE_):
			print "(** Read-Read-Relative data **)";

	if (_GC_FILTERING_):
			print "(Exons with mean % GC content between " + str(_GC_FILTER_MIN_) + "-" + str(_GC_FILTER_MAX_) + "%)";
			 
	if (sys.argv[1] == "PEARSON-CALIBRATION"):
		 
		doPearsonCalib();
	 
	elif (sys.argv[1] == "PEARSON-TABLE1"):	 
		 
		doPearson();
		 
	elif (sys.argv[1] == "PEARSON-TABLE2"):
			 
		outlier_count = 10;
		doPearson2(outlier_count);
		 
	elif (sys.argv[1] == "PEARSON-TABLE2-GC"):
		 
		outlier_count = 10;
		doPearson2(outlier_count);

	if (sys.argv[1] == "DEBUG"):
		print "DEBUG MODE:";
		#print GC_Lookup("00000002400000001265");
		print lstPearsonTblCache;
		print getPearsonCachedMotif("AAAA");
		print getPearsonCachedMotifQuartile("AAAA");

	elif (sys.argv[1] == "4MER-COUNT-HISTO"):

		## 4mers for Wild-type Drosophila
		## _4MER_JOB_PATH_ = "/home/local/mxba001/PhD/_working/4mer/";

		## 4mers for Mutant-r2-type Drosophila
		#_4MER_JOB_PATH_ = "/home/local/mxba001/PhD/_working/4mer/4mer-mut-r2-/";

		## 4mers for Mutant-r2-type Drosophila
		_4MER_JOB_PATH_ = "/home/local/mxba001/PhD/_working/4mer-mut-r2-run2/";


		_4MER_HIST_OUT_PATH_ = "/home/local/mxba001/Dropbox/PySpark/Results/Drosophila/wild-type/count-histograms/";

		Bases = "ATGC";
		lstResult = [];
		for a in Bases:
			for b in Bases:
				for c in Bases:
					for d in Bases:
						fourmer=a+b+c+d;

						FIRST_QUARTILE = calcFirstQuartile(_4MER_JOB_PATH_ + fourmer + '/herc-final-' + fourmer +
'.csv');
						FIRST_QUARTILE = 20;
						_4MER_HIST_OUTFILE_ = _4MER_HIST_OUT_PATH_ + fourmer + '.png';
						print "Outputting Histogram: ", fourmer, FIRST_QUARTILE;
						Histogram4mer(_4MER_JOB_PATH_ + fourmer + '/herc-final-' + fourmer + '.csv', _4MER_HIST_OUTFILE_, fourmer);

	elif (sys.argv[1] == "PLOT-4MER-SCATTERMATRIX"):
		print "Plotting 4mer scatter for various GC parameters and spacings";
		_GC_FILTERING_ = False;

		lstRawCorrels = [];
		Bases = "ATGC";
		lstResult = [];
		for a in Bases:
			for b in Bases:
				for c in Bases:
					for d in Bases:
						fourmer=a+b+c+d;
						print fourmer;
						job_prefix = "herc-final-";
						_NMER_JOB_FULL_PATH_ = _NMER_JOB_PATH_ + fourmer + "/" + job_prefix + fourmer + ".csv";

						FIRST_QUARTILE = calcFirstQuartile(_NMER_JOB_FULL_PATH_);
						print FIRST_QUARTILE;
						print _NMER_JOB_FULL_PATH_;
						
						# Obtain set of tuple results, 4 tuples, 1 for each spacing
						lstMotifCorrels = rawCorrelationsData(fourmer, _NMER_JOB_FULL_PATH_, False);
						#for m in lstMotifCorrels:
						#	print m;
						
						# Motif GC test
						#lstMotifCorrels = rawMotifData(fourmer, _NMER_JOB_FULL_PATH_, False);
						
						lstRawCorrels.extend(lstMotifCorrels);
		
		#for c in lstRawCorrels:
		#	print c;
		
		FileName = saveRawCorrelationsData(lstRawCorrels);
		print "Raw correlations: ", len(lstRawCorrels);
		print "Raw correlation records saved to: ", FileName;
		

						
	elif (sys.argv[1] == "PLOT-4MER-GC-CONTENT"):
		print "Plotting 4mer GC content for various motif-GC content and spacings";

		# We'll be cycling the filter range for the computations
		_GC_FILTERING_=True;

		# GC content bins
		lstREAD_GC_bins  = [[30,40], [40,50], [50,60], [60,70]];
		 
		# Store results
		dictResults = {30:{0:[], 1:[], 2:[], 3:[], 4:[]},  # Store [MOTIF, 10,50,100,200bp correlations], grouped by %GC content of motif. indexed by Read GC
					   40:{0:[], 1:[], 2:[], 3:[], 4:[]},  # Store [MOTIF, 10,50,100,200bp correlations], grouped by %GC content of motif. indexed by Read GC
					   50:{0:[], 1:[], 2:[], 3:[], 4:[]},  # Store [MOTIF, 10,50,100,200bp correlations], grouped by %GC content of motif. indexed by Read GC
					   60:{0:[], 1:[], 2:[], 3:[], 4:[]}}; # Store [MOTIF, 10,50,100,200bp correlations], grouped by %GC content of motif. indexed by Read GC

		 
		if os.path.isfile(_MOTIF_GC_CUBE_FILE_):
			# Load intermediate data from the cube
			print "Loading intermediate data from cube:", _MOTIF_GC_CUBE_FILE_;
			dictData = LoadData(_MOTIF_GC_CUBE_FILE_);
			# Fix dictionary key, in saving simplejson converts it to a string, we want an int key
			dictData = {int(k):dict(v) for k,v in dictData.items()};
			for gc_range in lstREAD_GC_bins:
				dictData[gc_range[0]] = {int(k):list(v) for k,v in dictData[gc_range[0]].items()};
			dictResults = dictData;
			time.sleep(2);
		else:
			print "Computing data for the cube:", _MOTIF_GC_CUBE_FILE_;
			for a in bases:
				for b in bases:
					for c in bases:
						for d in bases:

							fourmer=a+b+c+d;
							print fourmer;
							fourmer_gc = GC_Content(fourmer);
							bin_index_motif_gc = int(fourmer_gc / 25);
							print fourmer_gc, bin_index_motif_gc;
							#dictMOTIF_GC_bins[bin_index_motif_gc].append(fourmer);

							job_prefix = "herc-final-";
							read_rel_file = "";
							_NMER_JOB_FULL_PATH_ = _NMER_JOB_PATH_ + fourmer + "/" + job_prefix + fourmer + read_rel_file + ".csv";
							
							if (_PEARSON_TABLE_CACHE_):
								FIRST_QUARTILE = getPearsonCachedMotifQuartile(fourmer);
							else:
								FIRST_QUARTILE = calcFirstQuartile(_NMER_JOB_FULL_PATH_);
							
							for gc_range in lstREAD_GC_bins:
								_GC_FILTER_MIN_ = gc_range[0];
								_GC_FILTER_MAX_ = gc_range[1];
								bin_index_exon_gc = gc_range[0];
								 
								print "Motif GC content bin:", bin_index_motif_gc , " Median Exon GC content range: ", gc_range[0], '-', gc_range[1];
								dictResults[bin_index_exon_gc][bin_index_motif_gc].append( _calcPearson(fourmer, _NMER_JOB_FULL_PATH_, True) );
							 
			SaveData(_MOTIF_GC_CUBE_FILE_, dictResults);
		
		#print "--------------------------------";
		#print dictResults[30];
		#print "--------------------------------";
		#print dictResults[40];
		#print dictResults[50];
		#print dictResults[60];

		lstPts = []; # Points for the plots
		lstSpacings = ["dummy", "10bp", "50bp", "100bp", "200bp"];
		for bin_index_motif_gc in range(0,4+1):
			motif_gcp = bin_index_motif_gc * 25;
			for spacing_index in range(1,4+1):		  
			 
				# Start new plot
				lstPts = [];
			 
				for exon_gc in lstREAD_GC_bins:
					exon_gc_lbl = str(exon_gc[0]) + "-" + str(exon_gc[1]);
					bin_index_exon_gc = exon_gc[0];
					 
					print "motif GC:", motif_gcp,"spacing:", spacing_index, " exon GC:", bin_index_exon_gc;
					lstCorells = dictResults[bin_index_exon_gc][bin_index_motif_gc];
					#print lstCorells;
					for corel in lstCorells:
						lstPts.append([exon_gc[0], corel[spacing_index]]);
				 
				# Do plotting here
				#print lstPts;
				#print getCol(lstPts, 0);
				plotTitle = "motif GC content " + str(motif_gcp) + "%, " + str(lstSpacings[spacing_index]) + " spaced motif pairs";
				plotFileName = "motif-gc-" + str(motif_gcp) + "pc-" + str(lstSpacings[spacing_index]) + "-spaced.png";
				PlotGC(lstPts, plotTitle, plotFileName, True);	  

	elif (sys.argv[1] == "PLOT-4MER-GC-BOX"):
		print "Plotting Box and Whisker of 4mer GC content for various motif-GC content and spacings";

		# We'll be cycling the filter range for the computations
		_GC_FILTERING_=True;

		_MOTIF_GC_COMPARE_WITH_GC_ = 50;

		# GC content bins
		lstREAD_GC_bins  = [[30,40], [40,50], [50,60], [60,70]];
		
		
		# hold plot data for exonGC
		dictPts = {35:[], 45:[], 55:[], 65:[]};
		
		pCounter = -1; # To ensure a 0 start
		lstSpacings = [10, 50, 100, 200];
		for spacing in lstSpacings:
			spacingStr = str(spacing) + "bp"
			print '\n' + '\n' + "[ motif spacing: ", spacingStr+ " ]";
			
			dictLatexReadGCblocks = {30:[], 40:[], 50:[], 60:[]};
			
			latexMotifSpacingTableStart(spacingStr);
			
			for bin_index_motif_gc in range(0,4+1):
				motifGC = bin_index_motif_gc * 25;
				
				for exon_gc in lstREAD_GC_bins:
					exonGC = float((exon_gc[1] - exon_gc[0]) / 2) + exon_gc[0];
					exon_gc_lbl = str(exon_gc[0]) + "-" + str(exon_gc[1]);
					bin_index_exon_gc = exon_gc[0];
					#print '\n' + "Exon median read GC content: ", exon_gc_lbl;
					
					# load correlations for given exonGC, motifGC, spacing
					lstA = getCorrels(_CORREL_FILE_, exonGC, motifGC, spacing);
					dictPts[int(exonGC)] = getCol(lstA, 3); # get correlations col
					
					print spacing, exon_gc_lbl, exonGC, motifGC, ":", len(lstA), "correlations";

				plotTitle = "motif GC content " + str(motifGC) + "%, " + str(spacingStr) + " spaced motif pairs";
				plotFileName = "box-motif-gc-" + str(motifGC) + "pc-" + str(spacingStr) + "-spaced.png";
				PlotGCBox(dictPts[35], dictPts[45], dictPts[55], dictPts[65], plotTitle, plotFileName, True);	
				print "Written: ", plotFileName;		
		 
	elif (sys.argv[1] == "4MER-GC-CONTENT-STATS"):
		print "Plotting 4mer GC content for various motif-GC content and spacings";
		
		# latexXYZ commands generate a Latex document of the results.
		latexStart();

		# We'll be cycling the filter range for the computations
		_GC_FILTERING_=True;

		_MOTIF_GC_COMPARE_WITH_GC_ = 50;

		# GC content bins
		lstREAD_GC_bins  = [[30,40], [40,50], [50,60], [60,70]];
 
		lstPts = []; # Points for the plots
		lstPvalsT = []; # T-Test p-values; 
		lstPvalsW = []; # W-Test p-values; 
		pCounter = -1; # To ensure a 0 start
		lstSpacings = [10, 50, 100, 200];
		for spacing in lstSpacings:
			spacingStr = str(spacing) + "bp"
			print '\n' + '\n' + "[ motif spacing: ", spacingStr+ " ]";
			
			dictLatexReadGCblocks = {30:[], 40:[], 50:[], 60:[]};
			
			latexMotifSpacingTableStart(spacingStr);
 
			for exon_gc in lstREAD_GC_bins:
				exonGC = float((exon_gc[1] - exon_gc[0]) / 2) + exon_gc[0];
				exon_gc_lbl = str(exon_gc[0]) + "-" + str(exon_gc[1]);
				bin_index_exon_gc = exon_gc[0];
				print '\n' + "Exon median read GC content: ", exon_gc_lbl;

				lstA = [];
				lstB = [];
				
				for bin_index_motif_gc in range(0,4+1):
					motifGC = bin_index_motif_gc * 25;
				
					lstA = getCol(getCorrels(_CORREL_FILE_, exonGC, motifGC, spacing), 3);
					#print spacing, exon_gc_lbl, exonGC, motifGC, ":", len(lstA), "correlations (list a)";
					
					lstB = getCol(getCorrels(_CORREL_FILE_, exonGC, _MOTIF_GC_COMPARE_WITH_GC_, spacing), 3);
					#print spacing, exon_gc_lbl, exonGC, motifGC, ":", len(lstB), "correlations (compared to, b)";									
					
					if (motifGC == _MOTIF_GC_COMPARE_WITH_GC_):
						continue; # skip self-self comparison
					
						
					print spacing, exon_gc_lbl, exonGC, _MOTIF_GC_COMPARE_WITH_GC_, "vs", motifGC;
					print " " + str(_MOTIF_GC_COMPARE_WITH_GC_) + "%" + " vector:", len(lstA), ", " + str(motifGC) + "% vector:", len(lstB), "correlations";

					
					if (len(lstA) > 2) and (len(lstB) > 2):

						pCounter += 1; # index which pair of p-values we're computing, for FDR correction

						t, p = stats.ttest_ind(lstA, lstB, equal_var=False);					
							 
						if p == float('NaN'):
							print "TEST!";

						print t, p;
						lstPvalsT.append(p); # Store T-test p value for FDR correction
						

						dictLatexReadGCblocks[exon_gc[0]].append( "\\multicolumn{2}{c?}{\\textbf{" + exon_gc_lbl + "\%}" + "}");
						dictLatexReadGCblocks[exon_gc[0]].append("--");
						dictLatexReadGCblocks[exon_gc[0]].append( str(_MOTIF_GC_COMPARE_WITH_GC_) +" & " + str(motifGC) );
						dictLatexReadGCblocks[exon_gc[0]].append( str(len(lstA)) +" & " + str(len(lstB)) );
						 
						strPcorrection = "(" + str(pCounter).zfill(3) + "T)";
						dictLatexReadGCblocks[exon_gc[0]].append("\\multicolumn{2}{p{3cm}?}{\small{" + latexFloat(p) + strPcorrection +  "}}"); # p(t)
						 
						print " t-test, t = %g  p = %g" % (t, p);
						lstA, lstB = sliceLargerBySmaller(lstA, lstB);
						w, p = stats.wilcoxon(lstA, lstB);
						print " wilcoxon pairs:", len(lstA);
						print " wilcoxon, w = %g  p = %g" % (w, p);
						 
						lstPvalsW.append(p); # Store T-test p value for FDR correction
						 
						strPcorrection = "(" + str(pCounter).zfill(3) + "W)";
						dictLatexReadGCblocks[exon_gc[0]].append("\\multicolumn{2}{p{3cm}?}{\small{" + latexFloat(p) + strPcorrection + "}}"); # p(w)
						dictLatexReadGCblocks[exon_gc[0]].append("--");


					else:
						print "Insufficient data!";
						dictLatexReadGCblocks[exon_gc[0]].append( "\\multicolumn{2}{c?}{\\textbf{" + exon_gc_lbl + "\%}" + "}");
						dictLatexReadGCblocks[exon_gc[0]].append("--");
						dictLatexReadGCblocks[exon_gc[0]].append( str(_MOTIF_GC_COMPARE_WITH_GC_) +" & " + str(motifGC) );
						dictLatexReadGCblocks[exon_gc[0]].append( str(len(lstA)) +" & " + str(len(lstB)) );
						 
						dictLatexReadGCblocks[exon_gc[0]].append("\\multicolumn{2}{c?}{insufficient}");
						dictLatexReadGCblocks[exon_gc[0]].append("\\multicolumn{2}{c?}{data}");
						dictLatexReadGCblocks[exon_gc[0]].append("--");
			
			
			latexMotifSpacingTable(dictLatexReadGCblocks);
			latexMotifSpacingTableEnd(_JOB_FRIENDLY_TITLE_, spacingStr);
			
			lstPvalsTCorrected = fdrcorrection(lstPvalsT);
			lstPvalsWCorrected = fdrcorrection(lstPvalsW);
			
			rstr = "*";
			for i in range(0, len(lstPvalsT)):
				bReject = lstPvalsTCorrected[0][i];
				if (bReject):
					rstr = "*";
				else:
					rstr = "";
				latexReplace( "(" + str(i).zfill(3) + "T)", "(" + latexFloat(lstPvalsTCorrected[1][i]) + ")" + rstr );
			 
			rstr = "*";
			for i in range(0, len(lstPvalsW)):
				bReject = lstPvalsWCorrected[0][i];
				if (bReject):
					rstr = "*";
				else:
					rstr = "";
				latexReplace( "(" + str(i).zfill(3) + "W)", "(" + latexFloat(lstPvalsWCorrected[1][i]) + ")" + rstr );

		
		# Save FDR corrected p-values for plotting density plot
		FileName = _CORREL_FILE_;
		FileName = FileName.replace(".csv", ".pt");
		saveNumArray(FileName, lstPvalsTCorrected[1]);
		FileName = _CORREL_FILE_;
		FileName = FileName.replace(".csv", ".pw");
		saveNumArray(FileName, lstPvalsWCorrected[1]);
		 
		latexEnd();
		
	elif (sys.argv[1] == "PLOT-MOTIF"):
		_4MER_JOB_PATH_ = sys.argv[2];
		_MOTIF_ = sys.argv[3];


else:
	print "Usage Hercules-graphit.py <herc-final-data.csv> <OPERATION>";

