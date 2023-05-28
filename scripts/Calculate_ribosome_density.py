#!/bin/python
import argparse
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def getTranscriptBasicInformation(mrna_lncrna_fa):
	transCategory,transLen,utr3Len,utr5Len,cdsLen,transORF = {},{},{},{},{},{}
	trans_fh = open(mrna_lncrna_fa) 
	for line in trans_fh:
		if line[0]==">":
			a = line.rstrip().split()
			tID = a[0][1:]
			tCat = a[2].split("=")[1]
			transCategory[tID] = tCat
			if tCat == "mRNA":
				b = a[3].split("=")[1].split("-")
				transORF[tID] = [int(b[0]), int(b[1])+3] 
				cdsLen[tID] = int(b[1]) - int(b[0]) + 3
				utr5Len[tID] = int(b[0]) 
		else:
			line = line.rstrip()
			transLen[tID] = len(line) 
			if tCat == "mRNA":
				utr3Len[tID] = transLen[tID] - transORF[tID][1] 
	trans_fh.close()
	return transCategory,transLen,transORF,utr5Len,utr3Len,cdsLen


def getTranscriptExpression(rsem_isoforms_results):
	fh = open(rsem_isoforms_results)
	fh.readline()
	transRPKM = {}
	for line in fh:
		a = line.split()
		tid, fpkm = a[0], float(a[6])
		transRPKM[tid] = fpkm
	fh.close()
	return transRPKM


def check_read(read):
	read_len = len(read)
	nt = {"A":0, "C":0, "G":0, "T":0, "N":0}
	for i in range(0,read_len):
		if read[i] in nt:
			nt[read[i]] += 1
	flag = True
	for k in nt:
		if nt[k]/float(read_len) > 0.85:
			flag = False
	return flag


def batchRiboseqAlignmentByRead(riboseq_sam,transRPKM,transLen):
	rid = None
	batch = []
	riboProfile = {}
	for tid in transRPKM:
		riboProfile[tid] = [0.]*transLen[tid]
	sam_fh = open(riboseq_sam)
	for line in sam_fh:
		if line[0]=="@": continue
		a = line.split()
		if rid != a[0]:
			countAlignmentFromBatch(batch,transRPKM,transLen,riboProfile)
			batch[:] = []
		rid = a[0]
		if check_read(a[9]): batch.append(line)
	countAlignmentFromBatch(batch,transRPKM,transLen,riboProfile)
	sam_fh.close()
	return riboProfile


def countAlignmentFromBatch(batch,transRPKM,transLen,riboProfile):
	batch_size = len(batch)
	if batch_size == 0: return
	total_rpkm = 0.
	for line in batch:
		a = line.split()
		total_rpkm += transRPKM[a[2]]

	for line in batch:
		a = line.split()
		tid = a[2]
		w = 1./batch_size
		if total_rpkm > 0: w = transRPKM[tid]/total_rpkm
		align_len = int(a[5][:-1])
		align_pos = int(a[3])
		riboProfile[tid][align_pos] += w


def calculateRiboDensity(transRPKM,transLen,transORF,transCategory,riboProfile,prefix):
	cdsDensity,lncrnaDensity,utr3Density = {},{},{}
	total_reads = 0.
	for tid in transRPKM:
		total_reads += sum(riboProfile[tid])
		if transRPKM[tid] < 1: continue  #Remove low-expressed transcripts
		if transCategory[tid] == "mRNA":
			cdsDensity[tid] = sum(riboProfile[tid][transORF[tid][0]:transORF[tid][1]])/float(transORF[tid][1]-transORF[tid][0])
			if transLen[tid]-transORF[tid][1] >= 100: #3UTR with length of >= 100bp will be considered as a negative control
				utr3Density[tid] = sum(riboProfile[tid][transORF[tid][1]:transLen[tid]])/float(transLen[tid]-transORF[tid][1])
		else:
			lncrnaDensity[tid] = sum(riboProfile[tid])/float(transLen[tid])
	scale_factor = float(10**9)/total_reads
	pseudo_number = 10.**-10
	for tid in cdsDensity:
		cdsDensity[tid] = math.log(cdsDensity[tid]*scale_factor/transRPKM[tid]+pseudo_number,2)
	for tid in utr3Density:
		utr3Density[tid] = math.log(utr3Density[tid]*scale_factor/transRPKM[tid]+pseudo_number,2)
	for tid in lncrnaDensity:
		lncrnaDensity[tid] = math.log(lncrnaDensity[tid]*scale_factor/transRPKM[tid]+pseudo_number,2) 

	utr3_cutoff = np.percentile(utr3Density.values(), 90)
	fh = open(prefix+"/riboDensity.csv", "w")
	fh.write("#3utr_ribosome_density_90percentile: %f\n"%utr3_cutoff)
	fh.write("#transcript_id,region,ribosome_density\n")
	for tid in cdsDensity:
		fh.write(tid+",CDS,"+str(cdsDensity[tid])+"\n")
	for tid in utr3Density:
		fh.write(tid+",3UTR,"+str(utr3Density[tid])+"\n")
	for tid in lncrnaDensity:
		fh.write(tid+",lncRNA,"+str(lncrnaDensity[tid])+"\n")
	fh.close()

	return cdsDensity,lncrnaDensity,utr3Density


def getRiboDensity(prefix):
	cdsDensity,lncrnaDensity,utr3Density = {},{},{}
	fh = open(prefix+"/riboDensity.csv")
	for line in fh:
		if line[0]=="#": continue
		a = line.rstrip().split(",")
		if a[1] == "CDS":
			cdsDensity[a[0]] = float(a[2])
		elif a[1] == "3UTR":
			utr3Density[a[0]] = float(a[2])
		else:
			lncrnaDensity[a[0]] = float(a[2])
	fh.close()
	return cdsDensity,lncrnaDensity,utr3Density


def generatePlotOfRiboDensity(cdsDensity,lncrnaDensity,utr3Density,prefix): 
	utr3_cutoff = np.percentile(utr3Density.values(), 90)
	plt.figure(figsize=(8, 4), dpi=150)
	sns.kdeplot(utr3Density.values(), color="black", label="3\' UTR(n=%d)"%(len(utr3Density)), shade=True)
	sns.kdeplot(cdsDensity.values(), color="green", label="CDS(n=%d)"%(len(cdsDensity)), shade=True)
	sns.kdeplot(lncrnaDensity.values(), color="red", label="lncRNA(n=%d)"%(len(lncrnaDensity)), shade=True)
	plt.axvline(utr3_cutoff, c='black', alpha=1.0, linestyle='dashed', linewidth = 0.8)
	plt.xlim(-42, 10)
	plt.xlabel("Log2(ribosome density + 10e-10)", size=16)
	plt.ylabel("Kernel density", size=16)
	plt.legend(loc='upper left',frameon=False, fontsize=16)
	plt.tight_layout()
	plt.savefig(prefix+"/riboDensity.png")

	fh = open(prefix+"/list_of_riboLncRNA.txt", "w")
	for tid in lncrnaDensity:
		if lncrnaDensity[tid] >= utr3_cutoff:
			fh.write(tid+"\n")
	fh.close()

parser = argparse.ArgumentParser(description="")  
parser.add_argument('--mrna-lncrna-fa', help='Path to mRNA.lncRNA.fa file', default="/Users/chaozeng/Desktop/ribo-lncRNA/ref/mRNA.lncRNA.fa")
parser.add_argument('--rsem-isoforms-results', help='Path to RSEM isoform results', default="/Users/chaozeng/Desktop/ribo-lncRNA/raw/HeLa_rna.isoforms.results")
parser.add_argument('--ribo-sam', help='Path to Ribo-seq alignment(SAM)', default="/Users/chaozeng/Desktop/ribo-lncRNA/raw/HeLa_ribo.modified.sam")
parser.add_argument('--prefix', help='Path to save the ribosome density plot and the associated raw data', default="/Users/chaozeng/Desktop/ribo-lncRNA")
args = vars(parser.parse_args())

# transCategory,transLen,transORF,utr5Len,utr3Len,cdsLen = getTranscriptBasicInformation(args["mrna_lncrna_fa"])
# transRPKM = getTranscriptExpression(args["rsem_isoforms_results"])
# riboProfile = batchRiboseqAlignmentByRead(args["ribo_sam"],transRPKM,transLen)
# cdsDensity,lncrnaDensity,utr3Density = calculateRiboDensity(transRPKM,transLen,transORF,transCategory,riboProfile,args["prefix"])
cdsDensity,lncrnaDensity,utr3Density = getRiboDensity(args["prefix"])
generatePlotOfRiboDensity(cdsDensity,lncrnaDensity,utr3Density,args["prefix"])
