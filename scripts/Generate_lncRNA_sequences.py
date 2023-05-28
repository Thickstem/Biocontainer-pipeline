#!/bin/python 
import pandas as pd


import argparse
import os.path
import sys
import re
from string import maketrans


lncRNA_biotype = [
	"3prime_overlapping_ncRNA",
	"antisense",
	"lincRNA",
	"non_coding",
	"TEC",
	"non_stop_decay",
	"processed_transcript",
	"retained_intron",
	"sense_intronic",
	"sense_overlapping",
	"bidirectional_promoter_lncRNA",
	"macro_lncRNA",
]


def revcomp(seq):
	if len(seq)>0:
		return seq.translate(maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]
	else: 
		return ""

def parse_gtf(gtf_f):
	tid2info = {}
	with open(gtf_f) as gtf_fh:
		for line in gtf_fh:
			if line[0]=="#": continue
			a = line.split("\t")
			(seqname,feature,start,end,strand,attr_str) = (a[0],a[2],int(a[3]),int(a[4]),a[6],a[8])
			start = int(start)
			end = int(end)
			attr = {}
			for b in attr_str.split(";"):
				m = re.findall('\s*(.*)\s*\"(.*)\"', b)
				if m:
					name = m[0][0].replace(" ","")
					val = m[0][1].replace(" ","")
					if name == "tag":
						if "tag" in attr: 
							attr["tag"].append(val)
						else:
							attr["tag"] = [val]
					else:
						attr[name] = val
			
			if feature=="transcript":
				tid = attr["transcript_id"]
				if "transcript_version" in attr: tid += "."+attr["transcript_version"]
				tname = tid
				if "transcript_name" in attr: tname = attr["transcript_name"].upper()
				tid_tname = tid + "_" + tname
				tid2info[tid] = {"tid":tid_tname}

				gid = attr["gene_id"]
				if "gene_version" in attr: gid += "."+attr["gene_version"]
				gname = gid
				if "gene_name" in attr: gname = attr["gene_name"].upper()
				gid_gname = gid + "_" + gname
				tid2info[tid].update({"gid":gid_gname}) 

				gtype = "" 
				ttype = ""
				if "gene_biotype" in attr: gtype = attr["gene_biotype"]
				if "gene_type" in attr: gtype = attr["gene_type"]
				if "transcript_biotype" in attr: ttype = attr["transcript_biotype"]
				if "transcript_type" in attr: ttype = attr["transcript_type"]
				tcat = "other"
				if ttype=="protein_coding":
					tcat = "mRNA"
					if("cds_start_NF" in attr["tag"] 
						or "cds_end_NF" in attr["tag"]): 
						tcat = "mRNA_incomplete"
				elif(gtype in lncRNA_biotype
					and (ttype=="pseudogene" or ttype in lncRNA_biotype)):
					tcat = ttype
				tid2info[tid].update({"tcat": tcat}) 
				if strand=="+":
					strand=True
				else:
					strand=False
				tid2info[tid].update({"strand": strand}) 
				tid2info[tid].update({"seqname": seqname}) 
				tid2info[tid].update({"exon": []}) 
				tid2info[tid].update({"start_codon": []})
				tid2info[tid].update({"stop_codon": []})
			elif feature=="exon": 
				tid2info[tid]["exon"].append([start, end])
			elif feature=="start_codon": 
				tid2info[tid]["start_codon"].append([start, end])
			elif feature=="stop_codon": 
				tid2info[tid]["stop_codon"].append([start, end])
	return tid2info


def generate_transcript(tid2info,seqname,seq,ribo_noribo_lncrna,prelncrna_seq):
	for tid in ribo_noribo_lncrna:
		tid = "_".join(tid.split("_")[0:-1])
		if tid2info[tid]["seqname"]==seqname:	
			header = tid2info[tid]["tid"]
			header += "\t"+tid2info[tid]["seqname"]
			strand = "-"
			if tid2info[tid]["strand"]: strand = "+"
			header += "\t"+strand
			exon_seq = ""
			for x in tid2info[tid]["exon"]:
				exon_seq += ":".join(map(str,x)) + ","
			header += "\t"+exon_seq
			rna = ""
			for i in range(len(tid2info[tid]["exon"])):
				a = tid2info[tid]["exon"][i]
				if tid2info[tid]["strand"]: 
					if i > 0:
						b = tid2info[tid]["exon"][i-1]
						rna += seq[b[1]:a[0]-1].lower()
					rna += seq[a[0]-1:a[1]].upper()
				else:
					if i > 0:
						b = tid2info[tid]["exon"][i-1]
						rna = seq[a[1]:b[0]-1].lower() + rna
					rna = seq[a[0]-1:a[1]].upper() + rna
			if not tid2info[tid]["strand"]: rna = revcomp(rna)
			prelncrna_seq[tid2info[tid]["tid"]] = header + "\t" + rna + "\n"


def read_genome(tid2info,ribo_noribo_lncrna,genome_f,prelncrna_seq):
	seqname = None
	seq = None
	with open(genome_f) as genome_fh:
		for line in genome_fh:
			if line[0]=="#":
				continue
			if line[0]==">":
				if(seqname is not None): generate_transcript(tid2info,seqname,seq,ribo_noribo_lncrna,prelncrna_seq)
				seqname = line.split()[0][1:]
				seq = ""
			else:
				seq += line.rstrip()
		else:
			generate_transcript(tid2info,seqname,seq,ribo_noribo_lncrna,prelncrna_seq)


## Main
parser = argparse.ArgumentParser(description="Extracting transcript sequences from GTF and genomic sequence files")  
parser.add_argument("--genome", help="genomic sequence file (FASTA)", default="../ref/hg19.fa")
parser.add_argument("--gtf", help="gene annotation file (GTF) of Ensembl or GENCODE", default="../ref/gencode.v25lift37.annotation.gtf")
parser.add_argument("--rna", help="mRNA and lncRNA sequences file (FASTA)", default="../ref/mRNA.lncRNA.fa")
parser.add_argument("--xlsx", help="excel file for ribo-lncRNAs and noribo-lncRNAs", default="../ref/human_ribo_noribo_lncrna.xlsx")
parser.add_argument("--prefix", help="prefix of output files", default="../")
args = vars(parser.parse_args())
if args["prefix"][-1] != "/": args["prefix"] += "/"

## lncRNA class
df = pd.read_excel(args["xlsx"], sheet_name="human")
ribo_noribo_lncrna = []
fh_out = open(args["prefix"]+"lncrna_class.txt","w")
for i in range(0,len(df)):
	if df.iloc[i]["lncRNA_Class"] != "other":
		ribo_noribo_lncrna.append(df.iloc[i]["lncRNA_ID"])
		fh_out.write(df.iloc[i]["lncRNA_ID"] + "\t" + df.iloc[i]["lncRNA_Class"] + "\n")
fh_out.close()

## lncRNA sequences
fh = open(args["rna"])
line_num, tid_flag = 0, ""
lncrna_seq = {}
for line in fh:
	line_num += 1
	if line_num%2 == 1:
		tid = line.split()[0][1:]
		if tid in ribo_noribo_lncrna:
			tid_flag = tid
		else:
			tid_flag = ""
	elif tid_flag != "":
		line = line.rstrip().upper()
		lncrna_seq[tid_flag] = ">" + tid_flag + "\t" + str(len(line)) + "\n" + line + "\n"
fh.close()
fh_out = open(args["prefix"]+"lncrna.fa","w")
for tid in ribo_noribo_lncrna:
	fh_out.write(lncrna_seq[tid])
fh_out.close()

## Premature lncRNA sequences
tid2info = parse_gtf(gtf_f=args["gtf"])	
prelncrna_seq = {}
read_genome(tid2info,ribo_noribo_lncrna,args["genome"],prelncrna_seq)
fh_out = open(args["prefix"]+"prelncrna.txt","w")
fh_out.write("#lncRNA_ID\tchr\tstrand\texon\tseq\n")
for tid in ribo_noribo_lncrna:
	fh_out.write(prelncrna_seq[tid])
fh_out.close()
