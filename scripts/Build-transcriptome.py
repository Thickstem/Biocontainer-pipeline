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


def generate_transcript(tid2info,seqname,seq,only_codRNA):
	for tid in tid2info:
		if tid2info[tid]["seqname"]==seqname:
			if tid2info[tid]["tcat"] == "mRNA":
				header = ">"+tid2info[tid]["tid"]
				header += "\t"+"gene_id="+tid2info[tid]["gid"]
			
				start_codon = -1
				stop_codon = -1
				for a in tid2info[tid]["start_codon"]:
					if tid2info[tid]["strand"]:
						if start_codon==-1 or start_codon<a[0]:
							start_codon = a[0]
					else:
						if start_codon==-1 or start_codon>a[1]:
							start_codon = a[1]
				for a in tid2info[tid]["stop_codon"]:
					if tid2info[tid]["strand"]:
						if stop_codon==-1 or stop_codon<a[0]:
							stop_codon = a[0]
					else:
						if stop_codon==-1 or stop_codon>a[1]:
							stop_codon = a[1]
				
				rna = ""
				cds_start = 0
				cds_stop = 0
				for a in tid2info[tid]["exon"]:
					if tid2info[tid]["strand"]: 
						rna += seq[a[0]-1:a[1]]
						
						if start_codon>a[1]: 
							cds_start += a[1]-a[0]+1
						elif(start_codon>=a[0]
							and start_codon<=a[1]):
							cds_start += start_codon - a[0] 
					
						if stop_codon>a[1]: 
							cds_stop += a[1]-a[0]+1
						elif(stop_codon>=a[0]
							and stop_codon<=a[1]):
							cds_stop += stop_codon - a[0] 
					else:
						rna = seq[a[0]-1:a[1]] + rna
						
						if start_codon<a[0]: 
							cds_start += a[1]-a[0]+1
						elif(start_codon>=a[0]
							and start_codon<=a[1]):
							cds_start += a[1] - start_codon 
					
						if stop_codon<a[0]: 
							cds_stop += a[1]-a[0]+1
						elif(stop_codon>=a[0]
							and stop_codon<=a[1]):
							cds_stop += a[1] - stop_codon 
				if not tid2info[tid]["strand"]: rna = revcomp(rna)	
				if start_codon>=0 and stop_codon>=0:
					header += "\t"+"tcat="+tid2info[tid]["tcat"]
					header += "\t"+"cds=" + str(cds_start) + "-" + str(cds_stop)
					rna = rna[0:cds_start].lower() + rna[cds_start:cds_stop+3].upper() + rna[cds_stop+3:].lower()
				elif start_codon>=0:
					if tid2info[tid]["tcat"]=="mRNA": 
						header += "\t"+"tcat="+tid2info[tid]["tcat"]+"_incomplete"
						tid2info[tid]["tcat"]= tid2info[tid]["tcat"]+"_incomplete"
					else: header += "\t"+"tcat="+tid2info[tid]["tcat"]
					header += "\t"+"cds=" + str(cds_start) + "-?" 
					rna = rna[0:cds_start].lower() + rna[cds_start:].upper()
				elif stop_codon>=0:
					if tid2info[tid]["tcat"]=="mRNA": 
						header += "\t"+"tcat="+tid2info[tid]["tcat"]+"_incomplete"
						tid2info[tid]["tcat"]= tid2info[tid]["tcat"]+"_incomplete"
					else: header += "\t"+"tcat="+tid2info[tid]["tcat"]
					header += "\t"+"cds=" + "?-" + str(cds_stop)
					rna = rna[0:cds_stop+3].upper() + rna[cds_stop+3:].lower()
				else:
					if tid2info[tid]["tcat"]=="mRNA": 
						header += "\t"+"tcat="+tid2info[tid]["tcat"]+"_incomplete"
						tid2info[tid]["tcat"]= tid2info[tid]["tcat"]+"_incomplete"
					else: header += "\t"+"tcat="+tid2info[tid]["tcat"]
					rna = rna.lower()

				if tid2info[tid]["tcat"] == "mRNA":
					print header
					print rna

			elif tid2info[tid]["tcat"] != "other" and tid2info[tid]["tcat"] != "mRNA_incomplete":
				if not only_codRNA:
					header = ">"+tid2info[tid]["tid"]
					header += "\t"+"gene_id="+tid2info[tid]["gid"]
					header += "\t"+"tcat="+tid2info[tid]["tcat"]
					rna = ""
					for a in tid2info[tid]["exon"]:
						if tid2info[tid]["strand"]: 
							rna += seq[a[0]-1:a[1]]
						else:
							rna = seq[a[0]-1:a[1]] + rna
					if not tid2info[tid]["strand"]: rna = revcomp(rna)
					print header
					print rna.lower()


def read_genome(tid2info,genome_f,only_codRNA=False):
	seqname = None
	seq = None
	with open(genome_f) as genome_fh:
		for line in genome_fh:
			if line[0]=="#":
				continue
			if line[0]==">":
				if(seqname is not None): generate_transcript(tid2info,seqname,seq,only_codRNA)
				seqname = line.split()[0][1:]
				seq = ""
			else:
				seq += line.rstrip()
		else:
			generate_transcript(tid2info,seqname,seq,only_codRNA)

## Main
parser = argparse.ArgumentParser(description="Extracting transcript sequences from GTF and genomic sequence files")  
parser.add_argument("--genome", help="genomic sequence file (FASTA)", required=True)
parser.add_argument("--gtf", help="gene annotation file (GTF) of Ensembl or GENCODE", required=True)
args = vars(parser.parse_args())
tid2info = parse_gtf(gtf_f=args["gtf"])
read_genome(tid2info,genome_f=args["genome"],only_codRNA=False)





















