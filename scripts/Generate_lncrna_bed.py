#!/bin/python
import os

fhOut = open("../lncrna.bed", "w") 
fh = open("../prelncrna.txt")
fhOut.write("#chromsome\texonStart\texonEnd\trna\tclass\tstrand\n")
for line in fh:
	a = line[:-1].split()
	[tid, chro, strand, exon] = a[:4]
	lncrnaType = "2"
	exonList = exon.split(",")[:-1]
	for x in exonList:
		b = x.split(":")
		tmp = [chro, b[0], b[1], tid, lncrnaType, strand]
		fhOut.write("\t".join(tmp)+"\n")
fhOut.close()
fh.close()
os.system("sort -k1,1 -k2,2n ../lncrna.bed > ../lncrna.sorted.bed && rm ../lncrna.bed")


fhOut = open("../repeatmask.bed.tmp", "w")
fh = open("../repeatmask.out")
fh.next()
fh.next()
fh.next()
for line in fh:
	a = line.split()
	[chrom,strat,end,strand,repeatElement,repeatType] = [a[4],a[5],a[6],a[8],a[9],a[10]]
	if strand == "C": strand = "-"
	if "/" in repeatType:
		[repeatClass, repeatFamily] = repeatType.split("/")
		tmp = [chrom,strat,end,repeatElement,strand,repeatFamily,repeatClass]
		fhOut.write("\t".join(tmp)+"\n")
fh.close()
fhOut.close()
os.system("echo '# chromosome\tstart\tend\tlncrna\tlncrna_class\tstrand\tTE\tTE_family\tTE_class' >> ../lncrna_repeatmask.bed")
os.system("sort -k1,1 -k2,2n ../repeatmask.bed.tmp > ../repeatmask.bed.sorted.tmp")
fLncrna = "../lncrna.sorted.bed"
cmd = "bedtools intersect -wb -a ../lncrna.sorted.bed -b ../repeatmask.bed.sorted.tmp | awk \'{if($6==$11){print $0}}\' | cut -f1-6,10,12,13 >> ../lncrna_repeatmask.bed"
os.system(cmd)
os.system("rm ../repeatmask.bed.tmp ../repeatmask.bed.sorted.tmp")