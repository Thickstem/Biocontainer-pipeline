#!/bin/python 

lncranLen = {}
fh = open("../ref/human_ribolncrna_noribolncrna.fa")
for line in fh:
	if line[0]==">":
		a = line[1:-1].split()
		lncranLen[a[0]] = int(a[1])
fh.close()
fh_out = open("../ref/human_ribolncrna_noribolncrna.align","w")
fh_out.write( "## compiled from the mapping results of Blast" + "\n")
fh_out.write( "\t".join(["# seq1","seq1_len","seq2","seq2_len","match","match%"]) + "\n")
fh = open("../ref/human_ribolncrna_noribolncrna.blast")
for line in fh:
	a = line.split()
	tid1, tid2, matches = a[0], a[1], int(a[3])
	iden1 = "%.2f"%(100*float(matches)/lncranLen[tid1])
	iden2 = "%.2f"%(100*float(matches)/lncranLen[tid2])
	if tid1 == tid2: continue
	if max([float(iden1), float(iden2)]) < 60: continue
	if lncranLen[tid1] <  lncranLen[tid2]:
		fh_out.write("\t".join([tid2, str(lncranLen[tid2]), tid1, str(lncranLen[tid1]), str(matches), iden1]) + "\n")
	else:
		fh_out.write("\t".join([tid1, str(lncranLen[tid1]), tid2, str(lncranLen[tid2]), str(matches), iden2]) + "\n")
fh.close()
fh_out.close()


# deletedLncrna = set([])
# fh = open("../ref/human_ribolncrna_noribolncrna.align")
# for line in fh:
# 	if line[0] == "#":
# 		continue
# 	a = line.split()
# 	deletedLncrna.add(a[2])
# fh.close()

# fh = open("%s_lncrna_features.txt"%species)
# fhOut = open("%s_lncrna_features.reduced.txt"%species, "w")
# for line in fh:
# 	if line[0] == "#":
# 		fhOut.write(line)
# 	else:
# 		a = line.split()
# 		if a[1] == "trans-lncRNA":
# 			deletedLncrna.add(a[0])
# 		if a[0] not in deletedLncrna:
# 			fhOut.write(line)
# fh.close()
# fhOut.close()

# fh = open("%s_lncrna_orfs.txt"%species)
# fhOut = open("%s_lncrna_orfs.reduced.txt"%species, "w")
# for line in fh:
# 	if line[0] == "#":
# 		fhOut.write(line)
# 	elif line.split()[0] in deletedLncrna:
# 		None
# 	else:
# 		fhOut.write(line)
# fh.close()
# fhOut.close()

