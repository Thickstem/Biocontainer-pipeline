#!/bin/python

deletedLncrna = set([])
fh = open("../lncrna_align.txt")
for line in fh:
	if line[0] == "#":
		continue
	a = line.split()
	deletedLncrna.add(a[2])
fh.close()

fh = open("../lncrna_features.txt")
fhOut = open("../lncrna_features.reduced.txt", "w")
for line in fh:
	if line[0] == "#":
		fhOut.write(line)
	else:
		a = line.split()
		if a[1] == "trans-lncRNA":
			deletedLncrna.add(a[0])
		if a[0] not in deletedLncrna:
			fhOut.write(line)
fh.close()
fhOut.close()

fh = open("../lncrna_orfs.txt")
fhOut = open("../lncrna_orfs.reduced.txt", "w")
for line in fh:
	if line[0] == "#":
		fhOut.write(line)
	elif line.split()[0] in deletedLncrna:
		None
	else:
		fhOut.write(line)
fh.close()
fhOut.close()

