#!/bin/python
import os

os.system("mkdir -p ../tmp")
lncrnaFastaFile = "../lncrna.fa"
lncrnaFastaFH = open(lncrnaFastaFile)
a, n, fileList, lenList, tidList = None, 1, [], [], []
for line in lncrnaFastaFH:
	if len(line) < 1: continue
	if line[0] == ">":
		a = line[1:-1].split()
		tidList.append(a[0])
	else:
		f = "../tmp/%d.fa"%(n)
		fh = open(f,"w")
		fh.write(">" + a[0] + "\t" + a[1] + "\n")
		fh.write(line)
		fh.close()
		n += 1
		fileList.append(f)
		lenList.append(int(a[1])-1)
lncrnaFastaFH.close()

# jobFH = open("../parasoR.sh", "w")
# for i in range(len(fileList)):
# 	cmd = "ParasoR --pre --stem --input %s > %s.sp"%(fileList[i],fileList[i])
# 	jobFH.write(cmd + "\n")
# jobFH.close()
# os.system("bash ../parasoR.sh")

# fh = open("../lncrna_stemprob_parasor.txt", "w")
# for i,tid in enumerate(tidList):
# 	tmp = [tid]
# 	f1 = fileList[i]
# 	f2 = fileList[i] + ".sp"
# 	fh2 = open(f2)
# 	for line in fh2:
# 		if len(line)<1 or line[0]=="#": continue 
# 		prob = float(line.split()[-1])
# 		tmp.append("%.3f"%prob)
# 	fh2.close()
# 	fh.write("\t".join(tmp) + "\n")
# fh.close()

jobFH = open("../qgrs.sh", "w")
for i in range(len(fileList)):
	cmd = "qgrs -i %s > %s.g4"%(fileList[i],fileList[i])
	jobFH.write(cmd + "\n")
jobFH.close()
os.system("bash ../qgrs.sh")

fh = open("../lncrna_g4_qgrs.txt", "w")
for i,tid in enumerate(tidList):
	tmp = [tid]
	f1 = fileList[i]
	f2 = fileList[i] + ".g4"
	fh2 = open(f2)
	fh2.next()
	fh2.next()
	fh2.next()
	for line in fh2:
		a = line[:-1].split()
		tmp.append("\t".join(a[1:]))
		fh.write("\t".join(tmp) + "\n")
	fh2.close()
fh.close()