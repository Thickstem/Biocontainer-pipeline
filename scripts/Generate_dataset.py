import math
import numpy as np
from statistics import mean, median
import random

def load_lncrna_length():
	lncrnaLen = {}
	fh = open("../lncrna.fa")
	for line in fh:
		if line[0] == ">": 
			a = line[1:-1].split()
			lncrnaLen[a[0]] = float(a[1])
	fh.close()
	return lncrnaLen

def load_sequence():
	lncrnaSeq = {}
	fh = open("../lncrna_orfs.reduced.txt")
	fh.next()
	for line in fh:
		a = line.rstrip().split()
		lncrnaSeq[a[0]] = {"pOrfStartContext":a[1],"fOrfStartContext":a[2],"uOrfStartContext":a[3],
			"pOrfSeq":a[4],"fOrfSeq":a[5],"uOrfSeq":a[6]}
	fh.close()
	return lncrnaSeq

def load_lncrna_features():
	lncrnaLen = load_lncrna_length()
	fh = open("../lncrna_features.reduced.txt")
	colName = fh.readline()[2:-1].split("\t")
	df = {}
	lncrnaNum = 0
	TElist = set([])
	for colname in colName:
		df[colname] = []
	for line in fh:
		lncrnaNum += 1
		a = line[:-1].split("\t")
		for i,x in enumerate(a):
			try:
				df[colName[i]].append(float(x))
			except: 
				df[colName[i]].append(x)
		if df["TE"][-1] != "-":
			TElist.update(df["TE"][-1].split(";"))
	fh.close()

	for x in sorted(TElist):
		colName.append(x)
		df[x] = [0] * lncrnaNum

	extendList = ["m6aNearTIS_log","m6aNearTTS_log","m6aNearpORFstart_log","m6aNearfORFstart_log",
		"m6aNearuORFstart_log","m6aNearpORFend_log","m6aNearfORFend_log","m6aNearuORFend_log",
		"g4NearTIS_log","g4NearTTS_log","g4NearpORFstart_log","g4NearfORFstart_log",
		"g4NearuORFstart_log","g4NearpORFend_log","g4NearfORFend_log","g4NearuORFend_log",
		"m6aNearTIS_%","m6aNearTTS_%","m6aNearpORFstart_%","m6aNearfORFstart_%",
		"m6aNearuORFstart_%","m6aNearpORFend_%","m6aNearfORFend_%","m6aNearuORFend_%",
		"g4NearTIS_%","g4NearTTS_%","g4NearpORFstart_%","g4NearfORFstart_%","g4NearuORFstart_%",
		"g4NearpORFend_%","g4NearfORFend_%","g4NearuORFend_%"]
	for x in extendList:
		colName.append(x)
		df[x] = ["-"] * lncrnaNum

	# cleaning data
	for i in range(lncrnaNum):
		if df["class"][i] == "ribo-lncRNA": 
			df["class"][i] = 1
		else:
			df["class"][i] = 0

		if df["pOrfLen"][i] == "-": df["pOrfLen"][i] = 0
		if df["fOrfLen"][i] == "-": df["fOrfLen"][i] = 0
		if df["uOrfLen"][i] == "-": df["uOrfLen"][i] = 0
		df["pOrfLen"][i] = round(math.log((1+df["pOrfLen"][i]),10),3)
		df["fOrfLen"][i] = round(math.log((1+df["fOrfLen"][i]),10),3)
		df["uOrfLen"][i] = round(math.log((1+df["uOrfLen"][i]),10),3)

		if df["TE"][i] != "-":
			for x in df["TE"][i].split(";"):
				df[x][i] = 1

		tid = df["lncrna"][i]
		if df["m6A"][i] != "-":
			b = np.array(map(float,str(df["m6A"][i]).split(";")))
			df["m6aNearTIS_log"][i] = round(math.log((min(b)+1),10),3)
			df["m6aNearTTS_log"][i] = round(math.log((abs(lncrnaLen[tid]-max(b))+1),10),3)
			df["m6aNearTIS_%"][i] = round(min(b)/float(lncrnaLen[tid]),3)
			df["m6aNearTTS_%"][i] = round(abs(lncrnaLen[tid]-max(b))/float(lncrnaLen[tid]),3)
			if df["pOrfStart"][i] != "-":
				df["m6aNearpORFstart_log"][i] = round(math.log(min(abs(b-int(df["pOrfStart"][i]))+1),10),3)
				df["m6aNearpORFstart_%"][i] = round(min(abs(b-int(df["pOrfStart"][i])))/float(lncrnaLen[tid]),3)
			if df["pOrfEnd"][i] != "-":
				df["m6aNearpORFend_log"][i] = round(math.log(min(abs(b-int(df["pOrfEnd"][i]))+1),10),3)
				df["m6aNearpORFend_%"][i] = round(min(abs(b-int(df["pOrfEnd"][i])))/float(lncrnaLen[tid]),3)
			if df["fOrfStart"][i] != "-":
				df["m6aNearfORFstart_log"][i] = round(math.log(min(abs(b-int(df["fOrfStart"][i]))+1),10),3)
				df["m6aNearfORFstart_%"][i] = round(min(abs(b-int(df["fOrfStart"][i])))/float(lncrnaLen[tid]),3)
			if df["fOrfEnd"][i] != "-":
				df["m6aNearfORFend_log"][i] = round(math.log(min(abs(b-int(df["fOrfEnd"][i]))+1),10),3)
				df["m6aNearfORFend_%"][i] = round(min(abs(b-int(df["fOrfEnd"][i])))/float(lncrnaLen[tid]),3)
			if df["uOrfStart"][i] != "-":
				df["m6aNearuORFstart_log"][i] = round(math.log(min(abs(b-int(df["uOrfStart"][i]))+1),10),3)
				df["m6aNearuORFstart_%"][i] = round(min(abs(b-int(df["uOrfStart"][i])))/float(lncrnaLen[tid]),3)
			if df["uOrfEnd"][i] != "-":
				df["m6aNearuORFend_log"][i] = round(math.log(min(abs(b-int(df["uOrfEnd"][i]))+1),10),3)
				df["m6aNearuORFend_%"][i] = round(min(abs(b-int(df["uOrfEnd"][i])))/float(lncrnaLen[tid]),3)
		if df["G4"][i] != "-":
			b = np.array(map(float,str(df["G4"][i]).split(";")))
			df["g4NearTIS_log"][i] = round(math.log((min(b)+1),10),3)
			df["g4NearTTS_log"][i] = round(math.log((abs(lncrnaLen[tid]-max(b))+1),10),3)
			df["g4NearTIS_%"][i] = round(min(b)/float(lncrnaLen[tid]),3)
			df["g4NearTTS_%"][i] = round(abs(lncrnaLen[tid]-max(b))/float(lncrnaLen[tid]),3)
			if df["pOrfStart"][i] != "-":
				df["g4NearpORFstart_log"][i] = round(math.log(min(abs(b-int(df["pOrfStart"][i]))+1),10),3)
				df["g4NearpORFstart_%"][i] = round(min(abs(b-int(df["pOrfStart"][i])))/float(lncrnaLen[tid]),3)
			if df["pOrfEnd"][i] != "-":
				df["g4NearpORFend_log"][i] = round(math.log(min(abs(b-int(df["pOrfEnd"][i]))+1),10),3)
				df["g4NearpORFend_%"][i] = round(min(abs(b-int(df["pOrfEnd"][i])))/float(lncrnaLen[tid]),3)
			if df["fOrfStart"][i] != "-":
				df["g4NearfORFstart_log"][i] = round(math.log(min(abs(b-int(df["fOrfStart"][i]))+1),10),3)
				df["g4NearfORFstart_%"][i] = round(min(abs(b-int(df["fOrfStart"][i])))/float(lncrnaLen[tid]),3)
			if df["fOrfEnd"][i] != "-":
				df["g4NearfORFend_log"][i] = round(math.log(min(abs(b-int(df["fOrfEnd"][i]))+1),10),3)
				df["g4NearfORFend_%"][i] = round(min(abs(b-int(df["fOrfEnd"][i])))/float(lncrnaLen[tid]),3)
			if df["uOrfStart"][i] != "-":
				df["g4NearuORFstart_log"][i] = round(math.log(min(abs(b-int(df["uOrfStart"][i]))+1),10),3)
				df["g4NearuORFstart_%"][i] = round(min(abs(b-int(df["uOrfStart"][i])))/float(lncrnaLen[tid]),3)
			if df["uOrfEnd"][i] != "-":
				df["g4NearuORFend_log"][i] = round(math.log(min(abs(b-int(df["uOrfEnd"][i]))+1),10),3)
				df["g4NearuORFend_%"][i] = round(min(abs(b-int(df["uOrfEnd"][i])))/float(lncrnaLen[tid]),3)

	colName.remove("chr")
	del df["chr"]
	colName.remove("strand")
	del df["strand"]
	colName.remove("TE")
	del df["TE"]
	colName.remove("m6A")
	del df["m6A"]
	colName.remove("G4")
	del df["G4"]
	colName.remove("pOrfStart")
	del df["pOrfStart"]
	colName.remove("pOrfEnd")
	del df["pOrfEnd"]
	colName.remove("fOrfStart")
	del df["fOrfStart"]
	colName.remove("fOrfEnd")
	del df["fOrfEnd"]
	colName.remove("uOrfStart")
	del df["uOrfStart"]
	colName.remove("uOrfEnd")
	del df["uOrfEnd"]

	# imputation
	for x in colName:
		if x == "lncrna" or x == "class": continue
		tmp = []
		for i in range(lncrnaNum):
			if df[x][i] != "-":
				tmp.append(df[x][i])
		
		m = round(mean(tmp),3)
		for i in range(lncrnaNum):
			if df[x][i] == "-":
				df[x][i] = m

	# save imputated data
	# fh = open("train/%s_imputated.txt"%, "w")
	# fh.write("\t".join(colName) + "\n")
	# for i in range(lncrnaNum):
	# 	fh.write("\t".join([str(df[x][i]) for x in colName]) + "\n")
	# fh.close()

	return lncrnaNum, colName, df


def add_context_trimer_hexamer(data,seq,):
	h = {"context":{}, "trimer":{}, "hexamer":{}}
	fh = open("../5000_ratio.txt")
	fh.next()
	for line in fh:
		a = line.rstrip().split()
		if a[1] != "-":
			a[1] = int(a[1]) + 6
		if a[1] not in h[a[0]]: h[a[0]][a[1]] = {}
		h[a[0]][a[1]][a[2]] = float(a[-1])
	fh.close()

	data2 = []
	for x in data:
		score = []
		for y in ["pOrfStartContext", "fOrfStartContext", "uOrfStartContext"]:
			if seq[x[0]][y] == "-":
				score.append(0)
			else:
				n,s = 0,0
				for i in range(10):
					n += 1
					nt = seq[x[0]][y][i]
					if nt in h["context"][i]: s+=h["context"][i][nt]
				score.append(round(s/n,3))
		for y in ["pOrfSeq","fOrfSeq","uOrfSeq"]:
			L = len(seq[x[0]][y])
			if seq[x[0]][y] == "-" or L <= 6:
				score.append(0)
				score.append(0)
			else:
				s1,n1,s2,n2 = 0,0,0,0
				for i in range(0,L,3):
					n1 += 1
					trimer = seq[x[0]][y][i:i+3]
					if trimer in h["trimer"]["-"]: s1 += h["trimer"]["-"][trimer]
					if i+6 < L:
						hexamer = seq[x[0]][y][i:i+6]
						n2 += 1
						if hexamer in h["hexamer"]["-"]: s2 += h["hexamer"]["-"][hexamer]
				score.append(round(s1/n1,3))
				score.append(round(s2/n2,3))
		data2.append(score)
	col2 = ["pOrfStartContext", "fOrfStartContext", "uOrfStartContext", "pOrfSeqTrimer", "pOrfSeqHexamer", 
		"fOrfSeqTrimer", "fOrfSeqHexamer", "uOrfSeqTrimer", "uOrfSeqHexamer"]
	return data2, col2


def save_data():
	N, col1, df = load_lncrna_features()
	seq = load_sequence()
	filePath = "../train_data.txt"

	data1 = []
	for i in range(N):
		tmp = [df[x][i] for x in col1]
		data1.append(tmp)

	data2, col2 = add_context_trimer_hexamer(data1,seq,)

	fh = open(filePath, "w")
	fh.write("\t".join(col1[1:]) + "\t" + "\t".join(col2) + "\n")
	for i in range(len(data1)):
		fh.write("\t".join(map(str,data1[i][1:])) + "\t" + "\t".join(map(str,data2[i])) + "\n")
	fh.close()	


save_data()


