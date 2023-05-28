import numpy as np
import math

def main():
	fh = open("../lncrna.fa")
	lncrnaClass = load_lncrna_class()
	prelncrnaFeatures = extract_prelncrna_features()
	stemProb = load_stem_prob_file()
	m6Asramp,g4qgrs,repeatRmsk = load_m6A_g4_rmsk_files(lncrnaClass)
	
	fhOut1 = open("../lncrna_features.txt", "w")
	fhOut2 = open("../lncrna_orfs.txt", "w")
	head1 = ["# lncrna", "class", "chr", "strand", "fLen", "gc", "nE", 
		"fELen", "minELen", "maxELen", "avgELen", "fEgc", "minEgc", "maxEgc", "avgEgc", 
		"fILen", "minILen", "maxILen", "avgILen", "fIgc", "minIgc", "maxIgc", "avgIgc",
		"pOrfStart","pOrfEnd","pOrfLen","pOrfCov","pOrfSp",
		"pOrf5utrLen","pOrf5utrCov","pOrf5utrSp","pOrf5utrSpFC",
		"pOrf3utrLen","pOrf3utrCov","pOrf3utrSp","pOrf3utrSpFC",
		"fOrfStart","fOrfEnd","fOrfLen","fOrfCov","fOrfSp",
		"fOrf5utrLen","fOrf5utrCov","fOrf5utrSp","fOrf5utrSpFC",
		"fOrf3utrLen","fOrf3utrCov","fOrf3utrSp","fOrf3utrSpFC",
		"uOrfStart","uOrfEnd","uOrfLen","uOrfCov","uOrfSp",
		"uOrf5utrLen","uOrf5utrCov","uOrf5utrSp","uOrf5utrSpFC",
		"uOrf3utrLen","uOrf3utrCov","uOrf3utrSp","uOrf3utrSpFC",
		"m6A","G4","TE"]
	fhOut1.write("\t".join(head1) + "\n")
	head2 = ["# lncrna", "pOrfStartContext", "fOrfStartContext", "uOrfStartContext", "pOrfSeq", "fOrfSeq", "uOrfSeq"]
	fhOut2.write("\t".join(head2) + "\n") 
	for line in fh:
		if line[0]==">":
			tid = line[1:].split()[0]
		else:
			pORF, fORF, uORF = None, None, None
			pORF = find_orf(line[:-1])
			if pORF is not None:
				fORF = find_orf(line[:-1], bound=pORF[1]-1, longest=False)
				uORF = find_orf(line[:-1], bound=pORF[0], startCodon=["CTG","GTG","TTG"])
			pORFfeatures, pORFseq, pORFStartContext = extract_orf_features(pORF, line[:-1], stemProb[tid])
			fORFfeatures, fORFseq, fORFStartContext = extract_orf_features(fORF, line[:-1], stemProb[tid])
			uORFfeatures, uORFseq, uORFStartContext = extract_orf_features(uORF, line[:-1], stemProb[tid])	
			otherfeatures = "\t".join([get_m6A_g4_rmsk_str(m6Asramp, tid), get_m6A_g4_rmsk_str(g4qgrs, tid),
				get_m6A_g4_rmsk_str(repeatRmsk, tid)])
			fhOut1.write(tid + "\t" + lncrnaClass[tid] + "\t" + prelncrnaFeatures[tid] + "\t" + 
				pORFfeatures + "\t" + fORFfeatures + "\t" + uORFfeatures + "\t" + otherfeatures + "\n")
			fhOut2.write(tid + "\t" + pORFStartContext  + "\t" + fORFStartContext + "\t" + uORFStartContext + "\t" + 
				pORFseq + "\t" + fORFseq + "\t" + uORFseq + "\n")
			
	fh.close()
	fhOut1.close()
	fhOut2.close()

def load_lncrna_class():
	lncrnaClass = {}
	fh = open("../lncrna_class.txt")
	for line in fh:
		a = line[:-1].split()
		lncrnaClass[a[0]] = a[1]
	fh.close()
	return lncrnaClass

def find_orf(seq, bound=0, longest=True, startCodon=["ATG"], stopCodon=['TAG','TAA','TGA']):
	seqLen = len(seq)
	if bound == 0:
		bound = seqLen
	seq = seq.upper()
	orfStart, orfEnd, orfLen = 0,0,0
	for i in range(bound - 6):
		c1 = seq[i:i+3]
		if c1 in startCodon:
			for j in range(i+3, bound-3, 3):
				c2 = seq[j:j+3]
				if c2 in stopCodon:
					if longest:
						if (j-i) > orfLen:
							orfStart = i
							orfEnd = j
							orfLen = j - i
							break
					else:
						orfStart = i
						orfEnd = j
						orfLen = j - i
						break
		if not longest and orfLen>0:
			break
	if orfLen > 0:
		orfEnd += 2
		orfLen += 3
		return [orfStart, orfEnd, orfLen]
	else:
		return None

def load_stem_prob_file():
	stemProb = {}
	f = open("../lncrna_stemprob_parasor.txt")
	for line in f:
		if line[0]=="#": continue
		a = line[:-1].split()
		stemProb[a[0]] = map(float,a[1:])
	f.close()
	return stemProb

def extract_prelncrna_features():
	prelncrna = {}
	fh = open("../prelncrna.txt")
	fh.next()
	for line in fh:
		[tid,chrosome,strand,exonStr,seq] = line[:-1].split()
		e = exonStr.split(",")[:-1]
		if strand == "-":
			e.reverse()
		nE = len(e)
		eList = []
		for a in e:
			b = map(int,a.split(":"))
			eList.append(b[0])
			eList.append(b[1])
		eLenList = []
		for i in range(1,2*nE,2):
			eLenList.append(eList[i]-eList[i-1]+1)
		iLenList = []
		for i in range(2,2*nE,2):
			iLenList.append(eList[i]-eList[i-1]-1)
		if strand == "-":
			eLenList.reverse()
			iLenList.reverse()
		eGcList, iGcList = [], []
		k = 0
		rna = ""
		for i, elen in enumerate(eLenList):
			exon = seq[k:k+elen]
			rna += exon
			eGcList.append(float(exon.count("G")+exon.count("C"))/elen)
			k += elen
			try:
				ilen = iLenList[i]
				intron = seq[k:k+ilen]
				iGcList.append(float(intron.count("g")+intron.count("c"))/ilen)
				k += ilen
			except:
				None
		fLen = len(rna)
		gc = "%.3f"%(float(rna.count("G")+rna.count("C"))/fLen)
		fILen, minILen, maxILen, avgILen, fIgc, minIgc, maxIgc, avgIgc = 0, 0, 0, 0, 0, 0, 0, 0
		if len(iLenList) > 0:
			fILen = "%.3f"%math.log(iLenList[0]+1,10)
			minILen = "%.3f"%math.log(min(iLenList)+1,10)
			maxILen = "%.3f"%math.log(max(iLenList)+1,10)
			avgILen = "%.3f"%math.log(np.mean(iLenList)+1,10)
			fIgc = "%.3f"%iGcList[0]
			minIgc = "%.3f"%(min(iGcList))
			maxIgc = "%.3f"%(max(iGcList))
			avgIgc = "%.3f"%(np.mean(iGcList))
		fELen = "%.3f"%math.log(eLenList[0]+1,10)
		minELen = "%.3f"%math.log(min(eLenList)+1,10)
		maxELen = "%.3f"%math.log(max(eLenList)+1,10)
		avgELen = "%.3f"%math.log(np.mean(eLenList)+1,10)
		fEgc = "%.3f"%eGcList[0]
		minEgc = "%.3f"%(min(eGcList))
		maxEgc = "%.3f"%(max(eGcList))
		avgEgc = "%.3f"%(np.mean(eGcList))
		fLen = "%.3f"%math.log(fLen, 10)
		prelncrna[tid] = "\t".join(map(str, [chrosome, strand, fLen, gc, nE, fELen, minELen, maxELen, avgELen, fEgc, minEgc, maxEgc, avgEgc, 
			fILen, minILen, maxILen, avgILen, fIgc, minIgc, maxIgc, avgIgc]))
	fh.close()
	return prelncrna

def extract_orf_features(orfInfo, seq, stemProb):
	if orfInfo is not None:
		[orfStart, orfEnd, orfLen] = orfInfo
		rnaLen = float(len(seq))
		orfCov = "%.3f"%(orfLen/rnaLen)
		orfSp = "-" 
		if len(stemProb) > 0: orfSp = np.mean(stemProb[orfStart:orfEnd+1]) 
		utr5Len = orfStart
		utr5Cov = "%.3f"%(utr5Len/rnaLen)
		utr5Sp, utr5SpFC = "-", "-"
		if utr5Len >= 30 and len(stemProb) > 0: 
			utr5Sp = np.mean(stemProb[0:orfStart+1]) 
			utr5SpFC = math.log(utr5Sp/orfSp,2)
			utr5Sp = "%.3f"%(utr5Sp)
			utr5SpFC = "%.3f"%(utr5SpFC)
		utr3Len = int(rnaLen - orfEnd - 1) 
		utr3Cov = "%.3f"%(utr3Len/rnaLen)

		utr3Sp, utr3SpFC = "-", "-"
		if utr3Len >= 30 and len(stemProb) > 0: 
			utr3Sp = np.mean(stemProb[orfEnd+1:]) 
			utr3SpFC = math.log(utr3Sp/orfSp,2)
			utr3Sp = "%.3f"%(utr3Sp)
			utr3SpFC = "%.3f"%(utr3SpFC)
		orfStartContext = "-"
		if utr5Len >= 6:
			orfStartContext = seq[orfStart-6:orfStart+4]
		orfSeq = seq[orfStart:orfEnd+1]
		if orfSp != "-": orfSp = "%.3f"%(orfSp)
		utr5Len = "%.3f"%math.log(utr5Len + 1)
		utr3Len = "%.3f"%math.log(utr3Len + 1)
		return "\t".join(map(str, [orfStart, orfEnd, orfLen, orfCov, orfSp,
			utr5Len, utr5Cov, utr5Sp, utr5SpFC, utr3Len, utr3Cov, utr3Sp, utr3SpFC])), orfSeq, orfStartContext
	else:
		return "\t".join(["-"]*13), "-", "-"

def load_m6A_g4_rmsk_files(lncrnaClass):
	tid2tidFull = {}
	for tid in lncrnaClass:
		tid2tidFull[tid.split(".")[0]] = tid
	m6Asramp = {}
	fHm6Asramp = open("../lncrna_m6a_sramp.txt")
	fHm6Asramp.readline()
	for line in fHm6Asramp:
		if line[0]=="#": continue
		if "igh confidence" not in line: continue # High or very high confidence site will be selected
		a = line.split()
		tid = tid2tidFull[a[0].split(".")[0]]
		if tid not in m6Asramp:
			m6Asramp[tid] = []
		m6Asramp[tid].append(a[1])
	fHm6Asramp.close()

	g4qgrs = {}
	fHg4 = open("../lncrna_g4_qgrs.txt")
	for line in fHg4:
		if line[0]=="#": continue
		a = line.split()
		tid = a[0]
		if int(a[6]) < 30: continue
		if tid not in g4qgrs:
			g4qgrs[tid] = []
		g4qgrs[tid].append(a[1])
	fHg4.close()

	repeatRmsk = {}
	fHrmsk = open("../lncrna_repeatmask.bed")
	for line in fHrmsk:
		if line[0]=="#": continue
		a = line.rstrip().split()
		tid = a[3]
		teClass = a[-1]
		if "?" in teClass:
			teClass = teClass[:-1]
		if tid not in repeatRmsk:
			repeatRmsk[tid] = set([])
		repeatRmsk[tid].add(teClass)
	fHrmsk.close()

	return m6Asramp,g4qgrs,repeatRmsk

def get_m6A_g4_rmsk_str(data, tid):
	if tid not in data:
		return "-"
	return ";".join(data[tid])

main()
