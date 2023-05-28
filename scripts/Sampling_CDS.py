import numpy as np
import random
import math

fh = open("../ref/mrna.lncrna.fa")
data = []
for line in fh:
	if line[0]==">":
		a = line[1:-1].split()
	else:
		if len(a)==4:
			c = map(int,a[-1][4:].split("-"))
			if c[0]-6>0:
				tid = a[0]
				context = line[:-1][c[0]-6:c[0]+4].upper()
				cds = line[:-1][c[0]:c[1]+3]
				if(cds[:3]=="ATG" and cds[-3:] in ["TAG","TGA","TAA"] and (len(cds)-6)%6==0):
					flag = 1
					for i in range(0,len(cds)-3,3):
						if cds[i:i+3] in ["TAG","TGA","TAA"]:
							flag = 0
							break
					if flag == 1:
						data.append([tid, context, cds[3:-3]])
fh.close()

datIndex = np.random.choice(len(data), 5000, replace=False)

fh = open("../5000_cds.txt", "w")
fh.write("\t".join(["# lncrna","context sequecne","cds"]) + "\n")
for i in datIndex:
	fh.write("\t".join(data[i]) + "\n")
fh.close()

context = [{},{}]
trimer = [{},{}]
hexamer = [{},{}]
fh = open("../5000_random.txt", "w")
fh.write("\t".join(["# lncrna","context sequecne","cds"]) + "\n")
for i in datIndex:
	a = data[i][1][:-4] + data[i][2]
	b = list(a)
	random.shuffle(b)
	a2 = "".join(b[:6]) + "ATG" + b[6]
	b2 = "".join(b[6:])
	fh.write("\t".join([data[i][0], a2, b2]) + "\n")

	for j in range(10):
		x1 = data[i][1][j]
		x2 = a2[j]
		if j not in context[0]: context[0][j] = {}
		if j not in context[1]: context[1][j] = {}
		if x1 not in context[0][j]: context[0][j][x1] = 0 
		if x2 not in context[1][j]: context[1][j][x2] = 0 
		context[0][j][x1] += 1
		context[1][j][x2] += 1

	L = len(data[i][2])
	for j in range(0, L, 3):
		c1 = data[i][2][j:j+3]
		c2 = b2[j:j+3]
		if c1 not in trimer[0]: trimer[0][c1] = 0
		if c2 not in trimer[1]: trimer[1][c2] = 0
		trimer[0][c1] += 1
		trimer[1][c2] += 1
		if j+6 < L:
			h1 = data[i][2][j:j+6]
			h2 = b2[j:j+6]
			if h1 not in hexamer[0]: hexamer[0][h1] = 0
			if h2 not in hexamer[1]: hexamer[1][h2] = 0
			hexamer[0][h1] += 1
			hexamer[1][h2] += 1
		
fh.close()

fh = open("../5000_ratio.txt", "w")
fh.write("\t".join(["# feature","position","sequence","cds","random","log2((cds+1)/(random+1))"]) + "\n")
for i in range(10):
	x = set(context[0][i].keys())
	x.update(context[1][i].keys())
	for k in sorted(x):
		c1, c2 = 0, 0
		if k in context[0][i]: c1 = context[0][i][k]
		if k in context[1][i]: c2 = context[1][i][k]
		v = round(math.log((c1+1.)/(c2+1.),2), 4)
		fh.write("\t".join(map(str,["context", i-6, k, c1, c2, v])) + "\n")
x = set(trimer[0].keys())
x.update(trimer[1].keys())
for k in sorted(x):
	c1, c2 = 0, 0
	if k in trimer[0]: c1 = trimer[0][k]
	if k in trimer[1]: c2 = trimer[1][k]
	v = round(math.log((c1+1.)/(c2+1.),2), 4)
	fh.write("\t".join(map(str,["trimer", "-", k, c1, c2, v])) + "\n")
x = set(hexamer[0].keys())
x.update(hexamer[1].keys())
for k in sorted(x):
	c1, c2 = 0, 0
	if k in hexamer[0]: c1 = hexamer[0][k]
	if k in hexamer[1]: c2 = hexamer[1][k]
	v = round(math.log((c1+1.)/(c2+1.),2), 4)
	fh.write("\t".join(map(str,["hexamer", "-", k, c1, c2, v])) + "\n")
fh.close()


