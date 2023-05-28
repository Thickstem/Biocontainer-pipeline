import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
import math


fh = open("../train_data.txt")
col = fh.readline().rstrip().split()
data, dataScaled = {}, {}
for i in range(len(col)):
	data[col[i]] = [] 
	dataScaled[col[i]] = [] 
for line in fh:
	for i,v in enumerate(map(float,line.rstrip().split())):
		data[col[i]].append(v)
		if i == 0:
			dataScaled[col[i]].append(int(v))
fh.close()

# feature values are scaled to a range of 0 to 1
for i in range(1, len(col)):
	maximum = max(data[col[i]])
	minimum = min(data[col[i]])
	for x in data[col[i]]:
		v = round((x-minimum)/(maximum-minimum), 3)
		dataScaled[col[i]].append(v)
fh = open("../train_data_scaled.txt","w")
fh.write("\t".join(col) + "\n")
for k in range(len(dataScaled[col[0]])):
	tmp = []
	for i in range(len(col)):
		tmp.append(dataScaled[col[i]][k])
	fh.write("\t".join(map(str,tmp)) + "\n")
fh.close()


highCorrelatedT = 0.8
fh = open("../redundant_features.txt","w")
fh.write("\t".join(["feature_1","feature_2(removed)","correlation(r)"]) + "\n")
dat2 = {}
for i in range(1, len(col)):
	dat2[col[i]] = dataScaled[col[i]]
df2 = pd.DataFrame.from_dict(dat2)
fig = plt.figure(figsize=(20,20), dpi=100)
sns.set_style("whitegrid")
df2_corr = df2.corr(method='pearson')
redundantFeatures = []
for i in range(1, len(col)):
	for j in range(i+1, len(col)):
		if df2_corr.at[col[i], col[j]] > highCorrelatedT or df2_corr.at[col[i], col[j]] < -highCorrelatedT:
			fh.write("\t".join(map(str, [col[i], col[j], round(df2_corr.at[col[i], col[j]],3)])) + "\n")
			redundantFeatures.append(col[j])
fh.close()
g = sns.heatmap(df2_corr , vmax=1, vmin=-1, square=True,xticklabels=True,yticklabels=True, cmap="RdBu_r")
g.set_yticklabels(g.get_xticklabels(), rotation =0)
g.set_xticklabels(g.get_xticklabels(), rotation =90)
plt.title("Correlations(r) of features")
fig.tight_layout()
plt.savefig("../feature_corr.png")
plt.close(fig)

dat3 = {}
for i in range(1, len(col)):
	if col[i] not in redundantFeatures:
		dat3[col[i]] = dataScaled[col[i]]
df3 = pd.DataFrame.from_dict(dat3)
fig = plt.figure(figsize=(20,20), dpi=100)
sns.set_style("whitegrid")
df3_corr = df3.corr(method='pearson')
g = sns.heatmap(df3_corr , vmax=1, vmin=-1, square=True,xticklabels=True,yticklabels=True,cmap="RdBu_r")
g.set_yticklabels(g.get_xticklabels(), rotation =0)
g.set_xticklabels(g.get_xticklabels(), rotation =90)
plt.title("Correlations(r) of low-redundant features")
fig.tight_layout()
plt.savefig("../feature_corr_low-redundant.png")
plt.close(fig)

