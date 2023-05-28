import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from itertools import cycle
import math
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_selection import RFE, f_regression
from sklearn.linear_model import (LinearRegression, Ridge, Lasso, RandomizedLasso, LogisticRegression, LogisticRegressionCV)
from sklearn.ensemble import (RandomForestRegressor,GradientBoostingRegressor)
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics as mtr
from sklearn.metrics import (accuracy_score,roc_curve, auc)
from sklearn.svm import l1_min_c

def rescale(data):
    minmax = MinMaxScaler()
    vals = minmax.fit_transform(np.array([data]).T).T[0]
    vals = map(lambda x: round(x,2), vals)
    return vals

def make_kde(*args, **kwargs):
	sns.kdeplot(*args, cmap=next(make_kde.cmap_cycle), **kwargs)

def annotate_colname(x, **kws):
  ax = plt.gca()
  ax.annotate(x.name, xy=(0.05, 0.9), xycoords=ax.transAxes, fontweight='bold')

def main(randNum, vertical):

	redundantFeatures = []
	fh = open("../redundant_features.txt")
	fh.next()
	for line in fh:
		redundantFeatures.append(line.split()[1])
	fh.close()

	fh = open("../train_data_scaled.txt")
	feature = fh.readline().rstrip().split()[1:]
	X, y = [], []
	for line in fh:
		a = line.rstrip().split()
		y.append(int(a[0]))
		X.append(map(float,a[1:]))
	fh.close()

	featureLowRedundant = []
	idx = []
	for i,x in enumerate(feature):
		if x not in redundantFeatures:
			idx.append(i)
			featureLowRedundant.append(x)
	XlowRedundant = []
	for x in X:
		XlowRedundant.append([x[i] for i in idx])

	fh = open("../nonredundant_features.txt", "w")
	for x in featureLowRedundant:
		fh.write(x + "\n")
	fh.close()

	x_train, x_test, y_train, y_test = train_test_split(XlowRedundant, y, test_size=0.2, random_state=randNum)

	coef = []	
	acc = []
	aauc = []
	coef2 = [ [] for i in range(len(featureLowRedundant))]

	C = [x/1000. for x in range(9, 1001, 1)]
	model = LogisticRegression(penalty="l1", C=1.)
	for c in C:
		model.set_params(C=c)
		model.fit(x_train, y_train)
		coef.append([round(x,3) for x in model.coef_[0]])
		y_predict = model.predict(x_test)
		acc.append(accuracy_score(y_test, y_predict))
		fpr, tpr, thresholds = roc_curve(y_test, y_predict)
		roc_auc = auc(np.array(fpr), np.array(tpr))
		aauc.append(roc_auc)
		for i,x in enumerate(model.coef_[0]):
			coef2[i].append(round(x,3))
	fh = open("../logistic_l1_ft_coef.txt", "w")
	tmp = ["Accuracy","C"]
	tmp.extend(featureLowRedundant)
	fh.write("\t".join(map(str, tmp)) + "\n")
	for i in range(len(C)):
		tmp = [acc[i],C[i]]
		tmp.extend(coef[i])
		fh.write("\t".join(map(str,tmp)) + "\n")
	fh.close()

	fig, ax2 = plt.subplots(figsize=(15,5.5))
	ax2.plot(np.log10(np.array(C)), acc, linestyle="dashed", color="b")
	ax2.set_ylabel('Accuracy on testing dataset', color="blue")
	ax2.set_xlabel(r"log$_{10}$C")
	plt.tick_params(axis="y", labelcolor="b")

	ax1 = ax2.twinx( )
	b = []
	for i, f in enumerate(featureLowRedundant):
		n = 0
		for x in coef2[i]:
			if x != 0:
				n+=1
		b.append(n)

	fNum = 0
	fCoef = {}
	fCoef2 = {}
	verticalACC = 0
	for i in range(len(C)):
		if C[i] == vertical:
			verticalACC = acc[i]
			for k,j in enumerate(coef[i]):
				if abs(j) > 0:
					fCoef[featureLowRedundant[k]] = j
					fNum += 1
		if i == len(C)-1:
			verticalACC = acc[i]
			for k,j in enumerate(coef[i]):
				if abs(j) > 0:
					fCoef2[featureLowRedundant[k]] = j

	n = 0
	bb = np.argsort(b)[::-1]
	for i in bb:
		if b[i] == 0: continue
		f = featureLowRedundant[i]
		n += 1

	fRank2 = sorted(fCoef2, key=lambda k: fCoef2[k])[::-1]
	for x in fRank2:
		i = 0
		for y in featureLowRedundant:
			if y == x: break
			i += 1
		f = featureLowRedundant[i]
		ax1.plot(np.log10(np.array(C)), coef2[i], label=f)

	ax1.plot([], [], label="Accuracy", color="b", linestyle="dashed")
	ax1.set_ylabel("Coefficients")
	ax1.tick_params('y', colors='black')
	plt.title("Feature selection by $\mathcal{L}1$-logistic regression")
	plt.legend(ncol=2, bbox_to_anchor=(1.08, 1), loc=2, borderaxespad=0.)
	plt.axvline(x=math.log(vertical,10), color="black", ls="dashed")
	

	fRank = sorted(fCoef, key=lambda k: abs(fCoef[k]))[::-1]
	tmp = ["Selected %d features are sorted by importance (|coefficient|):"%fNum]
	for i, x in enumerate(fRank):
		tmp.append("%d, %s (%.3f)"%(i+1, x, fCoef[x]))
	ymin, ymax = ax1.get_ylim()
	xmin, xmax = ax1.get_xlim()
	plt.text(math.log(vertical,10), ymin, "%d features with |coefficient|>0 \n(C=%.3f, acc=%.3f) "%(fNum, vertical, verticalACC),
		horizontalalignment='right', verticalalignment='bottom')
	plt.text(xmin+0.02, ymax-0.1, "\n".join(tmp), horizontalalignment='left', verticalalignment='top')
	plt.tight_layout()
	plt.savefig("../logistic_l1_path.png")
	plt.close(fig)
	



main(1, 0.204)









