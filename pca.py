from biopipe import *

# merged data created by preprocess.py
# run LD filter
# run PCA

PlinkData("./data/merged"). \
	indep_pairwise_filter(out="./data/merged_noLD"). \
	pca()

