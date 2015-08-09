from biopipe import *

# convert plink to admixture
# use 8 threads
# use 5-fold cross validation
# run unsupervised admixture for 2-6 populations
admixture = PlinkData("./data/merged_noLD"). \
				to_admixture(). \
				threads(8). \
				cv(5). \
				unsupervised([2,3,4,5,6])