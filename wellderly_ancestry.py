from biopipe import *
import os
import itertools as it

WELLDERLY_DIR = "./data/wellderly/"
OTG_DIR = "./data/1000genomes/"

def preprocess(vcf):

	"""
	Uses biopipe to clean VCF. Steps are as follows:

		- keep only SNPs (discard complex variants)
		- reset variant IDs to consistent nomenclature
		- switch to using VCFtools
		- quality filter
		- generate new VCF

		- convert VCF to plink format
		- keep SNPs with call rate >90%
		- keeps SNPs with maf > 5%
		- create final bed file
	"""
	
	base_name = re.sub("\.vcf.*$", "", vcf)

	vcf_cleaned = VCFData(vcf). \
					snps_only(). \
					reset_ids(). \
					to_vcftools(). \
					min_Q(100). \
					recode()

	bed_cleaned = vcf_cleaned. \
			to_plink(). \
			make_bed(). \
			geno(0.1). \
			maf(0.05). \
			make_bed(base_name + "_cleaned")

	return bed_cleaned

def merge_vcfs(directory, out=None):
	""" 
	Preprocess and merge all VCF files in a directory
	"""

	files = get_raw_vcf_files(directory)
	cleaned = [preprocess(f) for f in files]

	merged = multi_outer_join(cleaned, out=out)
	return merged

def get_raw_vcf_files(directory):
	files = map(lambda x: directory + x, os.listdir(directory))
	files = filter(lambda x: "vcf" in x, files)
	files = filter(lambda x: "recode" not in x, files)
	return files

if __name__ == "__main__":

	# clean and merge 1000genomes VCF files
	otg = merge_vcfs(OTG_DIR, out="./data/1000genomes")
	well = merge_vcfs(WELLDERLY_DIR, out="./data/wellderly")

	# merge wellderly and 1000 genomes
	# filter for LD
	# run PCA
	data = inner_join(otg, well, out="./data/merged"). \
			indep_pairwise_filter(out="./data/merged_noLD"). \
			pca()

	# convert plink to admixture
	# use 8 threads
	# use 5-fold cross validation
	# run unsupervised admixture for 2-6 populations
	admixture = data. \
					to_admixture(). \
					threads(8). \
					cv(5). \
					unsupervised([2,3,4,5,6])




