from biopipe import *
import os
import itertools as it
import sys

def merge_vcfs(directory, out=None):
	""" 
	Preprocess and merge all VCF files in a directory
	"""


	files = get_raw_vcf_files(directory)
	cleaned = [preprocess(f) for f in files]

	print "Cleaning finished\n\nBeginning merge\n"

	merged = multi_outer_join(cleaned, out=out)
	return merged

def get_raw_vcf_file(directory):
	files = map(lambda x: directory + x, os.listdir(directory))
	files = filter(lambda x: "vcf" in x, files)
	files = filter(lambda x: "recode" not in x, files)
	files = filter(lambda x: "cleaned" in x, files)
	return files

if __name__ == "__main__":
	