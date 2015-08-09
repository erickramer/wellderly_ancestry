import os
import re
import gzip

N_LINES = 100
TEMPDIR = "./data/temp"

WELLDERLYDIR = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_results"
OTGDIR = "/gpfs/home/ekramer/data/1000genomes/vcf"

if __name__ == "__main__":

	# copy wellderly
	wellderly_files = os.listdir(WELLDERLYDIR)
	for f in wellderly_files:		
		if re.search("chr[0-9]+", f):
			ch = re.search("chr[0-9]+", f).group(0)


		input_file = "%s/%s" % (WELLDERLYDIR, f)
		output_file = "%s/wellderly_%s.vcf" % (TEMPDIR, ch)
		o = open(output_file, "w")

		i = 0
		for line in file(input_file, "r"):
			o.write(line)
			i += 1
			if i > N_LINES: break

	# copy 1000genomes
	otg_files = os.listdir(OTGDIR)
	for f in otg_files:		
		if re.search("vcf.gz$", f):
			if re.search("chr", f):
				ch = re.search("chr[0-9]+", f).group(0)
				
				input_file = "%s/%s" % (OTGDIR, f)
				output_file = "%s/1000genomes_%s.vcf" % (TEMPDIR, ch)
				o = open(output_file, "w")

				i = 0
				for line in gzip.GzipFile(input_file, "r"):
					o.write(line)
					i += 1
					if i > N_LINES: break

