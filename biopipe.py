import os
import sys
import re
import sh
from temp import temp_file
import copy
import gzip

## TO DO:
# VCF input for plink
# VCFtools filters

class BioData(object):

	def __init__(self):
		self.cmd = None

	def execute(self):
		raise NotImplementedError

	def duplicate(self):
		return copy.deepcopy(self)

class VCFToolsData(BioData):

	def __init__(self, filename):
		self.filename = os.path.abspath(filename)
		self.out_prefix = None
		self.result = None

		if re.search("gz$", filename):
			self.cmd = sh.vcftools.bake("--gzvcf", filename)
		else:
			self.cmd = sh.vcftools.bake("--vcf", filename)

	## filters

	def min_DP(self, dp=5):
		self.cmd = self.cmd.bake("--min-meanDP", dp)
		return self

	def max_DP(self, dp=100):
		self.cmd = self.cmd.bake("--max-meanDP", dp)
		return self

	def min_Q(self, q=100):
		self.cmd = self.cmd.bake("--minQ", q)
		return self

	def maf(self, m=0.05):
		self.cmd = self.cmd.bake("--maf", m)
		return self

	## outputs

	def recode(self, out=None):
		if out is None:
			self.out_prefix = temp_file()
		else:
			self.out_prefix = out
		self.result = self.out_prefix + ".recode.vcf"

		self.cmd = self.cmd.bake("--recode")
		return self.execute()

	## executions

	def execute(self):
		if self.out_prefix is None:
			self.out_prefix = self.filename
		if self.result is None:
			self.result = self.filename

		self.cmd = self.cmd.bake("--out", self.out_prefix)
		self.cmd()
		return VCFToolsData(self.result)
		
	## converstions

	def to_plink(self):
		return PlinkData(self.filename)

	def to_vcf(self):
		return VCFData(self.filename)


class VCFData(BioData):

	def __init__(self, filename):
		self.filename = filename

		if re.search("gz$", filename):
			self.f = gzip.GzipFile(filename, "r")
		else:
			self.f = open(filename, "r")

	def snps_only(self, out=None):
		if out is None: out = temp_file() + ".vcf.gz"
		o = gzip.GzipFile(out, "w")

		bases = ["A", "T", "G", "C", "a", "t", "g", "c"]
		header = None
		for line in self.f:

			if re.search("##", line):
				o.write(line)
			elif re.search("#", line):
				o.write(line)
				header = re.split("\t", line)
			else:
				info = re.split("\t", line)
				if info[header.index("REF")] in bases:
					if info[header.index("ALT")] in bases:
						o.write(line)
		o.close()
		return VCFData(out)

	def reset_ids(self, out=None):
		if out is None: out = temp_file() + ".vcf.gz"
		o = gzip.GzipFile(out, "w")

		header = None
		for line in self.f:
			if re.search("##", line):
				o.write(line)
			elif re.search("#", line):
				o.write(line)
				header = re.split("\t", line)
			else:
				info = re.split("\t", line)
				new_id = (info[header.index("#CHROM")],
							info[header.index("POS")],
							info[header.index("REF")],
							info[header.index("ALT")])

				info[header.index("ID")] = "%s_%s_%s_%s" % new_id
				o.write("\t".join(info))

		o.close()

		return VCFData(out)

	## converstions

	def to_vcftools(self):
		return VCFToolsData(self.filename)

	def to_plink(self):
		return self.to_vcftools().to_plink()


class PlinkData(BioData):

	def __init__(self, filename):
		self.filename = os.path.abspath(filename)
		self.filename = os.path.realpath(self.filename)

		if os.path.exists(self.filename + ".bed"):
			self.cmd = sh.plink.bake("--bfile", self.filename)
		elif os.path.exists(self.filename + ".ped"):
			self.cmd = sh.plink.bake("--file", self.filename)
		elif re.search("vcf", filename):
			self.cmd = sh.plink.bake("--vcf", self.filename, 
									"--biallelic-only", "strict",
									"--vcf-half-call", "m",
									"--vcf-filter")
		else:
			raise IOError('Did not recognize Plink Data type')

	## file accessors

	def bed(self):
		return self._file_check("bed")

	def snplist(self):
		return self._file_check("snplist")

	def prune_in(self):
		return self._file_check("prune.in")

	def bim(self):
		return self._file_check("bim")

	def fam(self):
		return self._file_check("fam")

	def pca(self):
		return self._file_check("pca")

	def eigenvec(self):
		return self._file_check("eigenvec")

	def eigenval(self):
		return self._file_check("eigenval")

	def _file_check(self, ext):
		filename = "%s.%s" % (self.filename, ext)
		if not os.path.exists(filename): 
			raise IOError(".%s file not found" % ext)
		return filename

	## filtering

	def maf(self, maf=0.05):
		self.cmd = self.cmd.bake("--maf", maf)
		return self

	def chr(self, c):
		self.cmd = self.cmd.bake("--chr", c)
		return self

	def geno(self, geno=0.01):
		self.cmd = self.cmd.bake("--geno", geno)
		return self

	def keep(self, keep):
		self.cmd = self.cmd.bake("--keep", keep)
		return self

	def extract(self, extract):
		self.cmd = self.cmd.bake("--extract", extract)
		return self

	def remove(self, remove):
		self.cmd = self.cmd.bake("--remove", remove)
		return self

	def exclude(self, exclude):
		self.cmd = self.cmd.bake("--exclude", exclude)
		return self

	def pheno(self, pheno):
		self.cmd = self.cmd.pheno("--pheno", pheno)
		return self

	## reflexive output

	def pca(self):
		self.cmd = self.cmd.bake("--pca")
		return self.rexecute()

	def write_snplist(self):
		self.cmd = self.cmd.bake("--write-snplist")
		return self.rexecute()

	def indep_pairwise(self, window=100, step=10, r=0.1):
		self.cmd = self.cmd.bake("--indep-pairwise", 
								window,
								step,
								r)
		return self.rexecute()


	## transitive output

	def make_bed(self, out=None):
		self.cmd = self.cmd.bake("--make-bed")
		return self.execute(out)

	def indep_pairwise_filter(self, window=100, step=10, r=0.1, out=None):
		self.indep_pairwise(window, step, r)
		return self.extract(self.prune_in()).make_bed(out)

	def bmerge(self, d, out):
		self.cmd = self.cmd.bake("--bmerge", d.bed(), d.bim(), d.fam())
		return self.make_bed(out)

	## executions

	def execute(self, out=None):
		if out is None: out = temp_file()
		self.cmd("--out", out)
		return PlinkData(out)

	def rexecute(self):
		self.cmd("--out", self.filename)
		return PlinkData(self.filename)

	## conversions
	def to_admixture(self):
		return AdmixtureData(self.bed())


class AdmixtureData(BioData):
	
	def __init__(self, filename):
		self.filename = os.path.abspath(filename)
		self.base = re.sub("\..+$", "", self.filename)
		self.dir = os.path.dirname(self.filename)
		self.cmd = sh.admixture

	## file accessors

	def P(self):
		return self._file_check("P")

	def Q(self):
		return self._file_check("Q")

	def log(self):
		return self._file_check("log")

	def _file_check(self, ext):
		filename = "%s.%i.%s" % (self.base, self.k, ext)
		if not os.path.exists(filename):
			raise IOError("%s file not found" % ext)
		return filename

	## options

	def cv(self, cv=5):
		self.cmd = self.cmd.bake("--cv=%i" % cv)
		return self

	def threads(self, j=4):
		self.cmd = self.cmd.bake("-j%i" % j)
		return self

	## outputs

	def unsupervised(self, k):
		if isinstance(k, list):
			results = []
			for k_i in k: 
				results.append(self.duplicate().unsupervised(k_i))
			return results
		else:
			self.cmd = self.cmd.bake(self.filename, k)
			self.out = "%s.%i.log" % (self.base, k)
			print self.filename
			print self.out
			return self.execute()

	def supervised(self, pop):
		raise NotImplementedError("Supervised admixture not supported yet")

	## execute

	def execute(self):
		self.cmd(_cwd=self.dir, _out=self.out)
		return AdmixtureData(self.filename)

def inner_join(a, b, out=None):

	a_snps = open(a.write_snplist().snplist()).readlines()
	b_snps = open(b.write_snplist().snplist()).readlines()

	union_snps = temp_file()
	o = open(union_snps, "w")

	u = list(set(a_snps) & set(b_snps))
	for s in u:
		o.write(s)
	o.close()

	a_reduced = PlinkData(a.filename).extract(union_snps).make_bed()
	b_reduced = PlinkData(b.filename).extract(union_snps).make_bed()

	return a_reduced.bmerge(b_reduced, out)

def outer_join(a, b, out=None):
	return a.bmerge(b, out)

def multi_outer_join(data, out=None):
	start = data.pop()
	end = data.pop()
	mid = reduce(outer_join, data, start)

	outer_join(end, mid, out=out)

def multi_inner_join(data, out=None):
	start = data.pop()
	end = data.pop()
	mid = reduce(inner_join, data, start)

	outer_join(end, mid, out=out)

if __name__ == "__main__":
	pass