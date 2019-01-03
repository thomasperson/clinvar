"""
Alternate implementation of master.py.  Removing pysam, pypez, and non clinvar related functions.
Making as python as can so can run on windows and *nix based systems.  Making complient for
python3 as well as python2

Run with -h to see all options.
"""


import os
import sys
from distutils import spawn
from distutils.dir_util import mkpath
import shutil
import operator

sys.path.insert(0, 'src'+os.sep)

import parse_clinvar_xml as pcx
import group_by_allele as gba


try:
	# For Python 3.0 and later
	from urllib.request import urlopen
except ImportError:
	# Fall back to Python 2's urllib2
	from urllib2 import urlopen

try:
	import configargparse
	import pandas   # make sure all dependencies are installed
except ImportError as e:
	sys.exit("ERROR: Python module not installed. %s. Please run 'pip install -r requirements.txt' " % e)

for executable in ['tabix', 'vt', 'bgzip']:  #NOTE:working on pure python implementation to remove need for these *nix based programs
	assert spawn.find_executable(executable), "Command %s not found, see README" % executable

pysam_installed=False
try:
	import pysam
	import normalize
	pysam_installed=True
except ImportError as e:
	print("ERROR: Python module pysam not installed. normalize.py will not be run." )  ###NOTE:  normalize.py uses pysam for fasta access.  Maybe I can pull out just that portion or reimplement or find another python module.

def rreplace(s, old, new):
	li = s.rsplit(old, 1) #Split only once
	return new.join(li)

def checkExists(fileAndPath):
	exists = os.path.isfile(fileAndPath)
	if exists:
		return True
	else:
		return False

def download_file(url,local_filename):
	response = urlopen(url)
	local_file=open(local_filename, 'wb')
	local_file.write(response.read())
	local_file.close()
	return

def downloadClinVarFiles(cli_args):
	normalize_py="https://raw.githubusercontent.com/thomasperson/minimal_representation/master/normalize.py"
	clinvar_xml="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz"
	clinvar_tsv="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

	if cli_args.download_new:
		print("Downloading normalize.py")
		download_file(normalize_py, "src"+os.sep+"normalize.py")
		print("Downloading normalize.py complete")

		print("Downloading "+  clinvar_xml)
		download_file(clinvar_xml, cli_args.output_tmp+"ClinVar.xml.gz")
		print("Downloading of ClinVarFullRelease_00-latest.xml.gz complete")

		print("Downloading " + clinvar_tsv)
		download_file(clinvar_tsv, cli_args.output_tmp+"ClinVar.tsv.gz")
		print("Downloading variant_summary.txt.gz complete")

	else:
		if not checkExists(cli_args.xml_file):
			print("No ClinVar XML file specified or specified file does not exist.  Downloading latest ClinVar XML file to use")
			download_file(clinvar_xml, cli_args.output_tmp+"ClinVar.xml.gz")
			print("Downloading ClinVarFullRelease_00-latest.xml.gz complete")

		if not checkExists(cli_args.tsv_file):
			print("No ClinVar TSV file specified or specified file does not exist.  Downloading latest ClinVar TSV file to use")
			download_file(clinvar_tsv, cli_args.output_tmp+"ClinVar.tsv.gz")
			print("Downloading variant_summary.txt.gz complete")

		if not checkExists("src"+os.sep+"normalize.py"):
			print("Downloading normalize.py")
			download_file(normalize_py, "src"+os.sep+"normalize.py")
			print("Downloading normalize.py complete")
	return

def sortRawTextFile(unSortedFile, sortedFile):
	to_sort=open(unSortedFile,'r')
	sorted_file=open(sortedFile,'w')
	header=0
	file_lines=[]
	for line in to_sort:
		if header==0:
			sorted_file.write(line)
			header+=1
		elif line.strip()=="":
			continue
		else:
			file_lines.append([x.strip() for x in line.split('\t')])
	to_sort.close()

	for line in sorted(file_lines, key=operator.itemgetter(8)):
		sorted_file.write("\t".join(line)+"\n")

	sorted_file.close()
	return

def createDirectories(cli_args):

	if cli_args.download_new:
		shutil.rmtree(cli_args.output_tmp)

	mkpath(cli_args.output_tmp)

	if cli_args.b38fasta is not None:
		mkpath(cli_args.output_dir+os.sep+"GRCh38"+os.sep+"multi")
		mkpath(cli_args.output_dir+os.sep+"GRCh38"+os.sep+"single")
	if cli_args.b37fasta is not None:
		mkpath(cli_args.output_dir+os.sep+"GRCh37"+os.sep+"multi")
		mkpath(cli_args.output_dir+os.sep+"GRCh37"+os.sep+"single")

	return

def parseArguments():
	dir_path = os.path.dirname(os.path.realpath(__file__))
	p = configargparse.getArgParser()
	g = p.add_argument_group('main args')
	g.add("--b37-genome", help="GRCh37 .fa genome reference file", default=None, required=False, dest="b37fasta", type=str)   ##REQUIRED for normalization, but normalization not requried to run script!
	g.add("--b38-genome", help="GRCh38 .fa genome reference file. NOTE: chromosome names must be like '1', '2'.. 'X', 'Y', 'MT'.", default=None, required=False, dest="b38fasta", type=str)  ##REQUIRED for normalization, but normalization not requried to run script!
	g.add("-X", "--clinvar-xml", help="The local filename of the ClinVarFullRelase.xml.gz file. If not set, grab the latest from NCBI.", dest="xml_file", default=dir_path+os.sep+"output_tmp"+os.sep+"ClinVar.xml.gz", type=str)
	g.add("-S", "--clinvar-variant-summary-table", help="The local filename of the variant_summary.txt.gz file. If not set, grab the latest from NCBI.", dest="tsv_file", default=dir_path+os.sep+"output_tmp"+os.sep+"ClinVar.tsv.gz", type=str)
	g.add("-N", "--new", help='Download all New.  Causes ouput_tmp to be removed and recreated and latest ClinVar Files to be downloaded', action='store_true', default=False, dest="download_new")
	g.add("--output-prefix", help="Final output files will have this prefix", type=str, dest='output_prefix',default="")
	g.add("--output-dir", default=dir_path+os.sep+"output"+os.sep, help="Final output files will be located here", type=str, dest='output_dir')  #../output/ is directory not a prefix.
	g.add("--tmp-dir", default=dir_path+os.sep+"output_tmp"+os.sep, help="Temporary output direcotry for temp files ", dest="output_tmp", type=str)
	g.add("--rm-temp", default=True, help="Removes tempoary directories and temp files when finished.  Setting flag will cause temp files to NOT be removed. ")

	return p.parse_args()

def main():

	cli_args=parseArguments()
	print(cli_args)
	fasta_files=False
	if cli_args.b37fasta is None and cli_args.b38fasta is None:
		print("Genome reference files are required for normalization.  Normalization will be skipped.")
	else:   ####TODO: AUTOMATIC FASTA DOWNLOAD?
		fasta_files=True
		if cli_args.b37fasta is not None and not checkExists(cli_args.b37fasta):
			sys.exit("Genome reference: file not found:\t"+cli_args.b37fasta)
		if cli_args.b38fasta is not None and not checkExists(cli_args.b38fasta):
			sys.exit("Genome reference: file not found:\t"+cli_args.b38fasta)

	#creates temp and outupt directories
	createDirectories(cli_args)

	#downloads clinvar and normalize_py
	downloadClinVarFiles(cli_args)


	############TODO!!!!!!  Currntly only does one genome build at a time.  Can multiprocess, but should run both at same time so file doesn't have to be parsed twice,
	###						 if you want both!!!  Just select which you want to output!

	###NOTE  cut -f 9 clinvar_table_raw.single.GRCh38.sorted.tsv | sort -u | wc -l == cut -f 1-4 clinvar_table_raw.single.GRCh38.sorted.tsv | sort -u | wc -l
	###		SO, Will run normalize later or not at all as only being sorted is what is required for group_by_allele.  MIGHT end up required after checking missing.
	### 	ALSO means fasta file NOT REQUIRED!

	####NOTE Multi essentially done after parse and sort...  Grouping by allele doens't make much sense when you have different sites, spot check manual review seems like lots of spam too.  Need a new grouping stratagy for multi.
	####QUESTION: How many Multi show up only in multi and havie higher/lower stars?
	####NOTE gonna normalize or not multi and call done.
	####NOTE There are missing variants from the multi.
	####NOTE Multi stars are comming from tsv so a bit screwy, with indiviaul aleles contributing more.  Doubtful can get a correct star unless in xml.


	if cli_args.b37fasta is not None:
		pcx.parse_clinvar_tree(cli_args.xml_file, cli_args.output_tmp, 'GRCh37')   #NOTE  parse_clinvar_xml.py uses findall pretty extensivly rather than relying structure of the xml....  need to check this for accuracy.
		sortRawTextFile(cli_args.output_tmp+"clinvar_table_raw.single.GRCh37.tsv",cli_args.output_tmp+"clinvar_table_raw.single.GRCh37.sorted.tsv")
		sortRawTextFile(cli_args.output_tmp+"clinvar_table_raw.multi.GRCh37.tsv",cli_args.output_tmp+"clinvar_table_raw.multi.GRCh37.sorted.tsv")
		if pysam_installed and fasta_files:
			normalize.normalize_tab_delimited_file(open(cli_args.output_tmp+"clinvar_table_raw.multi.GRCh37.sorted.tsv",'r'),open(cli_args.output_dir+"GRCh38"+os.sep+"multi"+os.sep+cli_args.output_prefix+"clinvar_multi_allele_haplotype.GRCh37.tsv",'w'),cli_args.b38fasta)
			normalize.normalize_tab_delimited_file(open(cli_args.output_tmp+"clinvar_table_raw.single.GRCh37.sorted.tsv",'r'),open(cli_args.output_tmp+"clinvar_table_raw.single.GRCh37.sorted.norm.tsv",'w'),cli_args.b37fasta)
			gba.group_by_allele(cli_args.output_tmp+"sorted.clinvar_table_raw.single.GRCh37.sorted.norm.tsv", cli_args.output_tmp+"clinvar_alleles_grouped.single.GRCh37.tsv")
		else:
			shutil.copyfile(cli_args.output_tmp+"clinvar_table_raw.multi.GRCh37.sorted.tsv",cli_args.output_dir+"GRCh37"+os.sep+"multi"+os.sep+cli_args.output_prefix+"clinvar_multi_allele_haplotype.GRCh37.tsv")
		gba.group_by_allele(cli_args.output_tmp+"sorted.clinvar_table_raw.single.GRCh37.sorted.tsv", cli_args.output_tmp+"clinvar_alleles_grouped.single.GRCh37.tsv")

	if cli_args.b38fasta is not None:
		#pcx.parse_clinvar_tree(cli_args.xml_file, cli_args.output_tmp, 'GRCh38')   ##NOTE  ALSO, skipped sequences...  check other sequence locations?!?!  as each record can store multiple places!
		sortRawTextFile(cli_args.output_tmp+"clinvar_table_raw.single.GRCh38.tsv",cli_args.output_tmp+"clinvar_table_raw.single.GRCh38.sorted.tsv")
		sortRawTextFile(cli_args.output_tmp+"clinvar_table_raw.multi.GRCh38.tsv",cli_args.output_tmp+"clinvar_table_raw.multi.GRCh38.sorted.tsv")
		if pysam_installed and fasta_files:
			normalize.normalize_tab_delimited_file(open(cli_args.output_tmp+"clinvar_table_raw.multi.GRCh38.sorted.tsv",'r'),open(cli_args.output_dir+"GRCh38"+os.sep+"multi"+os.sep+cli_args.output_prefix+"clinvar_multi_allele_haplotype.GRCh38.tsv",'w'),cli_args.b38fasta)
			normalize.normalize_tab_delimited_file(open(cli_args.output_tmp+"clinvar_table_raw.single.GRCh38.sorted.tsv",'r'),open(cli_args.output_tmp+"clinvar_table_raw.single.GRCh38.sorted.norm.tsv",'w'),cli_args.b37fasta)
			gba.group_by_allele(cli_args.output_tmp+"clinvar_table_raw.single.GRCh38.sorted.norm.tsv", ccli_args.output_tmp+"clinvar_alleles_grouped.single.GRCh38.tsv")
		else:
			shutil.copyfile(cli_args.output_tmp+"clinvar_table_raw.multi.GRCh38.sorted.tsv",cli_args.output_dir+"GRCh38"+os.sep+"multi"+os.sep+cli_args.output_prefix+"clinvar_multi_allele_haplotype.GRCh38.tsv")
		gba.group_by_allele(cli_args.output_tmp+"clinvar_table_raw.single.GRCh38.sorted.tsv", cli_args.output_tmp+"clinvar_alleles_grouped.single.GRCh38.tsv")






	return


if __name__ == '__main__':
	main()
	exit()
