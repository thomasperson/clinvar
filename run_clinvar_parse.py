"""
Alternate implementation of master.py.  Removing pysam, pypez, and non clinvar related functions.
Making as python as can so can run on windows and *nix based systems.  Making complient for
python3 as well as python2

Run with -h to see all options.
"""


import os
import sys
from distutils import spawn
import shutil
import operator

sys.path.insert(0, 'src'+os.sep)

import parse_clinvar_xml as pcx

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
for executable in ['wget', 'tabix', 'vt', 'bgzip']:  #working on pure python implementation to remove need for these *nix based programs
	assert spawn.find_executable(executable), "Command %s not found, see README" % executable


def rreplace(s, old, new):
	li = s.rsplit(old, 1) #Split only once
	return new.join(li)


def checkExists(fileAndPath):
	exists = os.path.isfile(fileAndPath)
	if exists:
		return True
	else:
		return False

def sortRawTextFile(unSortedFile):
	to_sort=open(unSortedFile,'r')
	sorted_file=open(rreplace(unSortedFile,".",".sorted."),'wb')
	header=0
	file_lines=[]
	for line in to_sort:
		if header==0:
			sorted_file.write(line)
			header+=1
		else:
			file_lines.append([x.strip() for x in line.split('\t')])
	to_sort.close()

	for line in sorted(file_lines, key=operator.itemgetter(8)):
		sorted_file.write("\t".join(line)+"\n")

	sorted_file.close()
	return rreplace(unSortedFile,".",".sorted.")


def download_file(url,local_filename):
	response = urlopen(url)
	local_file=open(local_filename, 'wb')
	local_file.write(response.read())
	local_file.close()
	return

def parseArguments():
	dir_path = os.path.dirname(os.path.realpath(__file__))
	p = configargparse.getArgParser()
	g = p.add_argument_group('main args')
	g.add("--b37-genome", help="GRCh37 .fa genome reference file", default=None, required=False, dest="b37fasta", type=str)
	g.add("--b38-genome", help="GRCh38 .fa genome reference file. NOTE: chromosome names must be like '1', '2'.. 'X', 'Y', 'MT'.", default=None, required=False, dest="b38fasta", type=str)
	g.add("-X", "--clinvar-xml", help="The local filename of the ClinVarFullRelase.xml.gz file. If not set, grab the latest from NCBI.", dest="xml_file", default=dir_path+os.sep+"output_tmp"+os.sep+"ClinVar.xml.gz", type=str)
	g.add("-S", "--clinvar-variant-summary-table", help="The local filename of the variant_summary.txt.gz file. If not set, grab the latest from NCBI.", dest="tsv_file", default=dir_path+os.sep+"output_tmp"+os.sep+"ClinVar.tsv.gz", type=str)
	g.add("-N", help='Download Latest ClinVar Files to default location and run pipeline', action='store_true', default=False, dest="download_new")
	g.add("--output-prefix", default="clinvar_alleles", help="Final output files will have this prefix", type=str, dest='output_prefix')  #../output/ is directory not a prefix.
	g.add("--output-dir", default=dir_path+os.sep+"output"+os.sep, help="Final output files will be located here", type=str, dest='output_dir')  #../output/ is directory not a prefix.
	g.add("--tmp-dir", default=dir_path+os.sep+"output_tmp"+os.sep, help="Temporary output direcotry for temp files ", dest="output_tmp", type=str)
	g.add("--rm-temp", default=True, help="Removes tempoary directories and temp files when finished.  Setting flag will cause temp files to not be removed.")


	return p.parse_args()

def main():

	cli_args=parseArguments()
	print(cli_args)

	if cli_args.b37fasta is None and cli_args.b38fasta is None:
		sys.exit("At least one genome reference file is required")
	else:
		if cli_args.b37fasta is not None and not checkExists(cli_args.b37fasta):
			sys.exit("Genome reference: file not found:\t"+cli_args.b37fasta)
		if cli_args.b38fasta is not None and not checkExists(cli_args.b38fasta):
			sys.exit("Genome reference: file not found:\t"+cli_args.b38fasta)

	if cli_args.download_new:
		shutil.rmtree(cli_args.output_tmp)
		os.makedirs(cli_args.output_tmp)

	elif not os.path.exists(cli_args.output_tmp):

		os.makedirs(cli_args.output_tmp)

	normalize_py="https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py"
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
		if cli_args.xml_file is not None and not checkExists(cli_args.xml_file):
			print("No ClinVar XML file specified or specified file does not exist.  Downloading latest ClinVar XML file to use")
			download_file(clinvar_xml, cli_args.output_tmp+"ClinVar.xml.gz")
			print("Downloading ClinVarFullRelease_00-latest.xml.gz complete")

		if cli_args.tsv_file is not None and not checkExists(cli_args.tsv_file):
			print("No ClinVar TSV file specified or specified file does not exist.  Downloading latest ClinVar TSV file to use")
			download_file(clinvar_tsv, cli_args.output_tmp+"ClinVar.tsv.gz")
			print("Downloading variant_summary.txt.gz complete")

		if not checkExists("src"+os.sep+"normalize.py"):
			print("Downloading normalize.py")
			download_file(normalize_py, "src"+os.sep+"normalize.py")
			print("Downloading normalize.py complete")


	############TODO!!!!!!  Currntly only does one genome build at a time.  Can multiprocess, but should run both at same time so file doesn't have to be parsed twice!!!  Just select which you want to output!

	if cli_args.b37fasta is not None:
		if not checkExists(cli_args.output_tmp+"clinvar_table_raw.single.GRCh37.tsv") or cli_args.download_new:
			pcx.parse_clinvar_tree(cli_args.xml_file, cli_args.output_tmp, 'GRCh37')   #NOTE  parse_clinvar_xml.py uses findall pretty extensivly rather than relying structure of the xml....  need to check this for accuracy.
	if cli_args.b38fasta is not None:
		if not checkExists(cli_args.output_tmp+"clinvar_table_raw.single.GRCh38.tsv") or cli_args.download_new:
			pcx.parse_clinvar_tree(cli_args.xml_file, cli_args.output_tmp, 'GRCh38')   ##NOTE  ALSO, skipped sequences...  check other sequence locations?!?!  as each record can store multiple places!
	sortRawTextFile(cli_args.output_tmp+"clinvar_table_raw.multi.GRCh38.tsv")
	sortRawTextFile(cli_args.output_tmp+"clinvar_table_raw.single.GRCh38.tsv")
	return


if __name__ == '__main__':
	main()
	exit()
