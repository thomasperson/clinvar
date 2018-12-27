"""
Alternate implementation of master.py.  Removing pysam, pypez, and non clinvar related functions.  .

Run with -h to see all options.
"""


from datetime import datetime
import ftplib
import os
import sys
from distutils import spawn
sys.path.insert(0, 'src/')

import parse_clinvar_xml as pcx

try:
	import configargparse
	import pandas   # make sure all dependencies are installed
except ImportError as e:
	sys.exit("ERROR: Python module not installed. %s. Please run 'pip install -r requirements.txt' " % e)
for executable in ['wget', 'tabix', 'vt']:
	assert spawn.find_executable(executable), "Command %s not found, see README" % executable


def checkExists(fileAndPath):
	exists = os.path.isfile(fileAndPath)
	if exists:
		return True
	else:
		return False

def get_handle(path):
	if path[-3:] == '.gz':
		handle = gzip.open(path)
	else:
		handle = open(path)
	return handle

def parseArguments():
	p = configargparse.getArgParser()
	g = p.add_argument_group('main args')
	g.add("--b37-genome", help="GRCh37 .fa genome reference file", default=None, required=False, dest="b37fasta", type=str)
	g.add("--b38-genome", help="GRCh38 .fa genome reference file. NOTE: chromosome names must be like '1', '2'.. 'X', 'Y', 'MT'.", default=None, required=False, dest="b38fasta", type=str)
	g.add("-X", "--clinvar-xml", help="The local filename of the ClinVarFullRelase.xml.gz file. If not set, grab the latest from NCBI.", dest="xml_file", default="output_tmp/ClinVar.xml.gz", type=str)
	g.add("-S", "--clinvar-variant-summary-table", help="The local filename of the variant_summary.txt.gz file. If not set, grab the latest from NCBI.", dest="tsv_file", default="output_tmp/ClinVar.tsv.gz", type=str)
	g.add("-U", help='Download Latest ClinVar Files to default location', action='store_true', default=False, dest="download")
	g.add("--output-prefix", default="../output/", help="Final output files will have this prefix", type=str, dest='output_prefix')  #../output/ is directory not a prefix.
	g.add("--output-dir", default="../output/", help="Final output files will have this prefix", type=str, dest='output_dir')  #../output/ is directory not a prefix.
	g.add("--tmp-dir", default="./output_tmp/", help="Temporary output direcotry for temp files ", dest="output_tmp", type=str)
	g.add("--rm-temp", default=True, help="Removes toempoary directory and temp files when finished")


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

	os.system("mkdir -p " + cli_args.output_tmp)
	if not checkExists("src/normalize.py"):
		os.system("wget -O src/normalize.py https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py")


	if cli_args.xml_file is not None and not checkExists(cli_args.xml_file):
		print("No ClinVar XML file specified or specified file does not exist.  Downloadeding latest ClinVar XML file to use")
		os.system("wget -O output_tmp/ClinVar.xml.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz")
	elif cli_args.download:
		os.system("wget -O "+cli_args.output_tmp+"/ClinVar.xml.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz")
	if cli_args.tsv_file is not None and not checkExists(cli_args.tsv_file):
		print("No ClinVar TSV file specified or specified file does not exist.  Downloadeding latest ClinVar TSV file to use")
		os.system("wget -O "+cli_args.output_tmp+"/ClinVar.tsv.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz")
	elif cli_args.download:
		os.system("wget -O output_tmp/ClinVar.tsv.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz")

	############TODO!!!!!!  Currntly only does one genome build at a time.  Can multiprocess, but should run both at same time so file doesn't have to be parsed twice!!!  AND WHY is the file handls outside the method?
	################# WHAT IS THE POINT OF MULTI/SINGL ONLY?  That just skipps writing to file...

	dir_path = os.path.dirname(os.path.realpath(__file__))
	if cli_args.b37fasta is not None:
		pcx.parse_clinvar_tree(cli_args.xml_file, dir_path+cli_args.output_tmp, 'GRCh37')   #NOTE  parse_clinvar_xml.py uses findall pretty extensivly rather than relying structure of the xml....  need to check this for accuracy.
	if cli_args.b38fasta is not None:
		pcx.parse_clinvar_tree(cli_args.xml_file, cli_args.output_tmp, 'GRCh38')

	return


if __name__ == '__main__':
	main()
	exit()
