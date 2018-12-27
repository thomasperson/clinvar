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
	g.add("--b37-genome", help="GRCh37 .fa genome reference file", default=None, required=False, dest="b37fasta")
	g.add("--b38-genome", help="GRCh38 .fa genome reference file. NOTE: chromosome names must be like '1', '2'.. 'X', 'Y', 'MT'.", default=None, required=False, dest="b38fasta")
	g.add("-X", "--clinvar-xml", help="The local filename of the ClinVarFullRelase.xml.gz file. If not set, grab the latest from NCBI.", dest="xml_file", default="output_tmp/ClinVar.xml.gz")
	g.add("-S", "--clinvar-variant-summary-table", help="The local filename of the variant_summary.txt.gz file. If not set, grab the latest from NCBI.", dest="tsv_file", default="output_tmp/ClinVar.tsv.gz")
	g.add("--output-prefix", default="../output/", help="Final output files will have this prefix")
	g.add("--tmp-dir", default="./output_tmp", help="Temporary output files will have this prefix")
	g2 = p.add_mutually_exclusive_group() #WHY? DOES IT SAVE TIME?
	g2.add("--single-only", dest="single_or_multi", action="store_const", const="single", help="Only generate the single-variant tables")
	g2.add("--multi-only", dest="single_or_multi", action="store_const", const="multi", help="Only generate the multi-variant tables")

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

	os.system("mkdir -p " + cli_args.tmp_dir)
	if not checkExists("src/normalize.py"):
		os.system("wget -O src/normalize.py https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py")


	if cli_args.xml_file is not None and not checkExists(cli_args.xml_file):
		print("No ClinVar XML file specified or specified file does not exist.  Downloadeding latest ClinVar XML file to use")
		os.system("wget -O output_tmp/ClinVar.xml.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz")
	if cli_args.tsv_file is not None and not checkExists(cli_args.tsv_file):
		print("No ClinVar TSV file specified or specified file does not exist.  Downloadeding latest ClinVar TSV file to use")
		os.system("wget -O output_tmp/ClinVar.tsv.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz")

	############TODO!!!!!!  Currntly only does one genome build at a time.  Can multiprocess, but should run both at same time so file doesn't have to be parsed twice!!!  AND WHY is the file handls outside the method?
	if cli_args.b37fasta is not None:
		f = open(args.multi, 'w')
        parse_clinvar_tree(get_handle(args.xml_path), dest=args.out, multi=f, genome_build='GRCh37')
        f.close()
    if cli_args.b37fasta is not None:
        parse_clinvar_tree(get_handle(args.xml_path), dest=args.out, genome_build=args.genome_build)

	return


if __name__ == '__main__':
	main()
	exit()
