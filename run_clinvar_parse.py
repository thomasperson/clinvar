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
import join_variant_summary_with_clinvar_alleles as isec
import clinvar_table_to_vcf as vcf


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

for executable in ['tabix', 'bgzip', 'vcf-sort']:
	assert spawn.find_executable(executable), "Command %s not found, see README" % executable

pysam_installed=False
try:
	import pysam
	import normalize
	pysam_installed=True
except ImportError as e:
	print("ERROR: Python module pysam not installed. normalize.py will not be run." )  ###NOTE:  normalize.py uses pysam for fasta access. TODO: Maybe I can pull out just that portion or reimplement or find another python module.


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
	clinvar_xml="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz"
	clinvar_variant_summary_tsv="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
	clinvar_submission_summary_tsv="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"

	if cli_args.download_new:
		if cli_args.run_xml:
			print("Downloading "+  clinvar_xml)
			download_file(clinvar_xml, cli_args.xml_file)
			print("Downloading of ClinVarFullRelease_00-latest.xml.gz complete")
		if cli_args.run_tsv:
			print("Downloading "+  clinvar_submission_summary_tsv)
			download_file(clinvar_submission_summary_tsv, cli_args.S_tsv_file)
			print("Downloading of submission_summary.txt.gz complete")
		print("Downloading " + clinvar_variant_summary_tsv)
		download_file(clinvar_variant_summary_tsv, cli_args.V_tsv_file)
		print("Downloading of variant_summary.txt.gz complete")

	else:
		if not checkExists(cli_args.xml_file) and cli_args.run_xml:
			print("No ClinVar XML file specified or specified file does not exist.  Downloading latest file to use")
			download_file(clinvar_xml, cli_args.output_tmp+"ClinVarFullRelease_00-latest.xml.gz")
			print("Downloading ClinVarFullRelease_00-latest.xml.gz complete")

		if not checkExists(cli_args.V_tsv_file):
			print("No ClinVar variant_summary file specified or specified file does not exist.  Downloading latest file to use")
			download_file(clinvar_variant_summary_tsv, cli_args.V_tsv_file)
			print("Downloading variant_summary.txt.gz complete")

		if not checkExists(cli_args.S_tsv_file) and cli_args.run_tsv:
			print("No ClinVar submission_summary file specified or specified file does not exist.  Downloading latest file to use")
			download_file(clinvar_submission_Reading & Writing GZIP Files in Pythonsummary_tsv, cli_args.S_tsv_file)
			print("Downloading variant_summary.txt.gz complete")

	return

def sortRawXMLoutTextFile(unSortedFile, sortedFile):
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

	for line in sorted(file_lines, key=operator.itemgetter(8)):   #Sorts on variant ID rather than genomic positon.
		sorted_file.write("\t".join(line)+"\n")

	sorted_file.close()
	return

def createDirectories(cli_args):

	if cli_args.download_new and os.path.exists(cli_args.output_tmp):
		shutil.rmtree(cli_args.output_tmp)
	if cli_args.download_new and os.path.exists(cli_args.output_dir) and cli_args.rm_output:
		shutil.rmtree(cli_args.output_dir)

	mkpath(cli_args.output_tmp)
	if cli_args.b38fasta is not None:
		mkpath(cli_args.output_dir+os.sep+"GRCh38"+os.sep+"multi")
		mkpath(cli_args.output_dir+os.sep+"GRCh38"+os.sep+"single")
	if cli_args.b37fasta is not None:
		mkpath(cli_args.output_dir+os.sep+"GRCh37"+os.sep+"multi")
		mkpath(cli_args.output_dir+os.sep+"GRCh37"+os.sep+"single")

	return

def runXMLpipeLine(cli_args, genome_build_id,fasta):
	pcx.parse_clinvar_tree(cli_args.xml_file, cli_args.output_tmp, genome_build_id)   #NOTE  parse_clinvar_xml.py uses findall pretty extensivly rather than relying structure of the xml....  need to check this for accuracy.
	print("Sorting ouput Files:")
	sortRawXMLoutTextFile(cli_args.output_tmp+"clinvar_table_raw.single."+genome_build_id+".tsv",cli_args.output_tmp+"clinvar_table_raw.single."+genome_build_id+".sorted.tsv")
	sortRawXMLoutTextFile(cli_args.output_tmp+"clinvar_table_raw.multi."+genome_build_id+".tsv",cli_args.output_tmp+"clinvar_table_raw.multi."+genome_build_id+".sorted.tsv")
	if pysam_installed:
		normalize.normalize_tab_delimited_file(cli_args.output_tmp+"clinvar_table_raw.multi."+genome_build_id+".sorted.tsv",cli_args.output_dir+""+genome_build_id+""+os.sep+"multi"+os.sep+cli_args.output_prefix+"clinvar_multi_allele_haplotype."+genome_build_id+".tsv",fasta)
		normalize.normalize_tab_delimited_file(cli_args.output_tmp+"clinvar_table_raw.single."+genome_build_id+".sorted.tsv",cli_args.output_tmp+"clinvar_table_raw.single."+genome_build_id+".sorted.norm.tsv",fasta)
		gba.group_by_allele(cli_args.output_tmp+"clinvar_table_raw.single."+genome_build_id+".sorted.norm.tsv", cli_args.output_tmp+"clinvar_alleles_grouped.single."+genome_build_id+".tsv")
		pass
	else:
		shutil.copyfile(cli_args.output_tmp+"clinvar_table_raw.multi."+genome_build_id+".sorted.tsv",cli_args.output_dir+genome_build_id+os.sep+"multi"+os.sep+cli_args.output_prefix+"clinvar_multi_allele_haplotype."+genome_build_id+".tsv")
		gba.group_by_allele(cli_args.output_tmp+"clinvar_table_raw.single."+genome_build_id+".sorted.tsv", cli_args.output_tmp+"clinvar_alleles_grouped.single."+genome_build_id+".tsv")
	isec.join_variant_summary_with_clinvar_alleles(cli_args.V_tsv_file, cli_args.output_tmp+"clinvar_alleles_grouped.single."+genome_build_id+".tsv", genome_build_id,cli_args.output_dir+genome_build_id+os.sep+"single"+os.sep+cli_args.output_prefix+"clinvar_allele_trait_pairs.single."+genome_build_id+".tsv.gz")
	if spawn.find_executable('vcf-sort') is not None:
		vcf.table_to_vcf(cli_args.output_dir+genome_build_id+os.sep+"single"+os.sep+cli_args.output_prefix+"clinvar_allele_trait_pairs.single."+genome_build_id+".tsv.gz", cli_args.b38fasta, cli_args.output_tmp+"clinvar_allele_trait_pairs.single."+genome_build_id+".unsorted.vcf")
		os.system("vcf-sort "+cli_args.output_tmp+"clinvar_allele_trait_pairs.single."+genome_build_id+".unsorted.vcf > "+cli_args.output_dir+genome_build_id+os.sep+"single"+os.sep+cli_args.output_prefix+"clinvar_allele_trait_pairs.single."+genome_build_id+".sorted.vcf")
		if spawn.find_executable('bgzip') is not None:
			os.system("bgzip -f "+cli_args.output_dir+genome_build_id+os.sep+"single"+os.sep+cli_args.output_prefix+"clinvar_allele_trait_pairs.single."+genome_build_id+".sorted.vcf")
			if spawn.find_executable('tabix') is not None:
				os.system("tabix -p vcf -f "+cli_args.output_dir+genome_build_id+os.sep+"single"+os.sep+cli_args.output_prefix+"clinvar_allele_trait_pairs.single."+genome_build_id+".sorted.vcf.gz")
	else:
		vcf.table_to_vcf(cli_args.output_dir+genome_build_id+os.sep+"single"+os.sep+cli_args.output_prefix+"clinvar_allele_trait_pairs.single."+genome_build_id+".tsv.gz", cli_args.b37fasta, cli_args.output_dir+genome_build_id+os.sep+"single"+os.sep+cli_args.output_prefix+"clinvar_allele_trait_pairs.single."+genome_build_id+".unsorted.vcf")

	return

def parseArguments():
	dir_path = os.path.dirname(os.path.realpath(__file__))
	p = configargparse.getArgParser()
	g = p.add_argument_group('main args')
	g.add("--b37-genome", help="GRCh37 .fa genome reference file  NOTE: fai index also required. ", default=None, required=False, dest="b37fasta", type=str)   ##REQUIRED for normalization, but normalization not requried to run script!
	g.add("--b38-genome", help="GRCh38 .fa genome reference file. NOTE: Chromosome names must be like '1', '2'.. 'X', 'Y', 'MT'.", default=None, required=False, dest="b38fasta", type=str)  ##REQUIRED for normalization, but normalization not requried to run script!
	g.add("--run-XML-pipeline", default=False, action='store_true',  help="Run original XML parse and variant_summary merge pipeline.  Not run by default", dest="run_xml")
	g.add("--run-TSV-pipeline", default=True, action='store_false',  help="Run new pipeline that merges just the ClinVar tsv files: variant_summary and submission_summary.", dest="run_tsv")
	g.add("-X", "--clinvar-xml", help="The local filename of the ClinVarFullRelase.xml.gz file. If not set, will download the latest from ClinVar NCBI FTP site.", dest="xml_file", default=dir_path+os.sep+"output_tmp"+os.sep+"ClinVarFullRelease_00-latest.xml.gz", type=str)
	g.add("-V", "--clinvar-variant-summary-table", help="The local filename of the variant_summary.txt.gz file. If not set, will download the latest from ClinVar NCBI FTP site.", dest="V_tsv_file", default=dir_path+os.sep+"output_tmp"+os.sep+"variant_summary.tsv.gz", type=str)
	g.add("-S", "--clinvar-submission-summary-table", help="The local filename of the submission_summary.txt.gz file. If not set, will download the latest from ClinVar NCBI FTP site.", dest="S_tsv_file", default=dir_path+os.sep+"output_tmp"+os.sep+"submission_summary.txt.gz", type=str)
	g.add("--output-prefix", help="Final output files will have this prefix", type=str, dest='output_prefix',default="")
	g.add("--output-dir", default=dir_path+os.sep+"output"+os.sep, help="Final output files will be located here", type=str, dest='output_dir')
	g.add("-N", "--new", help='Download all New.  Causes ouput_tmp to be removed and recreated and latest ClinVar Files to be downloaded', action='store_true', default=False, dest="download_new")
	g.add("--tmp-dir", default=dir_path+os.sep+"output_tmp"+os.sep, help="Temporary output direcotry for temp files ", dest="output_tmp", type=str)
	g.add("--rm-temp", default=True, action='store_false',  help="Causes temporary directories and temp files to NOT be removed when finished.  Default: Automaticly removed. ", dest="rm_tmp")
	g.add("--rm-output", default=False, action='store_true', help="Causes ouput directories to be deleted and recreated.  Requires -N option to also be set.  Default: False.", dest="rm_output")

	return p.parse_args()

def main():

	cli_args=parseArguments()
	print("Arguments Used for run:")
	print(cli_args)
	print()
	if cli_args.b37fasta is None and cli_args.b38fasta is None:
		sys.exit("ERROR: At least one fasta file must be provided.")
	else:
		if cli_args.b37fasta is not None and not checkExists(cli_args.b37fasta):
			sys.exit("Genome reference: file not found:\t"+cli_args.b37fasta)
		if cli_args.b38fasta is not None and not checkExists(cli_args.b38fasta):
			sys.exit("Genome reference: file not found:\t"+cli_args.b38fasta)

	#creates temp and outupt directories
	createDirectories(cli_args)

	#downloads clinvar and normalize_py
	downloadClinVarFiles(cli_args)

	#Runs the original, XML pipeline.  NOT default.
	if cli_args.run_xml:

		##TODO!!!!!!  Currntly only does one genome build at a time.  Should just ouput both at same time cli_args.output_tmp so xml doesn't have to be parsed twice
		##NOTE Multi done after parse and sort and normalize...  Grouping by allele doens't make sense when everything in is haplotyped based.  Need a new grouping stratagy for multi.
		##QUESTION: How many Multi show up only in multi and havie higher/lower stars?
		##NOTE!!! There are missing variants in the multi.
		##NOTE!!!  Just gonna build a new pipeline to just use tsv files. The Submission summary table can be joined to the Variant summary table on allele ID
		if cli_args.b37fasta is not None:
			runXMLpipeLine(cli_args, 'GRCh37',cli_args.b37fasta)
		if cli_args.b38fasta is not None:
			runXMLpipeLine(cli_args, 'GRCh38',cli_args.b38fasta)

	if cli_args.run_tsv:
		#TODO!!!  Do all the things!
		#QUESTION what do the haplotypes look like in pure tsv?
		pass

	#remves temp files.
	if cli_args.rm_tmp:
		shutil.rmtree(cli_args.output_tmp)

	return


if __name__ == '__main__':
	main()
	exit()
