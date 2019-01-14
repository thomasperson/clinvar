import argparse
import collections
import gzip
import os
import re
import pandas as pd
import sys

#from parse_clinvar_xml import HEADER


def gzopen(path, mode='r', verbose=True):
	if path.endswith(".gz"):
		return gzip.open(path, mode)
	else:
		return open(path, mode)


def table_to_vcf(input_table_path, input_reference_genome, output_vcf):
	# validate args
	input_reference_genome_fai = input_reference_genome + ".fai"
	if not os.path.isfile(input_table_path):
		sys.exit("ERROR: %s not found" % input_table_path)
	if not os.path.isfile(input_reference_genome_fai):
		sys.exit("ERROR: %s (reference FASTA .fai) not found" %
				 input_reference_genome_fai)

	 outVCF=open(output_vcf,'w')

	# read input table. low_memory allows dtypes to be inferred
	t = pd.read_table(gzopen(input_table_path), low_memory=False)

	HEADER=list(t.columns.values)

	missing_columns = {"CHROM", "POS", "REF", "ALT"} - set(t.columns)
	if missing_columns:
		sys.exit("ERROR: %s is missing columns: %s" % (input_table_path, str(missing_columns)))

	outVCF.write("##fileformat=VCFv4.1\n##source=clinvar\n")

	descriptions = {    ##TODO!  UPDATE AND EXPAND!!!
		'GOLD_STARS': "Number of gold stars as shown on clinvar web pages to summarize review status. Lookup table described at http://www.ncbi.nlm.nih.gov/clinvar/docs/details/ was used to map the REVIEW_STATUS value to this number.",
	}
	for key in HEADER:
		outVCF.write("##INFO=<ID="+key.upper()"+,Number=1,Type=String,Description="+descriptions.get(key, key.upper())+">\n")
	with open(input_reference_genome_fai) as in_fai:
		for line in in_fai:
			chrom, length, _ = line.split("\t", 2)
			outVCF.write("##contig=<ID="+chrom.replace("chr", "")+",length="+length+">\n"
	outVCF.write("##reference="+input_reference_genome+"\n")

	outVCF.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])+"\n")
	for i, table_row in t.iterrows():
		vcf_row = []
		vcf_row.append(table_row["CHROM"])
		vcf_row.append(table_row["POS"])
		vcf_row.append('.')  # ID
		vcf_row.append(table_row["REF"])
		vcf_row.append(table_row["ALT"])
		vcf_row.append('.')  # QUAL
		vcf_row.append('.')  # FILTER

		info_field = collections.OrderedDict()

		# from VCF spec:
		#    INFO - additional information: (String, no white-space, semi-colons, or equals-signs permitted; commas are
		#    permitted only as delimiters for lists of values) INFO fields are encoded as a semicolon-separated series of short
		#    keys with optional values in the format: <key>=<data>[,data].
		loc_column = ['CHROM', 'POS', 'REF', 'ALT']
		for key in HEADER:
			if key not in loc_column:
				if pd.isnull(table_row[key]):
					continue
				value = str(table_row[key])
				value = re.sub('\s*[,]\s*', '..', value)  # replace , with ..
				value = re.sub('\s*[;]\s*', '|', value)  # replace ; with |
				value = value.replace("=", " eq ").replace(" ", "_")

				info_field[key.upper()] = value
		vcf_row.append(";".join([key+"="+value for key, value in info_field.items()]))

		outVCF.write("\t".join(map(str, vcf_row))+"\n")

	outVCF.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',dest='input_table_path', type=str,help="Tab-delimited input table")
	parser.add_argument('--fasta',dest='input_reference_genome', type=str, help="Reference FASTA used. NOTE: The associated .fai file, e.g. b38.fa.fai, is necessary for the VCF header generation")
	parser.add_argument('-o', '--outfile', dest='outfile', type=str, help="Output vcf file.")
	args = parser.parse_args()

	table_to_vcf(args.input_table_path, args.input_reference_genome, arg.outfile)
