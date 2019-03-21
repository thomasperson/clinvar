import argparse
import collections
import gzip
import os
import re
import pandas as pd
import sys
import time

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
		'VARIATION_ID': "The identifier assigned by ClinVar and used to build the URL, namely https://ncbi.nlm.nih.gov/clinvar/VariationID",
		'CLINICAL_SIGNIFICANCE': "Interpretation of the variation-condition relationship",
		'DATE_LAST_EVALUATED' : "The last date the variation-condition relationship was evaluated by this submitter",
		'DESCRIPTION' : "An optional free text description of the basis of the interpretation",
		'SUBMITTED_PHENOTYPE_INFO' :  "The name(s) or identifier(s)  submitted for the condition that was interpreted relative to the variant",
		'REPORTED_PHENOTYPE_INFO' : "The MedGen identifier/name combinations ClinVar uses to report the condition that was interpreted. 'na' means there is no public identifer in MedGen for the condition.",
		'REVIEW_STATUS' : "The level of review for this submission, namely http//www.ncbi.nlm.nih.gov/clinvar/docs/variation_report/#review_status",
		'COLLECTION_METHOD' : "The method by which the submitter obtained the information provided",
		'ORIGIN_COUNTS' : "The reported origin and the number of observations for each origin",
		'SUBMITTER' : "The submitter of this record",
		'SCV' :  "The accession and current version assigned by ClinVar to the submitted interpretation of the variation-condition relationship",
		'SUBMITTED_GENE_SYMBOL': "The gene symbol reported in this record",
		'EXPLANATION_OF_INTERPRETATION' : "More details if ClinicalSignificance is 'other' or 'drug response'",
		'GENE_SYMBOL' : "The offical symbol from HGNC for the gene in which this allele is found.  Will be null if the variant is intergenic or spans multiple genes.",
		'GENE_ID': "NCBI's identfier for the gene in which this allele is found.  Will be null if the variant is intergenic or spans multiple genes.",
		'ALLELE_ID' : "The identifier assigned by ClinVar to the simple allele. Will be null if the allele is complex.",
		'TYPE' : "The type of HGVS expression",
		'CLIN_PATH' : "1 for True, 0 for False.  If any Clinical Lab has returned the variant as a P/LP then set to 1.",



	}
	for key in HEADER:
		outVCF.write("##INFO=<ID="+key.upper()+",Number=1,Type=String,Description="+descriptions.get(key, key.upper())+">\n")
	with open(input_reference_genome_fai) as in_fai:
		for line in in_fai:
			chrom, length, _ = line.split("\t", 2)
			outVCF.write("##contig=<ID="+chrom.replace("chr", "")+",length="+length+">\n")
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

def table_to_vcf_oneline(input_table_path, input_reference_genome, output_vcf):
	""""This method creates a vcf with a subset of available tsv fields, merged
	into one vcf info field."""

	DayMonthYear =str(time.strftime("%d")+time.strftime("%b")+time.strftime("%Y")).upper()

	# validate args
	input_reference_genome_fai = input_reference_genome + ".fai"
	if not os.path.isfile(input_table_path):
		sys.exit("ERROR: %s not found" % input_table_path)
	if not os.path.isfile(input_reference_genome_fai):
		sys.exit("ERROR: %s (reference FASTA .fai) not found" %
				 input_reference_genome_fai)

	outVCF=open(output_vcf,'w')

	# read input table. low_memory allows dtypes to be inferred
	df = pd.read_table(gzopen(input_table_path), low_memory=False)

	missing_columns = {"CHROM", "POS", "REF", "ALT"} - set(df.columns)
	if missing_columns:
		sys.exit("ERROR: %s is missing columns: %s" % (input_table_path, str(missing_columns)))

	outVCF.write("##fileformat=VCFv4.2\n##source=clinvar\n")

	df = df[['VARIATION_ID','CHROM','POS','REF','ALT','ASSEMBLY','NAME','TYPE','GENE_SYMBOL','CLINICAL_SIGNIFICANCE','CLIN_SIG_SIMPLE','LAST_EVALUATED','RS_NUM_DBSNP','REVIEW_STATUS','NUMBER_OF_SUBMITTERS','GOLD_STARS','CONFLICTED','CLIN_PATH']]

	HEADER=list(df.columns.values)

	loc_column = ['CHROM', 'POS', 'REF', 'ALT']

	vcf_info_key="CLINVAR_"+DayMonthYear

	outVCF.write("##INFO=<ID=CLINVAR_"+DayMonthYear+",Number=.,Type=String,Description=\"ClinVar Annotations.  Format: '"+"|".join(HEADER)+"'\">\n")

	with open(input_reference_genome_fai) as in_fai:
		for line in in_fai:
			chrom, length, _ = line.split("\t", 2)
			outVCF.write("##contig=<ID="+chrom.replace("chr", "")+",length="+length+">\n")
	outVCF.write("##reference="+input_reference_genome+"\n")

	outVCF.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])+"\n")
	for i, table_row in df.iterrows():
		vcf_row = []
		vcf_row.append(table_row["CHROM"])
		vcf_row.append(table_row["POS"])
		vcf_row.append('.')  # ID
		vcf_row.append(table_row["REF"])
		vcf_row.append(table_row["ALT"])
		vcf_row.append('.')  # QUAL
		vcf_row.append('.')  # FILTER

		info_field = collections.OrderedDict()

		ClinVarOutInfoField=vcf_info_key+"="

		for key in HEADER:
			if pd.isnull(table_row[key]):
				continue
			value = str(table_row[key]).strip()
			if key =="NAME":
				value = value.replace(":", "").replace(";", "%3B").replace("=", "%3D").replace(" ", "").replace(",", "").strip()
			else:
				value = value.replace(":", "%3A").replace(";", "%3B").replace("=", "%3D").replace(" ", "_").replace(",", "").strip().upper()
			ClinVarOutInfoField+=value+"|"

		vcf_row.append(ClinVarOutInfoField.strip("|"))

		outVCF.write("\t".join(map(str, vcf_row))+"\n")

	outVCF.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',dest='input_table_path', type=str,help="Tab-delimited input table")
	parser.add_argument('--fasta',dest='input_reference_genome', type=str, help="Reference FASTA used. NOTE: The associated .fai file, e.g. b38.fa.fai, is necessary for the VCF header generation")
	parser.add_argument('-o', '--outfile', dest='outfile', type=str, help="Output vcf file.")
	args = parser.parse_args()

	table_to_vcf(args.input_table_path, args.input_reference_genome, arg.outfile)
