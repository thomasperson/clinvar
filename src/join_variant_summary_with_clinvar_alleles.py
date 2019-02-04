#!/usr/bin/env python
import sys
import pandas as pd
import argparse
import re

def join_variant_summary_with_clinvar_alleles(variant_summary_table, clinvar_alleles_table,genome_build_id,out_name):

	variant_summary = pd.read_csv(variant_summary_table, sep="\t", index_col=False, compression="gzip",low_memory=False)
	print ("variant_summary raw", variant_summary.shape)

	clinvar_alleles = pd.read_csv(clinvar_alleles_table, sep="\t",index_col=False, low_memory=False)

	print ("clinvar_alleles raw", clinvar_alleles.shape)

	FINAL_HEADER = list(clinvar_alleles.columns.values) + ['GOLD_STARS', 'CONFLICTED']

	# use lowercase names and replace . with _ in column names:
	variant_summary = variant_summary.rename(columns=dict((col, col.upper().strip().strip("#").replace(".", "_"))for col in variant_summary.columns))
	# rename first column to allele_id:
	variant_summary = variant_summary.rename(columns={'ALLELEID': 'ALLELE_ID', 'CLINICALSIGNIFICANCE': 'CLINICAL_SIGNIFICANCE','REVIEWSTATUS': 'REVIEW_STATUS','LASTEVALUATED':'LAST_EVALUATED'})
	#variant_summary = variant_summary.rename(columns={variant_summary.columns[0]: "ALLELE_ID"})

	# extract relevant columns for the correct assembly and
	# rename clinicalsignificance, reviewstatus, lastevaluated:
	variant_summary = variant_summary[variant_summary['ASSEMBLY'] == genome_build_id]
	variant_summary = variant_summary[['ALLELE_ID' ,'CLINICAL_SIGNIFICANCE','REVIEW_STATUS','LAST_EVALUATED']]


	# remove the duplicated records in variant summary due to alternative loci such as PAR but would be problematic for rare cases like translocation
	variant_summary=variant_summary.drop_duplicates()
	print ("variant_summary after filter", variant_summary.shape)

	# remove clinical_significance and review_status from clinvar_alleles:
	clinvar_alleles = clinvar_alleles.drop(['CLINICAL_SIGNIFICANCE', 'REVIEW_STATUS','LAST_EVALUATED'], axis=1)

	# pandas is sensitive to some rows having allele_id joined on ;, causing
	# an object dtype, with some entries being ints and others strs
	clinvar_alleles['ALLELE_ID'] = clinvar_alleles['ALLELE_ID'].astype(str)
	variant_summary['ALLELE_ID'] = variant_summary['ALLELE_ID'].astype(str)

	print ("clinvar_alleles after filter", clinvar_alleles.shape)

	# join the two tables on allele_id:
	df = pd.merge(clinvar_alleles, variant_summary,on='ALLELE_ID', how='inner')
	print ("merged raw", df.shape)

	# map review_status to gold starts:
	gold_star_map = {
		'no assertion provided': 0,
		'no assertion for the individual variant': 0,
		'no assertion criteria provided': 0,
		'criteria provided, single submitter': 1,
		'criteria provided, conflicting interpretations': 1,
		'criteria provided, multiple submitters, no conflicts': 2,
		'reviewed by expert panel': 3,
		'practice guideline': 4,
		'-':'-'
	}
	df['GOLD_STARS'] = df.REVIEW_STATUS.map(gold_star_map)

	# The use of expressions on clinical significance on ClinVar aggregate records (RCV) https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/#conflicts
	# conflicted = 1 if using "conflicting"
	df['CONFLICTED'] = df['CLINICAL_SIGNIFICANCE'].str.contains("onflicting", case=False).astype(int)

	# reorder columns just in case
	df = df.ix[:, FINAL_HEADER]
	print ("merged final", df.shape)

	df= df.sort_values(['CHROM', 'POS', 'REF', 'ALT'])

	df.to_csv(out_name, sep="\t", index=False, compression='gzip')

	return

def join_variant_summary_with_submission_summary(variant_summary_file, grouped_submission_summary, grouped_var_citations, genome_build_id, out_name):

	raw_grouped_var_citations = pd.read_csv(grouped_var_citations, sep="\t", index_col=False, low_memory=True)
	print ("raw_grouped_var_citations", raw_grouped_var_citations.shape)
	raw_grouped_submission_summary = pd.read_csv(grouped_submission_summary, sep="\t", index_col=False, low_memory=True)
	print ("raw_grouped_submission_summary", raw_grouped_submission_summary.shape)

	raw_grouped_submission_summary=raw_grouped_submission_summary.drop(['SUBMITTED_GENE_SYMBOL'], axis=1)

	raw_grouped_submission_summary['EXPLANATION_OF_INTERPRETATION'] = raw_grouped_submission_summary['EXPLANATION_OF_INTERPRETATION'].str.replace(r'(\-;)+','').str.replace(r'(;\-)+','').str.replace(r'^\-','').str.replace(r'\-$','')
	raw_grouped_submission_summary['DESCRIPTION'] = raw_grouped_submission_summary['DESCRIPTION'].str.replace(r'(\-;)+','').str.replace(r'(;\-)+','').str.replace(r'^\-','').str.replace(r'\-$','')

	raw_grouped_var_citations['VARIATION_ID'] = raw_grouped_var_citations['VARIATION_ID'].astype(int)
	raw_grouped_submission_summary['VARIATION_ID'] = raw_grouped_submission_summary['VARIATION_ID'].astype(int)

	# join the two tables on allele_id:
	df = pd.merge(raw_grouped_submission_summary, raw_grouped_var_citations,on='VARIATION_ID', how='left')

	del raw_grouped_submission_summary
	del raw_grouped_var_citations

	variant_summary = pd.read_csv(variant_summary_file, sep="\t", index_col=False, compression="gzip",low_memory=False)
	print ("variant_summary raw", variant_summary.shape)
	# use lowercase names and replace . with _ in column names:
	#variant_summary = variant_summary.rename(columns=dict((col, col.upper().replace("#", "").replace(".", "_").strip().strip("#").replace(" ", "_").replace("/","_").replace("(","").replace(")","")) for col in variant_summary.columns))
	#rename first column to allele_id:
	variant_summary = variant_summary.rename(columns={"#AlleleID": 'ALLELE_ID',
											"Type": 'TYPE',
											"Name": 'NAME',
											"GeneID": 'GENE_ID',
											"GeneSymbol": 'GENE_SYMBOL',
											"ClinicalSignificance": 'CLINICAL_SIGNIFICANCE',
											"ClinSigSimple" : "CLIN_SIG_SIMPLE",
											"LastEvaluated": 'LAST_EVALUATED',
											"RS# (dbSNP)" : "RS_NUM_DBSNP",
											"nsv/esv (dbVar)": "DBVAR_ID",
											"RCVaccession" : "RCV_ACCESSION",
											"PhenotypeIDS" : "PHENOTYPE_IDS",
											"PhenotypeList" : "PHENOTYPE_LIST",
											"Origin" : "ORIGIN",
											"OriginSimple" : "ORIGIN_SIMPLE",
											"Assembly" : "ASSEMBLY",
											"ChromosomeAccession" : "CHROMOSOME_ACCESSION",
											"Chromosome" : "CHROM",
											"Start" : "START",
											"Stop" : "STOP",
											"ReferenceAllele": "REF",
											"AlternateAllele": "ALT",
											"Cytogenetic": "CYTOGENETIC",
											"ReviewStatus" :'REVIEW_STATUS',
											"NumberSubmitters": "NUMBER_OF_SUBMITTERS",
											"Guidelines": "GUIDLINES",
											"TestedInGTR" : "TESTED_IN_GTR",
											"OtherIDs": "OTHER_IDS",
											"SubmitterCategories": "SUBMITTER_CATEGORIES",
											"VariationID": "VARIATION_ID"
																		})

	variant_summary = variant_summary[variant_summary['ASSEMBLY'] == genome_build_id]
	print ("variant_summary buildID", variant_summary.shape)
	variant_summary=variant_summary.drop_duplicates()
	print ("variant_summary drop_duplicates", variant_summary.shape)

	variant_summary['VARIATION_ID'] = variant_summary['VARIATION_ID'].astype(int)

	# join the two tables on VARIATION_ID:
	df = pd.merge(df, variant_summary,on='VARIATION_ID', how='inner')

	del variant_summary

	print ("merged raw", df.shape)

	# map review_status to gold starts:
	gold_star_map = {
		'no assertion provided': 0,
		'no assertion for the individual variant': 0,
		'no assertion criteria provided': 0,
		'criteria provided, single submitter': 1,
		'criteria provided, conflicting interpretations': 1,
		'criteria provided, multiple submitters, no conflicts': 2,
		'reviewed by expert panel': 3,
		'practice guideline': 4,
		'-':'-'
	}
	df['GOLD_STARS'] = df.REVIEW_STATUS.map(gold_star_map)
	df['POS']=df.START

	# The use of expressions on clinical significance on ClinVar aggregate records (RCV) https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/#conflicts
	# conflicted = 1 if using "conflicting"
	df['CONFLICTED'] = df['CLINICAL_SIGNIFICANCE'].str.contains("onflicting", case=False).astype(int)

	df= df.sort_values(['CHROM', 'POS', 'REF', 'ALT'])

	df.to_csv(out_name, sep="\t", index=False, compression='gzip')
	return

def clinStar():
	pass



def main():
	parser = argparse.ArgumentParser(description='Joins the grouped ClinVar xml parse output to the ClinVar TSV file')
	parser.add_argument("-S", "--clinvar-variant-summary-table", help="The local variant_summary.txt.gz file.", dest="variant_summary_table", default=dir_path+os.sep+"output_tmp"+os.sep+"ClinVar.tsv.gz", type=str)
	parser.add_argument('-i', '--infile', type=str, help="Sorted and grouped output of group_by_allele.py.", dest="clinvar_alleles_table")
	parser.add_argument('-o', '--outfile', type=str, help="Output file to store grouped alleles.", dest="out_name")
	parser.add_argument('-G', '--build', choices=['GRCh37', 'GRCh38'],dest='genome_build_id')
	args = parser.parse_args()

	assert args.out_name.endswith('.gz'), ("Provide a filename with .gz extension as the output will be gzipped")

	join_variant_summary_with_clinvar_alleles(args.variant_summary_table, args.clinvar_alleles_table, args.genome_build_id,args.out_name)

	return




if __name__ == "__main__":
	main()
	exit()
