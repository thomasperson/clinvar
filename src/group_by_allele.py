#!/usr/bin/env python

import argparse
import gzip
import sys
import io
from collections import OrderedDict


import parse_clinvar_xml as pcx
# recommended usage:   RUN ON A SORTED FILE!
# ./group_by_allele.py < clinvar_combined.tsv > clinvar_alleles.tsv

def group_submission_summary_file(infile_path, outfile_path):
	""""""
	inFile = None
	try:
		inFile = gzip.open(infile_path, 'rt')
	except ValueError:
		# Workaround for Python 2.7 under Windows
		inFile = gzip.open(infile_path, "r")
	submission_file=dict()
	header=[]
	header_replace={"VariationID": 'VARIATION_ID',
				"ClinicalSignificance": 'CLINICAL_SIGNIFICANCE_INDV_SUB',
				"DateLastEvaluated":'SUBMISSION_DATES',
				"Description":'DESCRIPTION',
				"SubmittedPhenotypeInfo":'SUBMITTED_PHENOTYPE_INFO',
				"ReportedPhenotypeInfo":'REPORTED_PHENOTYPE_INFO',
				"ReviewStatus":'REVIEW_STATUS_INDV_SUB',
				"CollectionMethod":'COLLECTION_METHOD',
				"OriginCounts":'ORIGIN_COUNTS',
				"Submitter":'SUBMITTER',
				"SubmittedGeneSymbol":'SUBMITTED_GENE_SYMBOL',
				"ExplanationOfInterpretation":'EXPLANATION_OF_INTERPRETATION'
		}

	ClinSigField = OrderedDict([ ('BENIGN' , 0) ,
								('LIKELY_BENIGN' , 0),
								('UNCERTAIN_SIGNIFICANCE',0),
								('LIKELY_PATHOGENIC',0),
								('PATHOGENIC',0),
								('DRUG_RESPONSE',0),
								('ASSOCIATION',0),
								('RISK_FACTOR',0),
								('PROTECTIVE',0),
								('AFFECTS',0),
								('CONFLICTING_DATA_FROM_SUBMITTERS',0),
								('OTHER',0),
								('NOT_PROVIDED',0)   ])

	for lineIndex,line in enumerate(inFile):
		if line.strip()=="":
			continue
		fields = [x.strip() for x in line.strip("#").split('\t')]
		if line.startswith("#VariationID") and lineIndex>4:
			for f in fields:
				if f in header_replace:
					header.append(header_replace[f])
				else:
					header.append(f)
		elif not line.startswith("#"):
			if fields[0] in submission_file:
				for i,f in enumerate(fields[1:]):
					submission_file[fields[0]][i]=submission_file[fields[0]][i]+";"+f
				if "/" in fields[1]:
					newField1=fields[1].split("/")[-1].upper().strip().replace(" ", "_")
					submission_file[fields[0]][-1][newField1]+=1
				else:
					submission_file[fields[0]][-1][fields[1].upper().strip().replace(" ", "_")]+=1
			else:
				value=[]
				if "/" in fields[1]:
					newField1=fields[1].split("/")[-1].upper().strip().replace(" ", "_")
					fields[1]=newField1
				for f in fields[1:]:
					value.append(f)
				V_ClinSigField=ClinSigField.copy()
				V_ClinSigField[fields[1].upper().strip().replace(" ", "_")]=1
				value.append(V_ClinSigField)
				submission_file[fields[0]]=value
	inFile.close()
	outFile = open(outfile_path,'w')
	outFile.write("\t".join(header)+"\t")
	outFile.write("\t".join(ClinSigField.keys())+"\n")
	try:
		for k,v in submission_file.iteritems():
			outFile.write(k+'\t')
			outFile.write("\t".join(v[:-1])+"\t")
			outFile.write("\t".join(str(e) for e in v[-1].values())+"\n")
	except AttributeError:
		for k,v in submission_file.items():
			outFile.write(k+'\t')
			outFile.write("\t".join(v[:-1])+"\t")
			outFile.write("\t".join(str(e) for e in v[-1].values())+"\n")
	outFile.close()

	return submission_file

def group_var_citations(infile_path, outfile_path):
	""""""
	inFile = open(infile_path, "r")
	header=['VARIATION_ID','CITATION_SOURCE:CITATION_ID']
	var_cite_file=dict()
	for line in inFile:
		if line.strip()=="":
			continue
		fields = [x.strip() for x in line.strip("#").split('\t')]
		if line.startswith("#"):
			continue
		if fields[1] in var_cite_file:
			var_cite_file[fields[1]].append(fields[-2]+":"+fields[-1])
		else:
			var_cite_file[fields[1]]=[fields[-2]+":"+fields[-1]]
	inFile.close()
	outFile = open(outfile_path,'w')
	outFile.write("\t".join(header)+"\n")
	try:
		for k,v in var_cite_file.iteritems():
			outFile.write(k+'\t')
			outFile.write(";".join(v)+"\n")
	except AttributeError:
		for k,v in var_cite_file.items():
			outFile.write(k+'\t')
			outFile.write(";".join(v)+"\n")
	outFile.close()

	return var_cite_file


def group_by_allele_rawXML(infile_path, outfile_path):
	"""Run through a SORTED clinvar_table.tsv file from the parse_clinvar_xml script, and make it unique on CHROM POS REF ALT

	Args:
		infile: Input file stream for reading clinvar_table_sorted.tsv
		outfile: Output file stream to write to.
	"""
	infile=open(infile_path,'r')
	outfile=open(outfile_path,'w')

	last_data = None
	last_unique_id = None
	counter = 0

	for line in infile:
		if line.strip()=="":
			continue
		data = dict(zip(pcx.HEADER, [x.strip() for x in line.strip().split('\t')]))
		unique_id = '_'.join([data['CHROM'], str(data['POS']), data['REF'], data['ALT']])
		if unique_id == last_unique_id:
			data = group_alleles_rawXML(last_data, data)
		elif last_data is not None:
			# note that using a comprehension instead of just data.values() preserves column order
			# the next line (data) is not duplicated as the current line(last_data) then just print last_data
			outfile.write('\t'.join([last_data[colname] for colname in pcx.HEADER])+'\n')
		last_data = data
		last_unique_id = unique_id
		counter += 1

	if last_data is not None:
		outfile.write('\t'.join([last_data[colname] for colname in pcx.HEADER])+'\n')
	else:
		raise ValueError("%s has 0 records" % infile)

	infile.close()
	outfile.close()

def group_alleles_rawXML(data1, data2):
	"""Group two variants with same genomic coordinates.

	Args:
		data1: dictionary of column-name, value pairs for table record #1
		data2: dictionary of column-name, value pairs for table record #2
	"""

	if (data1['CHROM'], data1['POS'], data1['REF'], data1['ALT']) != (data2['CHROM'], data2['POS'], data2['REF'], data2['ALT']):
		raise ValueError("data1 variant id != data2 variant id: %s != %s" % (data1, data2))

	combined_data = data1  # this sets defaults, now we fix it:

	# 'pathogenic', 'benign', 'conflicted', 'gold_stars',
	# concatenate columns that may have lists of values
	loc_column = ['CHROM','POS','REF','ALT']
	num_field = ['PATHOGENIC', 'LIKELY_PATHOGENIC','UNCERTAIN_SIGNIFICANCE','LIKELY_BENIGN', 'BENIGN']
	info_column = [x for x in pcx.HEADER if x not in loc_column and x not in num_field]

	for column_name in info_column:
		all_non_empty_values = filter(lambda s: s, data1[column_name].split(';') + data2[column_name].split(';'))
		# deduplicate values, while preserving order
		deduplicated_values = []
		for value in all_non_empty_values:
			if value not in deduplicated_values:
				deduplicated_values.append(value)

		combined_data[column_name] = ';'.join(deduplicated_values)

	for column_name in num_field:
		combined_data[column_name]=str(int(data1[column_name])+int(data2[column_name]))

	return combined_data

def main():
	parser = argparse.ArgumentParser(description='De-duplicate the output from parse_clinvar_xml.py')
	parser.add_argument('-i', '--infile', type=str, help="Sorted raw single output of parse_clinvar_xml.py.  i.e. clinvar_table_raw.single.GRCh38.sorted.tsv")
	parser.add_argument('-o', '--outfile', type=str, help="Output file to store grouped alleles.")
	args = parser.parse_args()

	group_by_allele_rawXML(args.infile, args.outfile)


if __name__ == '__main__':
	main()
