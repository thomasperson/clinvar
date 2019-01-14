#!/usr/bin/env python

import argparse
import gzip
import sys

import parse_clinvar_xml as pcx
# recommended usage:   RUN ON A SORTED FILE!
# ./group_by_allele.py < clinvar_combined.tsv > clinvar_alleles.tsv


def group_by_allele(infile_path, outfile_path):
    """Run through a SORTED clinvar_table.tsv file from the parse_clinvar_xml script, and make it unique on CHROM POS REF ALT

    Args:
        infile: Input file stream for reading clinvar_table_sorted.tsv
        outfile: Output file stream to write to.
    """
    infile=open(infile_path,'r')
    outfile=open(outfile_path,'w')
    #header = next(infile)
    #outfile.write("\t".join(pcx.HEADER)+"\n")
    #column_names = header.strip('\n').split('\t')

    last_data = None
    last_unique_id = None
    counter = 0

    for line in infile:
        if line.strip()=="":
            continue
        data = dict(zip(pcx.HEADER, [x.strip() for x in line.strip().split('\t')]))
        unique_id = '_'.join([data['CHROM'], str(data['POS']), data['REF'], data['ALT']])
        if unique_id == last_unique_id:
            data = group_alleles(last_data, data)
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

def group_alleles(data1, data2):
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

    group_by_allele(args.infile, args.outfile)


if __name__ == '__main__':
    main()
