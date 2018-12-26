"""
Alternate implementation of master.py.  Removing pysam, pypez, non clinvar related functions.  .

Run with -h to see all options.
"""

import configargparse
from datetime import datetime
import ftplib
import os
import sys
from distutils import spawn

try:
    import configargparse
    import pandas   # make sure all dependencies are installed
except ImportError as e:
    sys.exit("ERROR: Python module not installed. %s. Please run 'pip install -r requirements.txt' " % e)
for executable in ['wget', 'tabix', 'vt']:
    assert spawn.find_executable(executable), "Command %s not found, see README" % executable


def parseArguments():
    p = configargparse.getArgParser()
    g = p.add_argument_group('main args')
    g.add("--b37-genome", help="b37 .fa genome reference file", default=None, required=False)
    g.add("--b38-genome", help="b38 .fa genome reference file. NOTE: chromosome names must be like '1', '2'.. 'X', 'Y', 'MT'.", default=None, required=False)
    g.add("-X", "--clinvar-xml", help="The local filename of the ClinVarFullRelase.xml.gz file. If not set, grab the latest from NCBI.")
    g.add("-S", "--clinvar-variant-summary-table", help="The local filename of the variant_summary.txt.gz file. If not set, grab the latest from NCBI.")
    #g.add("-E", "--exac-sites-vcf",  help="ExAC sites vcf file. If specified, a clinvar table with extra ExAC fields will also be created.")
    #g.add("-GE", "--gnomad-exome-sites-vcf",  help="gnomAD exome sites vcf file. If specified, a clinvar table with extra gnomAD exome info fields will also be created.")
    #g.add("-GG", "--gnomad-genome-sites-vcf",  help="gnomAD genome sites vcf file. If specified, a clinvar table with extra gnomAD genome info fields will also be created.")
    g.add("--output-prefix", default="../output/", help="Final output files will have this prefix")
    g.add("--tmp-dir", default="./output_tmp", help="Temporary output files will have this prefix")
    g = p.add_mutually_exclusive_group()
    g.add("--single-only", dest="single_or_multi", action="store_const", const="single", help="Only generate the single-variant tables")
    g.add("--multi-only", dest="single_or_multi", action="store_const", const="multi", help="Only generate the multi-variant tables")

    return p.parse_args()


def main():
    

if __name__ == '__main__':
    main()
