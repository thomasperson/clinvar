## This code is currently being refactored so is currently in beta.  Goal is to make it as pure python as possible with support for Python2.7.x+ and Python3.6.x+ on Windows and Linux, though some functionality may be missing if C complied programs are not present.  
## NOTE: This Fork is not supported or authorized by the MacAurthur Lab in anyway, though I do thank them for their work building the initial framework. Also some functionality found in the original repo has being removed as has been deemed not needed.   

#### In 1 sentence

This repo provides a tool to create a VCF containing individual lab variant annotations from ClinVar.  The XML pipeline provides the original MacAurthur pipeline with some bug fixes. The new TSV only pipeline merges various TSV files provided by ClinVar into one file and also provides a new flag to indicate if any Clinical Lab submitted the variant as Pathogenic or Likely Pathogenic. This is to help with triaging of conflicted varaiants.  

#### Output Files

The full set of tables generated by the pipeline for the latest clinvar release is stored in the [output](output/) folder.
These tables are separated by genome build ([output/b37](output/GRCh37) or [output/GRCh38](output/GRCh38)) and whether they represent
simple mono-allelic variants (`../single/` sub-directory),
or more complex variants such as compound-heterozygotes, haplotypes, etc. (`../multi/` sub-directory) - for a total of 4 directories that contain tables.
The TSV pipeline does not processes compound-heterozygous variants.  Only directories for build selected for are created.  

The directories contains the following files:
* __clinvar_allele_trait_pairs.*.tsv.gz__: table where each row represents an allele-trait pair.
* __clinvar_alleles.*.tsv.gz__: table where each row represents a single variant allele. This is generated by grouping _clinvar_allele_trait_pairs.tsv.gz_ by allele.
* __clinvar_alleles.*.vcf.gz__: _clinvar_alleles.tsv.gz_ converted to VCF format.
* __clinvar_alleles_stats.*.txt__:  summary of the different columns in _clinvar_alleles.*.tsv.gz_, along with some basic stats on the different values that appear in each column.


#### Motivation

[ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) is a public database hosted by NCBI for the purpose of collecting assertions as to genotype-phenotype pairings in the human genome. One common use case for ClinVar is as a catalog of genetic variants that have been reported to cause Mendelian disease. By including lab specific variant classification, more nuanced determinations can be made about variants seen in the wild.

ClinVar makes its data available via [FTP](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/) in three formats: XML, TXT, and VCF. We found that none of these files were ideally suited for our purposes. The XML file is large and complex, with multiple entries for the same genomic variant, making it difficult to quickly look up a variant of interest. In addition, both the XML and TXT representations are not guaranteed to be unique on genomic coordinates, and also contain many genomic coordinates that have been parsed from HGVS notation, and therefore may be right-aligned (in contrast to left alignment, the standard for VCF) and may also be non-minimal (containing additional nucleotides of context to the left or right of a given variant).

#### Processing XML Pipeline

To create a flat representation of ClinVar XML suited for our purposes, we took several steps, encapsulated in the pipeline [./run_clinvar_parse.py](./run_clinvar_parse.py)

1. Download the latest XML and TXT dumps from ClinVar FTP.
2. Parse the XML file using [src/parse_clinvar_xml.py](src/parse_clinvar_xml.py) to extract fields of interest into a flat file.
3. Normalize using [our Python implementation](src/normalize.py) of [vt normalize](http://genome.sph.umich.edu/wiki/Variant_Normalization) (see [[Tan 2015]]).
4. Group the allele-trait records by allele using [src/group_by_allele.py](src/group_by_allele.py) to aggregate interpretations from multiple submitters by allele, independent of conditions.
5. Join the TXT file using [src/join_variant_summary_with_clinvar_alleles.py](src/join_variant_summary_with_clinvar_alleles.py) to aggregate interpretations from multiple submitters independent of conditions.
6. Generate the VCF file and other tables based on the file created in 5.


#### Processing TSV Pipeline

To create a merged TSV of ClinVar suited for our purposes, we took several steps, encapsulated in the pipeline [./run_clinvar_parse.py](./run_clinvar_parse.py)

1. Download the latest TSV dumps from ClinVar FTP.
2. Group the submission summary records by VaraintID using [src/group_by_allele.py](src/group_by_allele.py) to aggregate interpretations from multiple submitters by allele, independent of conditions.
3. Group the variant summary records by VaraintID using [src/group_by_allele.py](src/group_by_allele.py) to aggregate interpretations from multiple submitters by allele, independent of conditions.
4. Join the TSV files using [src/join_variant_summary_with_clinvar_alleles.py](src/join_variant_summary_with_clinvar_alleles.py) to aggregate interpretations from multiple submitters independent of conditions.
5. Add in the CLIN_PATH flag to indicate which of the Variants had a Clinical Lab submit the variant with a Pathogenic or Likely Pathogenic Clinical Significance
6. Normalize using [our Python implementation](src/normalize.py) of [vt normalize](http://genome.sph.umich.edu/wiki/Variant_Normalization) (see [[Tan 2015]]).
7. Generate the VCF file and other tables based on the final TSV created in 6.

&dagger;Because a ClinVar record may contain multiple assertions of Clinical Significance, we defined the following additional columns to represent the clinical significances(https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig):

+ `pathogenic`,`likely_pathogenic`,`uncertain_significance`,`likely_benign` and `benign` are the counts of submissions reported the variants as "Pathogenic","Likely pathogenic","Uncertain significance","Likely benign" and "Benign" respectively.
+ `conflicted` is `1` if the variant is reported as "Conflicting interpretation of pathogenicity" or "conflicting data from submitters".

#### Usage

The pipeline scripts expect the following programs to be available on your system (and in your `$PATH`):

[python](https://www.python.org/), [tabix](http://www.htslib.org/download/), [bgzip](http://www.htslib.org/download/ ),[vcf-sort](https://vcftools.github.io/index.html)

If not all are present some steps may be skipped.

To run the default TSV pipeline:

```
pip install --user --upgrade -r requirements.txt
python run_clinvar_parse.py --b38-genome /path/to/b38.fa
```

See `python run_clinvar_parse.py -h` for additional options.  The TSV pipeline runs by default.  To run *ONLY* the XML pipeline, add the flags `--run-TSV-pipeline` AND `--run-XML-pipeline`.


Additional helper scripts are available for users to use check the processing results:
[src/grab_interesting_variations.py](src/grab_interesting_variations.py) to extract the raw xml entry given a list of ClinVar variation IDs.
`python grab_interesting_variations.py <ClinVarFullRelease.xml.gz> <comma-separated list of variation IDs> <out.xml.gz> `
[src/diff_clinvar_alleles.py](src/diff_clinvar_alleles.py) to compare the differences of two ClinVar_alleles_*.tsv.gz output files.
`python diff_clinvar_alleles.py <clinvar_alleles.A.tsv.gz> <clinvar_alleles.B.tsv.gz>`

#### Usage notes

Because ClinVar XML contains a great deal of data complexity, a deliberate decision was made to *not* attempt to capture all fields in the resulting file. An effort to capture a subset of fields that were believed would be most useful for genome-wide filtering, and also included `variation_id` as a column to enable the user to look up additional details on the ClinVar website. For instance, the page for the variant with `variation_id` 7105 is located at [ncbi.nlm.nih.gov/clinvar/variation/7105/](http://www.ncbi.nlm.nih.gov/clinvar/variation/7105/).

#### Limitation

This tool has some limited functionality on Windows due to the difficulty in installing some tools/packages. For example, the `normalize.py` step requires the pysam package be available.  Additional tools such as [tabix](http://www.htslib.org/download/) and [bgzip](http://www.htslib.org/download/) are also required for compressing and indexing of vcf files.  These steps will also be skipped if these tools/packages are not installed.  

#### License, terms, and conditions

ClinVar data, as a work of the United States federal government, are in the public domain and are redistributed here under [the same terms](http://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/) as they are distributed by ClinVar itself. Importantly, note that ClinVar data are "not intended for direct diagnostic use or medical decision-making without review by a genetics professional". The code in this repository is distributed under an MIT license.

[Tan 2015]: http://www.ncbi.nlm.nih.gov/pubmed/25701572 "Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants. Bioinformatics. 2015 Jul 1;31(13):2202-4. doi: 10.1093/bioinformatics/btv112. Epub 2015 Feb 19. PubMed PMID: 25701572."

#### How to cite the XML pipeline
Zhang X, Minikel EV, O'Donnell-Luria AH et al. [ClinVar data parsing](https://wellcomeopenresearch.org/articles/2-33/v1) [version 1; referees: 2 approved]. Wellcome Open Res 2017, 2:33
(doi: 10.12688/wellcomeopenres.11640.1)

#### The TSV pipeline has not been published
Please message the repo with any feature requests, questions, bug complaints, or links to papers that might utilize it!
