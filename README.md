## scAI (single-cell allelic imbalance)
Scripts for allele-specific analyses of single-cell sequencing data.

### Generate allele-specific count matrices

`count_allelic.py` can be used to create allele-specific count matrices and .bam files. Required inputs are a filtered alignment, a phased .vcf file with heterozygous variants (for a single sample) as well as a cell index and regions file. Assumes diploid genome.

Requirements

- numpy (tested with 1.16.4)
- pysam (tested with 0.15.3)

      usage: count_allelic.py [-h] --cell_index CELL_INDEX --regions REGIONS --bam
                              BAM_FILE --vcf VCF_FILE --out_prefix OUT_PREFIX
                              [--save-sparse] [--output-bam]

      optional arguments:
        -h, --help            show this help message and exit
        --cell_index CELL_INDEX
                              Text file where each row corresponds to a cell barcode
                              to be considered.
        --regions REGIONS     Tab-separated file in .bed format with genomic regions
                              to be quantified.
        --bam BAM_FILE        Bam file with aligned reads. Must be sorted and
                              indexed.
        --vcf VCF_FILE        VCF file with phased, heterozygous variants. Must be
                              indexed.
        --out_prefix OUT_PREFIX
                              Prefix used for output files.
        --save-sparse         Output count matrices in sparse format. Highly
                              recommended!
        --output-bam          Output allele-specific bam files.
