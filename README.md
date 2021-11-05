## scAI (single-cell allelic imbalance)
Scripts for allele-specific analyses of single-cell sequencing data.

### Generate allele-specific count matrices

`count_allelic.py` can be used to create allele-specific count matrices and .bam files. Required inputs are a filtered alignment, a phased .vcf file with heterozygous variants (for a single sample) as well as a cell index and regions file. Cell identifiers are assumed to be stored under the 'CB' tag in the input .bam file. Only diploid genomes are supported. In order to obtain reliable estimates of allele-specific signals, we recommend filtering for mapping biases using [WASP](https://github.com/bmvdgeijn/WASP/tree/master/mapping).

Requirements

- numpy (tested with 1.16.4)
- pysam (tested with 0.15.3)

      usage: count_allelic.py [-h] --cell_index CELL_INDEX --regions REGIONS --bam
                              BAM_FILE --vcf VCF_FILE --out_prefix OUT_PREFIX
                              [--save_sparse] [--output_bam]

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
        --save_sparse         Output count matrices in sparse format. Highly
                              recommended!
        --output_bam          Output allele-specific bam files.
        
The script will produce three matrices

- `counts_allele1.txt.gz` Region-by-cell matrix of counts for allele one.
- `counts_allele2.txt.gz` Region-by-cell matrix of counts for allele two.
- `counts_total.txt.gz` Region-by-cell matrix of counts for all reads (allele-specific and unassigned).

If `--save_sparse`, these matrices are stored as tab-separated files using a sparse triplet format. In this case each line denotes a non-zero entry of the count matrix of the form (region ID, cell ID, count). If `--output_bam`, this script additionally produces matching .bam files containing the reads used for the quantifiations in each matrix.
