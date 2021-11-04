#!/bin/env python
#
# Copyright 2013-2014 Graham McVicker and Bryce van de Geijn
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
#
# Modifications copyright (C) 2021 Tobias Heinen

"""
This program reads BAM files and counts the number of reads that match
either allele for all regions in the provided region file and all cell barcodes
from the cell index file. 

This script is a heavily modified version of the code originally developped for
the WASP pipeline:

    https://github.com/bmvdgeijn/WASP/blob/master/CHT/bam2h5.py
"""

import sys
import os
import gzip
import warnings

import argparse
import numpy as np

import pysam

# codes used by pysam for aligned read CIGAR strings
BAM_CMATCH     = 0 # M
BAM_CINS       = 1 # I
BAM_CDEL       = 2 # D
BAM_CREF_SKIP  = 3 # N
BAM_CSOFT_CLIP = 4 # S
BAM_CHARD_CLIP = 5 # H
BAM_CPAD       = 6 # P
BAM_CEQUAL     = 7 # =
BAM_CDIFF      = 8 # X

BAM_CIGAR_DICT = {0 : 'M',
                  1 : 'I',
                  2 : 'D',
                  3 : 'N',
                  4 : 'S',
                  5 : 'H',
                  6 : 'P',
                  7 : '=',
                  8 : 'X'}
unimplemented_CIGAR = [0, set()]

partially_mapped_reads = 0


def open_file(file_name, mode):
    """Opens compressed and uncompressed files."""
    if os.path.splitext(file_name)[1] == '.gz':
        return gzip.open(file_name, mode)
    else:
        return open(file_name, mode)


def is_indel(variant):
    """Checks if variant is an indel."""
    if (len(variant.alleles[0]) != 1) or (len(variant.alleles[1])) != 1:
        return True


def choose_overlap_variant(read, vcf_file):
    '''Picks out a single variant from those that the read overlaps.
    Returns a tuple containing 3 elements:
        [0] VariantRecord of chosen variant
        [1] Offset into read sequence
        [2] Flag indicating whether read overlaps known indel
    If there are no overlapping variants or the read cannot be processed,
    (None, None, False) is returned instead.
    '''
    read_variant_offsets = []
    variants = []

    block_start_read = 0
    block_start_genome = read.reference_start

    n_match_segments = 0
    overlap_indel = False

    for cig in read.cigartuples:
        op, op_len = cig

        if op == BAM_CMATCH:
            # this is a block of match/mismatch in read alignment
            block_end_read = block_start_read + op_len
            block_end_genome = block_start_genome + op_len

            # get offsets of any variants that this read overlaps
            new_variants = list(vcf_file.fetch(
                read.reference_name,
                block_start_genome,
                block_end_genome))
            for v in new_variants:
                read_variant_offsets.append(v.start - block_start_genome + block_start_read)
            variants.extend(new_variants)

            block_start_read = block_end_read
            block_start_genome = block_end_genome

            n_match_segments += 1
        elif op == BAM_CREF_SKIP:
            # spliced read, skip over this region of genome
            block_start_genome += op_len
        elif op == BAM_CSOFT_CLIP:
            # end of read is soft-clipped, which means it is
            # present in read, but not used in alignment
            block_start_read += op_len
        elif op == BAM_CINS:
            # dealing with insertion
            block_start_read += op_len
        elif op == BAM_CDEL:
            # dealing with deletion
            block_start_genome += op_len
        elif op == BAM_CHARD_CLIP:
            # end of read is hard-clipped, so not present
            # in read and not used in alignment
            pass
        else:
            unimplemented_CIGAR[0] += 1
            unimplemented_CIGAR[1].add(BAM_CIGAR_DICT[op])

            return (None, None, overlap_indel)

    # are any of the variants indels? If so, discard.
    for v in variants:
        if is_indel(v):
            overlap_indel = True
            return (None, None, overlap_indel)

    n_overlap_variants = len(read_variant_offsets)
    if n_overlap_variants == 0:
        # no variants overlap this read
        return (None, None, overlap_indel)

    if n_overlap_variants == 1:
        return (variants[0], read_variant_offsets[0], overlap_indel)

    # choose ONE overlapping variant uniformly at random
    r = np.random.randint(0, n_overlap_variants)
    return (variants[r], read_variant_offsets[r], overlap_indel)


def get_read_allele(read, vcf_file):
    """Determines the read allele (1 or 2) or returns None if inconclusive."""
    global partially_mapped_reads

    # TODO 
    # if start < 1 or end > chrom.length:
    #     sys.stderr.write('WARNING: skipping read aligned past end of '
    #                      'chromosome. read: %d-%d, %s:1-%d\n' %
    #                      (start, end, chrom.name, chrom.length))
    #     return


    if read.query_alignment_length!= read.query_length:
        partially_mapped_reads += 1
        return None

    # look for variants that overlap mapped read position, and if there
    # are more than one, choose one at random
    variant, read_offset, overlap_indel = choose_overlap_variant(read, vcf_file)

    if overlap_indel or variant is None:
        return None

    base = read.query_sequence[read_offset]
    [allele1 , allele2] = variant.samples[0].alleles
    if base == allele1:
        return 1
    if base == allele2:
        return 2
    else:
        # no match
        return None


def print_sparse_row(x, row, file):
    """Prints row vector in sparse format.

    For each element e with index i in x, this function prints a tab-separated
    entry (row, i, e) to the given file.
    
    Args:
        x: List of row values.
        row: Row number.
        file: File object.
    """
    for i, e in enumerate(x):
        if e > 0:
            print('\t'.join(str(x) for x in [row, i, e]), file=file)


def parse_args():
    """Parses command line arguments."""

    parser = argparse.ArgumentParser()
    parser.add_argument('--cell_index',
                        help=(
                            'Text file where each row corresponds to a cell '
                            'barcode to be considered.'),
                        metavar='CELL_INDEX',
                        required=True)
    parser.add_argument('--regions',
                        help=(
                            'Tab-separated file in .bed format with genomic '
                            'regions to be quantified.'),
                        metavar='REGIONS',
                        required=True)
    parser.add_argument('--bam',
                        help=(
                            'Bam file with aligned reads. Must be sorted and '
                            'indexed.'),
                        metavar='BAM_FILE',
                        required=True)
    parser.add_argument('--vcf',
                        help=(
                            'VCF file with phased, heterozygous variants. Must '
                            'be indexed.'),
                        metavar='VCF_FILE',
                        required=True)
    parser.add_argument('--out_prefix',
                        help='Prefix used for output files.',
                        metavar='OUT_PREFIX',
                        required=True)
    parser.add_argument('--save-sparse',
                        help=(
                            'Output count matrices in sparse format. Highly '
                            'recommended!'),
                        action='store_true',
                        default=False)
    parser.add_argument('--output-bam',
                        help='Output allele-specific bam files.',
                        action='store_true',
                        default=False)
    return parser.parse_args()


def main():
    args = parse_args()

    # parse cell barcodes
    with open(args.cell_index, 'r') as f:
        barcodes = [x.strip() for x in f.readlines()]

    # dict mapping barcodes to integers
    cell_barcode_dict = {}
    for x, cell in enumerate(barcodes):
        cell_barcode_dict[cell] = x
    n_barcodes = len(barcodes)

    # parse regions file
    regions_file = open_file(args.regions, 'rt')

    # create a txt file to hold the counts
    counts_allele1 = gzip.open(args.out_prefix + 'counts_allele1.txt.gz', 'wt+')
    counts_allele2 = gzip.open(args.out_prefix + 'counts_allele2.txt.gz', 'wt+')
    counts_total = gzip.open(args.out_prefix + 'counts_total.txt.gz', 'wt+')

    if not args.save_sparse:
        print('region' + '\t' + '\t'.join(barcodes), file=counts_allele1)
        print('region' + '\t' + '\t'.join(barcodes), file=counts_allele2)
        print('region' + '\t' + '\t'.join(barcodes), file=counts_total)

    sam_file = pysam.AlignmentFile(args.bam, 'rb')
    vcf_file = pysam.VariantFile(args.vcf)

    if args.output_bam:
        bam_allele1 = pysam.AlignmentFile(args.out_prefix + 'reads_allele1.bam', 'wb', template=sam_file)
        bam_allele2 = pysam.AlignmentFile(args.out_prefix + 'reads_allele2.bam', 'wb', template=sam_file)
        bam_total = pysam.AlignmentFile(args.out_prefix + 'reads_total.bam', 'wb', template=sam_file)

    print('Counting allele-specific read counts per region ...')
    print(' Regions processed: 0', end='\r')

    # iterate over regions
    region_cnt = 0
    for line in regions_file:
        region_cnt += 1
        if region_cnt % 100 == 0:
            print(' Regions processed: %d' % region_cnt, end='\r')

        region = line.split()
        region_name = '_'.join(region)
        region_counts_allele1 = [0] * n_barcodes
        region_counts_allele2 = [0] * n_barcodes
        region_counts_total = [0] * n_barcodes

        # fetch reads for region
        for read in sam_file.fetch(region[0], int(region[1]) - 1, int(region[2])):

            # determine allele for this read
            allele = get_read_allele(read, vcf_file)

            cell_barcode = read.query_name.split(':')[0]
            if allele == 1:
                try:
                    region_counts_allele1[cell_barcode_dict[cell_barcode]] += 1
                    if args.output_bam:
                        bam_allele1.write(read)
                except KeyError:
                    pass
            if allele == 2:
                try:
                    region_counts_allele2[cell_barcode_dict[cell_barcode]] += 1
                    if args.output_bam:
                        bam_allele2.write(read)
                except KeyError:
                    pass
            try:
                region_counts_total[cell_barcode_dict[cell_barcode]] += 1
                if args.output_bam:
                    bam_total.write(read)
            except KeyError:
                pass

        if args.save_sparse:
            print_sparse_row(
                region_counts_allele1, row=region_cnt-1, file=counts_allele1)
            print_sparse_row(
                region_counts_allele2, row=region_cnt-1, file=counts_allele2)
            print_sparse_row(
                region_counts_total, row=region_cnt-1, file=counts_total)
        else:
            print(
                '\t'.join([region_name] + [str(x) for x in region_counts_allele1]),
                file=counts_allele1)
            print(
                '\t'.join([region_name] + [str(x) for x in region_counts_allele2]),
                file=counts_allele2)
            print(
                '\t'.join([region_name] + [str(x) for x in region_counts_total]),
                file=counts_total)

    if args.output_bam:
        bam_allele1.close()
        bam_allele2.close()
        bam_total.close()
    sam_file.close()
    vcf_file.close()

    print(' Regions processed: %d' % region_cnt)
    print('Done.')

    # check if any of the reads contained an unimplemented CIGAR
    if unimplemented_CIGAR[0] > 0:
        print(
            'WARNING: Encountered %d instances of unimplemented CIGAR codes:\n' % unimplemented_CIGAR[0]
            + '[' + ','.join([str(c) for c in unimplemented_CIGAR[1]]) + ']\n'
            + 'Reads with these CIGAR codes were skipped.')

    if partially_mapped_reads > 0:
        print(
            'WARNING: %d partially mapped reads were skipped.\n'
            + 'Handling of partially mapped reads is not implemented.')


if __name__ == '__main__':
    main()



