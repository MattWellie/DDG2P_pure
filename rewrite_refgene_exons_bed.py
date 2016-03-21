"""
Script to reprocess the refGene file from UCSC to have exons separately on each row

initial columns:
    ['transcript', 'chrom', 'strand', 'txStart', 'txEnd', 'exonCount', 'exonStarts', 'exonEnds', 'gene']

This format is downloaded straight from the UCSC table browser, and is tab-delimited
The exonStarts and exonEnds columns are comma-separated, listing all start and stop points in order
- Regardless of gene strand/direction, coordinates are ascending

This file starts off about 62,000 lines long
This will be in the hundreds of thousands after processing...
The aim will be to cross-reference the DDG2P genes against this file, taking only those genes
 and transcripts of interest. In theory this can be done either before or afterwards, though
 doing it before would be more efficient. It shouldn't take too long either way, as it will
 only be a simple text manipulation
"""

import csv, os, cPickle
import numpy as np
import pandas as pd

# Set file names for input and output
input_file = 'refgene_and_exons_hg19.bed'
output_file = 're-sorted_refgene_exons.bed'

# File created in create_annotated_dict.py
ddg2p_annotations = 'annotated_master.cpickle'
ddg2p_with_exons = 'ddg2p_with_exons.cpickle'
ddg2p_with_exons_binary = 'ddg2p_with_exons_binary.cpickle'


# Haplo-insufficiency predictions from DDD project
hi_in = 'hi_predictions.txt'


# Grab the input contents (BED file) and read them into a dataframe structure
df = pd.read_csv(input_file, sep='\t')

# Read in the annotation dictionary
# This is the DDG2P file which contains features such as the diseases implicated, mode of inheritance, or HI scores
with open(ddg2p_annotations, 'rb') as pickle_handle:
    ddg2p_dict = cPickle.load(pickle_handle)

# Open the output file as a handle to write into
with open(output_file, 'wb') as handle:
    """
    This loop will do two things:
    - Re-write the BED file to show each exon as a separate row
    - Create a dictionary of the BED file contents for use later on
    """

    output_dictionary = {}

    # export the row of titles as the first line
    handle.write('chrom\tstart\tstop\ttranscript\tgene\tstrand\ttxStart\ttxEnd\n')
    for index, row in df.iterrows():

        # Get values from the row using column header values
        chromosome = row['chrom'].replace('chr', '')
        transcript = row['transcript']
        strand = row['strand']
        gene = row['gene']
        exon_num = int(row['exonCount'])
        tranStart = row['txStart']
        tranEnd = row['txEnd']

        # A condition to see if it is a DDG2P gene, based on presence in DDG2P annotation file
        if gene in ddg2p_dict[chromosome].keys():

            # Handle single-exon genes separately
            if exon_num == 1:
                # Use a list structure to prevent need for further downstream modifications
                # The '-1' is used to remove the trailing comma present on all lists
                starts =[row['exonStarts'][:-1]]
                ends = [row['exonEnds'][:-1]]
            else:
                # The '-1' is used to remove the trailing comma present on all lists
                # There was originally a try-catch here, but it was never activated
                # The split into multi- and single-exon genes and removal of trailing comma works
                starts = row['exonStarts'][:-1].split(',')
                ends = row['exonEnds'][:-1].split(',')

            # For each exon, write a new row to the output BED file
            for x in range(exon_num):
                handle.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome,
                                                                       starts[x],
                                                                       ends[x],
                                                                       transcript,
                                                                       gene,
                                                                       strand,
                                                                       tranStart,
                                                                       tranEnd))

                # Now write the same content to the dictionary

                # These lines will be executed if the dictionary doesn't already have these keys
                if chromosome not in output_dictionary:
                    output_dictionary[chromosome] = {}
                if gene not in output_dictionary[chromosome].keys():
                    output_dictionary[chromosome][gene] = {}

                if transcript not in output_dictionary[chromosome][gene].keys():
                    output_dictionary[chromosome][gene][transcript] = {'exons'     : {},
                                                                       'start'     : int(tranStart),
                                                                       'stop'      : int(tranEnd),
                                                                       'exon_num'  : exon_num,
                                                                       'exon_list' : range(1, exon_num + 1),
                                                                       'strand'    : strand,
                                                                       'txStart'   : int(tranStart),
                                                                       'txEnd'     : int(tranEnd)
                                                                      }

                    # Add details from the annotations dictionary where available
                    # Should be available for all of them, as they are DDG2P
                    if gene in ddg2p_dict[chromosome]:
                        output_dictionary[chromosome][gene][transcript]['disease'] = \
                            ddg2p_dict[chromosome][gene]['diseases']

                # Add the exon coordinates under an appropriate index
                # Using range for Python indexing would mean exons start at 0, unintuitive
                output_dictionary[chromosome][gene][transcript]['exons'][x + 1] = {'start': int(starts[x]),
                                                                                   'stop' : int(ends[x])}

# Now add HI annotations where available
with open(hi_in, 'r') as handle:
    for line in handle:
        # Break up the details in the HI file, and add to appropriate genes in dict
        row_list = line.replace('_', ' ').split('\t')
        chromosome = row_list[0].replace('chr', '')
        detail_list = row_list[3].split('|')
        gene_symbol = detail_list[0]
        hi_score = float(detail_list[1])
        hs_score = detail_list[2]
        if chromosome in output_dictionary:
            if gene_symbol in output_dictionary[chromosome]:
                output_dictionary[chromosome][gene_symbol]['hi_score'] = hi_score
                output_dictionary[chromosome][gene_symbol]['hs_score'] = hs_score

print 'Exons covered: {}'.format(output_dictionary.keys())

# 2 separate writes, some binary conflict issues on reading files in Mac OSX
with open(ddg2p_with_exons_binary, 'wb') as handle:
    cPickle.dump(output_dictionary, handle)

with open(ddg2p_with_exons, 'w') as handle:
    cPickle.dump(output_dictionary, handle)
