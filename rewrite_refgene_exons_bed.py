"""
Script to reprocess the refGene file from UCSC to have exons separately on each row

initial columns:
    ['transcript', 'chrom', 'strand', 'exonCount', 'exonStarts', 'exonEnds', 'gene']

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


# Grab the input contents and read them into a dataframe structure
df = pd.read_csv(input_file, sep='\t')

# Read in the annotation dictionary
with open(ddg2p_annotations, 'rb') as pickle_handle:
    ddg2p_dict = cPickle.load(pickle_handle)

# Open the output file as a handle
with open(output_file, 'wb') as handle:
    # This loop will do two things:
        # Re-write the BED file to show each exon as a separate row
        # Create a dictionary of the BED file contents for use later on

    transcripts = set()
    dict = {}

    # export the row of titles as the first line
    handle.write('chrom\tstart\tstop\ttranscript\tgene\tstrand\ttxStart\ttxEnd\n')
    for index, row in df.iterrows():

        # Get values from the row
        chromosome = row['chrom'].replace('chr', '')
        transcript = row['transcript']
        strand = row['strand']
        gene = row['gene']
        exon_num = int(row['exonCount'])
        tranStart = row['txStart']
        tranEnd = row['txEnd']


        # Put in a condition here to see if it is a DDG2P gene
        if gene in ddg2p_dict[chromosome].keys():

            transcripts.add(transcript)
            if exon_num == 1:
                starts = row['exonStarts'][:-1]
                ends = row['exonEnds'][:-1]
            else:
                try:
                    starts = row['exonStarts'][:-1].split(',')
                    ends = row['exonEnds'][:-1].split(',')
                except AttributeError:
                    print 'starts: {}'.format(row['exonStarts'])
                    print 'ends: {}'.format(row['exonEnds'])
                    this = raw_input('Error finding exon boundaries!')

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
                if chromosome not in dict:
                    dict[chromosome] = {}
                if gene not in dict[chromosome].keys():
                    dict[chromosome][gene] = {}
                if transcript not in dict[chromosome][gene].keys():
                    dict[chromosome][gene][transcript] = {'exons': {},
                                                          'start': int(tranStart),
                                                          'stop': int(tranEnd),
                                                          'exon_num': exon_num,
                                                          'exon_list': range(1, exon_num + 1),
                                                          'strand': strand,
                                                          'txStart': int(tranStart),
                                                          'txEnd': int(tranEnd),
                                                          # Add details from the annotations dictionary
                                                          'disease': ddg2p_dict[chromosome][gene]['diseases']
                                                          }
                # Add the exon coordinates under an appropriate index
                dict[chromosome][gene][transcript]['exons'][x + 1] = {'start': int(starts[x]),
                                                                      'stop': int(ends[x])}

with open(hi_in, 'r') as handle:
    for line in handle:
        row_list = line.replace('_', ' ').split('\t')
        chromosome = row_list[0].replace('chr', '')
        detail_list = row_list[3].split('|')
        gene_symbol = detail_list[0]
        hi_score = float(detail_list[1])
        hs_score = detail_list[2]
        if chromosome in dict:
            if gene_symbol in dict[chromosome]:
                dict[chromosome][gene_symbol]['hi_score'] = hi_score
                dict[chromosome][gene_symbol]['hs_score'] = hs_score

# 2 separate writes, some binary conflict issues on reading files in Mac OSX
with open(ddg2p_with_exons_binary, 'wb') as handle:
    cPickle.dump(dict, handle)

with open(ddg2p_with_exons, 'w') as handle:
    cPickle.dump(dict, handle)
