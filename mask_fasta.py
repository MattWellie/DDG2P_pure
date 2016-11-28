"""
Script to mask parts of a reference genome

This should be flexible enough to to use any line length in the FastA
This should also be able to take a CSV file as input to determine the regions

This will succeed if the end index exceeds the length of the given chromosome
"""
import argparse
from argparse import ArgumentParser


arg_parser = argparse.ArgumentParser(description='Specify the files to look at')
arg_parser.add_argument('-f', help='FastA file to rewrite', required=True, dest='fasta')
arg_parser.add_argument('-m', help='CSV file containing masking regions as chr,start,stop', required=True, dest='csv')
arg_parser.add_argument('-o', help='Name of the new fasta to write', required=True, dest='out')
args = arg_parser.parse_args()

def parse_csv(csv):

    regions_to_mask = {}

    with open(csv) as handle:

        # iterate over each line in the csv file and find regions to mask
        for line in handle:

            if line.rstrip() == '':
                continue

            linelist = line.split(',')
            if len(linelist) != 3:

                print('There is a line here which doesn\'t have 3 values')
                print(line)
                print(linelist)
                exit()

            # set values based on the csv values
            chrom = linelist[0].lower()
            start = int(linelist[1].rstrip())
            stop  = int(linelist[2].rstrip())

            if chrom not in regions_to_mask:
                regions_to_mask[chrom] = {}

            regions_to_mask[chrom][start] = stop

    return regions_to_mask

# counter for the genomic bases seen so far
base_count = 0
masking_now = False
next_end = 0
next_start = 0

# csv into a dict
regions_to_mask = parse_csv(args.csv)

# now process the file

# Each time you reach a new chromosome, get the next start point
# Load that start and end point into memory
# Once you reach the end of that region, delete the index and go looking for the next one
# This is aimed at needing fewer loops per line, better for huge files
with open(args.fasta) as in_handle:

    with open(args.out, 'w') as outhandle:

        current_chrom = ''
        current_chrom_positions = []
        ignore = False

        for line in in_handle:

            # FastA headers - print and use to find regions to mask
            if line[0] == '>':

                # print the chromosome indicator to the output file
                print >>outhandle, line.rstrip()

                # isolate the chromosome number/letter and shift to lower case
                current_chrom = line.replace('>chr', '').lower().rstrip()
                # print 'Current chromsome: {}'.format(current_chrom)

                # check if this chromosome has any regions of interest
                if current_chrom in regions_to_mask:

                    # if it does, get an ordered list of the start sites for the masking
                    current_chrom_positions = sorted(regions_to_mask[current_chrom].keys())
                    if len(current_chrom_positions) >= 1:

                        # set the script to stop ignoring the base sequence
                        ignore = False
                        next_start = int(current_chrom_positions[0])
                        next_end = int(regions_to_mask[current_chrom][next_start])

                        # print 'Next place to replace: {} - {}'.format(next_start, next_end)

                    else:
                        # run out of indices
                        # this shouldn't happen in this loop
                        ignore = True

                    # count from 0 for this chromosome
                    base_count = 0

                # otherwise set the script to ignore the rest of this chromosome
                else:
                    ignore = True

            # Line of bases but no masking
            elif ignore:

                print >>outhandle, line.rstrip()

            elif masking_now:

                length = len(line)
                base_count += length

                # print 'This line ends on {}'.format(base_count)
                if (base_count) >= next_end:

                    print >>outhandle, ('N'*length).rstrip()

                    masking_now = False

                    # delete the entry in the dictionary for this start point
                    del(regions_to_mask[current_chrom][next_start])
                    # print regions_to_mask

                    # get the indexes now available
                    current_chrom_positions = sorted(regions_to_mask[current_chrom].keys())
                    if len(current_chrom_positions) >= 1:

                        next_start = current_chrom_positions[0]
                        next_end = regions_to_mask[current_chrom][next_start]

                    # if there are no indexes left for this chromosome
                    else:

                        ignore = True

                # Otherwise this line needs to be masked
                else:

                    print >>outhandle, ('N'*length).rstrip()

            # count the bases on each line, decide if masking should start
            else:

                # get the length of the line and decide what to do

                length = len(line)
                base_count += length

                # if the next region starts before or at the end of this line
                # start masking
                if base_count >= next_start:

                    print >>outhandle, ('N'*length).rstrip()
                    masking_now = True

                    # And check if the masking should only happen for this line:
                    if base_count >= next_end:

                        # delete the entry in the dictionary for this start point
                        del(regions_to_mask[current_chrom][next_start])
                        # print regions_to_mask

                        # get the indexes now available
                        current_chrom_positions = sorted(regions_to_mask[current_chrom].keys())
                        if len(current_chrom_positions) >= 1:

                            next_start = current_chrom_positions[0]
                            next_end = regions_to_mask[current_chrom][next_start]

                            if base_count >= next_start:

                                # keep masking_now as true
                                pass

                            else:
                                masking_now = False

                        # Otherwise if there are no more indexes on this chromsome
                        else:
                            masking_now = False
                            ignore = True


                # if you shouldn't start masking yet, print and increase the count
                else:

                    print >>outhandle, line.rstrip()

print 'After all regions have been masked and deleted, these regions remain: '
print regions_to_mask

# elements will be retained in the dictionary if the region to mask extends past the end of the chromosome
# I'm choosing not to fix this as this may be useful in debugging
