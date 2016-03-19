"""

This script will take three inputs;
  - A CSV file containing CNV details, exported directly from Cartagenia
  - A pickle file which contains the chromosomal coordinates of all DDG2P genes
  - A BED file contining the probe sites of the array
  
The contents of these files will be cross referenced in several steps:
  - CNVs will be looked up against gene positions to find where partial duplications occur
  - These CNVs will be further investigated relative to specific base positions to find
    candidate genes to investigate for likely frameshifts
  - Where appropriate, the probe sites will be used to add confidence in a call
  - Based on exon positions and sizes, the possibility of a framshift will be calculated
  
This will bring several changes relative to the previous approach:
  - The degree of multiplication present in the CNV will not be factored in directly    
    - Unless the CNV is intragenic (which would require both ends to be intronic), only one additional
      copy will be included in the calculation
    - This prevents a 3x increased CNV being multiplied by 3 before being calculated, where only one 
      additional stretch of exons is likely to be included in the transcript
  - The coordinates of the CNV will not be used directly
  - Reduces the number of assumptions made in the process
  
This is also intended to generalise the process, as the BED file from the array will be used directly, making it interchangeable if a new set of probes are used
  
"""

import csv, os, cPickle
import numpy as np
import pandas as pd

def work_out_frameshift(chrom, start, stop, gene_dict, copy_num, gene):
    """ 
    Method to work out whether gene is frameshifted by partial duplication
    
    """
    # Process the CNV, finding out if a frameshift is likely to occur, and if it can be predicted with confidence
    fs_dictionary = calculate_exon_shift(chrom, start, stop, gene_dict, gene)
            
    # Process to get text output   
    output_strings = process_fs_dictionary(fs_dictionary, copy_num, gene, chrom)
    return output_strings
    
def process_fs_dictionary(fs_dict, copies, gene, chrom):
    """
    This method will read the contents of the fs dictionary, created in 'calculate_exon_shift()'
    The results of this process will be that the output files is written to with summary details
    Depending on opinion, this may be more appropriate as an excel document
    """
    transcripts = fs_dict[gene].keys()
    output_strings = []
    # Start with X chromosome instances, which will have results conditional on gender  
    # Conditional only applies where the CNV is intragenic
    output_strings.append('Chromosome {}, Gene {}\n'.format(chrom, gene))
    if chrom == 'X':
        # Deal with intragenic differently - all copies taken into account
        # Present conditional impact, based on sex
        output_strings.append('\tIf this patient is Male:\n')
        for transcript in transcripts:
            # If the CNV had an intronic end:
            if 'fail' in fs_dict[gene][transcript].keys():
                # Print fail statement and move on
                output_strings.append('\t\t{}\t{}\n'.format(transcript, fs_dict[gene][transcript]['fail']))
            elif fs_dict[gene][transcript]['fs']:
                # Special set of conditions based on intragenic CNVs
                if fs_dict[gene][transcript]['type'] == 'intragenic':
                    # Assume baseline is one X chrom
                    if (copies - 1) % 3 != 0:
                        # Use pre-calculated 'fs' value (True/False)
                        if fs_dict[gene][transcript]['fs']:
                            output_strings.append('\t\t{}\t{}-base frameshift\n'.format(transcript, fs_dict[gene][transcript]['frame change']))
                            # Intragenic will have two different 'confidence' reports
                            if fs_dict[gene][transcript]['confidence'] == 2:
                                output_strings.append("\t\t\t3' End: \t{}\n".format(fs_dict[gene][transcript]['confidence'][0]))
                                output_strings.append("\t\t\t5' End: \t{}\n".format(fs_dict[gene][transcript]['confidence'][1]))
                            else:
                                for x in fs_dict[gene][transcript]['confidence']:
                                    output_strings.append("\t\t\t{}\n".format(x))
                        else:
                            output_strings.append('\t\t{}\tNo frameshift\n'.format(transcript))

                    # If there are 3 additional copies, no frameshift will occur (probably)
                    else:
                        output_strings.append('\t\t{}\tNo frameshift\n'.format(transcript))
                else:
                    output_strings.append('\t\t{}\t{}-base frameshift\n'.format(transcript, fs_dict[gene][transcript]['frame change']))
            else:
                output_strings.append('\t\t{}\tNo Frameshift\n'.format(transcript))
                            
        output_strings.append('\n\tIf this patient is Female:\n')
        #if female and not a 3x increase
        for transcript in transcripts:
            # If the CNV had an intronic end:
            if 'fail' in fs_dict[gene][transcript].keys():
                # Print fail statement and move on
                output_strings.append('\t\t{}\t{}\n'.format(transcript, fs_dict[gene][transcript]['fail']))
            elif fs_dict[gene][transcript]['fs']:
                # Special set of conditions based on intragenic CNVs
                if fs_dict[gene][transcript]['type'] == 'intragenic':
                    # Assume baseline is one X chrom
                    if (copies - 2) % 3 != 0:
                        if fs_dict[gene][transcript]['fs']:
                            output_strings.append('\t\t{}\t{}-base frameshift\n'.format(transcript, fs_dict[gene][transcript]['frame change']))
                            # Intragenic will have two different 'confidence' reports
                            if fs_dict[gene][transcript]['confidence'] == 2:
                                output_strings.append("\t\t3' End{}\n".format(fs_dict[gene][transcript]['confidence'][0]))
                                output_strings.append("\t\t5' End{}\n".format(fs_dict[gene][transcript]['confidence'][1]))
                            else:
                                for x in fs_dict[gene][transcript]['confidence']:
                                    output_strings.append("\t\t{}\n".format(x))
                        else:
                            output_strings.append('\t\tTranscript: {} - No frameshift\n'.format(transcript))
                    else:
                        if fs_dict[gene][transcript]['fs']:
                            output_strings.append('\t\t{}\t{}-base frameshift\n'.format(transcript, fs_dict[gene][transcript]['frame change']))
                        else:
                            output_strings.append('\t\tTranscript: {} - No frameshift\n'.format(transcript))
                else:
                    output_strings.append('\t\t{}\t{}-base frameshift\n'.format(transcript, fs_dict[gene][transcript]['frame change']))
            else:
                output_strings.append('\t\t{}\tNo Frameshift\n'.format(transcript))
    # Work on Y c'some separately
    elif chrom == 'Y':
        # For each transcript
        for transcript in transcripts:
            # If the CNV had an intronic end:
            if 'fail' in fs_dict[gene][transcript].keys():
                # Print fail statement and move on
                output_strings.append('\t{}\t{}\n'.format(transcript, fs_dict[gene][transcript]['fail']))
            elif fs_dict[gene][transcript]['fs']:
                # identify intragenic situations
                if fs_dict[gene][transcript]['type'] == 'intragenic':
                    # Rule out frameshift using copy number
                    if (copies - 1) % 3 == 0:
                        output_strings.append('\t{}\tIntragenic, no frameshift\n'.format(transcript))
                    else:
                        if fs_dict[gene][transcript]['fs']:
                            output_strings.append('\t{}\t{}-base frameshift cause by intragenic CNV\n'.format(transcript,
                                                                                        fs_dict[gene][transcript]['frame change']))
                            # Intragenic will have two different 'confidence' reports
                            if fs_dict[gene][transcript]['confidence'] == 2:
                                output_strings.append("\t\t3' End: {}\n".format(fs_dict[gene][transcript]['confidence'][0]))
                                output_strings.append("\t\t5' End: {}\n".format(fs_dict[gene][transcript]['confidence'][1]))
                            else:
                                for x in fs_dict[gene][transcript]['confidence']:
                                    output_strings.append("\t{}\n".format(x))
                        else:
                            output_strings.append('\t{}\tIntragenic, no frameshift\n'.format(transcript))
                # All other types of frameshift
                else:
                    if (copies - 1) % 3 == 0:
                        output_strings.append('\t{}\t{} CNV, no frameshift\n'.format(transcript, fs_dict[gene][transcript]['type']))
                    # Actual frameshifts
                    else:
                        output_strings.append('\t\t{}\t{}-base frameshift due to {} CNV\n'.format(transcript,
                                                                                                  fs_dict[gene][transcript]['frame change'],
                                                                                                  fs_dict[gene][transcript]['type']))
                        output_strings.append("\t\t\t{}\n".format(fs_dict[gene][transcript]['confidence'][0]))
            else:
                output_strings.append('\t{}\tNo Frameshift\n'.format(transcript))

    else:
        # Autosomes
        # For each transcript
        for transcript in transcripts:
            # If the CNV had an intronic end:
            if 'fail' in fs_dict[gene][transcript].keys():
                # Print fail statement and move on
                output_strings.append('\t{}\t{}\n'.format(transcript, fs_dict[gene][transcript]['fail']))
            elif fs_dict[gene][transcript]['fs']:
                # identify intragenic situations
                if fs_dict[gene][transcript]['type'] == 'intragenic':
                    # Rule out frameshift using copy number
                    if (copies - 2) % 3 == 0:
                        output_strings.append('\t{}\tIntragenic, no frameshift\n'.format(transcript))
                    else:
                        if fs_dict[gene][transcript]['fs']:
                            output_strings.append('\t{}\t{}-base frameshift\n'.format(transcript, fs_dict[gene][transcript]['frame change']))
                            # Intragenic will have two different 'confidence' reports
                            # assert len(fs_dict[gene][transcript]['confidence']) == 2, 'Wrong number of confidence values for intragenic dup: {}'.format(fs_dict[gene][transcript]['confidence'])
                            # Print both ends of the confidence statement
                            if fs_dict[gene][transcript]['confidence'] == 2:
                                output_strings.append("\t\t3' End{}\n".format(fs_dict[gene][transcript]['confidence'][0]))
                                output_strings.append("\t\t5' End{}\n".format(fs_dict[gene][transcript]['confidence'][1]))
                            else:
                                for x in fs_dict[gene][transcript]['confidence']:
                                    output_strings.append("\t\t{}\n".format(x))
                        else:
                            output_strings.append('\t{}\tIntragenic, no frameshift\n'.format(transcript))
                # All other types of frameshift
                else:
                    if (copies - 2) % 3 == 0:
                        output_strings.append('\tTranscript: {} - {} CNV, no frameshift\n'.format(transcript, fs_dict[gene][transcript]['type']))
                    # Actual frameshifts
                    else:
                        output_strings.append('\t{}\t{}-base frameshift due to {} CNV\n'.format(transcript,
                                                                                                fs_dict[gene][transcript]['frame change'],
                                                                                                fs_dict[gene][transcript]['type']))
                        if fs_dict[gene][transcript]['frame change'] == 0:
                            print 'base change 0'
                            print
                        output_strings.append("\t\t{}\n".format(fs_dict[gene][transcript]['confidence'][0]))
            else:
                output_strings.append('\t{}\tNo Frameshift\n'.format(transcript))
    return output_strings

    
def calculate_exon_shift(chrom, start, stop, gene_dict, gene):
    """
    At this point, the gene_dict object has the structure:
    gene_dict[transcript] = {'exons': {exon_number: {'start': starts[x],
                                                     'stop': ends[x]},
                             'start': tranStart,
                             'stop': tranEnd,
                             'exon_num': exon_num,
                             'exon_list': range(1, exon_num + 1),
                             'strand': strand,
                                                          # Add details from the annotations dictionary
                             'disease': ddg2p_dict[chromosome][gene]['diseases']}
    may also have gene_dict['hi_score'] and ['hs_score'] if haploinsufficiency is present from DDD

    :return:
    """
    # Specific method for calculating the number of bases affected by FS
    out_dict = {gene: {}}
    for transcript in [transcript for transcript in gene_dict.keys()
                                   if transcript not in ['hi_score', 'hs_score']]:
        # Create a variable to determine if there is an overlap for this transcript
        # Not all transripts will be affected by a CNV due to varying lengths
        any_overlap = False
        cnv_type = ''
        bases_affected = 0
        exon_list = gene_dict[transcript]['exon_list']
        gene_start = gene_dict[transcript]['txStart']
        gene_stop = gene_dict[transcript]['txEnd']
        exons_affected = []
        exonic_end = False
        gap_found = False
        confidence = []
        if (gene_start < start) and (gene_stop > start and gene_stop < stop):
            # 3' end, exons overlapping start
            cnv_type = "3'"
            # Situation; gene starts outside CNV boundary
            last_end = gene_dict[transcript]['exons'][1]['stop']
            for exon in exon_list:
                # Get details from the exon entry
                exon_start = gene_dict[transcript]['exons'][exon]['start']
                exon_stop = gene_dict[transcript]['exons'][exon]['stop']
                # Record the end coordinate of the previous exon
                # This will establish the region between adjacent exons for the FS calculation
                
                
                # If the end point of the CNV is exonic, quit 
                if start >= exon_start and start <= exon_stop:
                    any_overlap = True
                    # End point is exonic, don't calculate
                    exonic_end = True
                    out_dict[gene][transcript] = {'fail': 'CNV end point is exonic, Frameshift not calculated'}
                    break
                elif start <= exon_start:
                    any_overlap = True
                    exons_affected.append(str(exon))
                    bases_affected += (exon_stop - exon_start)-1
                    if not gap_found:
                        gap_found = True
                        gap_start = last_end
                        gap_end = exon_start
                        # Do something to find probes in this region
                        confidence.append(get_confidence(chrom, last_end, exon_start, start, stop, False)) 
                else:
                    # Moves the 'end of last exon marker' along 
                    last_end = exon_stop
            
        elif (gene_start > start and gene_start < stop) and (gene_stop > stop):
            # 5' end, exons overlapping stop
            cnv_type = "5'"
            # Situation; gene will start within CNV region
            
            # Set initial value for the gap point
            last_end = gene_dict[transcript]['exons'][1]['stop']
            
            for exon in exon_list:
                exon_start = gene_dict[transcript]['exons'][exon]['start']
                exon_stop = gene_dict[transcript]['exons'][exon]['stop']
                
                # If the end point of the CNV is exonic, quit 
                if stop >= exon_start and stop <= exon_stop:
                    any_overlap = True
                    # End point is exonic, don't calculate
                    exonic_end = True
                    out_dict[gene][transcript] = {'fail': 'CNV end point is exonic, FS not calculated'}
                    break
                
                elif stop >= exon_stop:
                    any_overlap = True
                    exons_affected.append(str(exon))
                    bases_affected += (exon_stop - exon_start)-1
                    
                else:
                    if not gap_found:
                        gap_found = True
                        # Do something to find probes in this region
                        confidence.append(get_confidence(chrom, last_end, exon_start, start, stop, True)) 

        elif (gene_start > start) and (gene_stop < stop):
            # intragenic - use both ends independently
            cnv_type = "intragenic"
            # Situation; starts outside CNV, enters CNV, ends outside CNV
            
            # Initialise Boolean values
            outside = True
            inside = False
            
            # Set initial value for the gap point
            last_end = gene_dict[transcript]['exons'][1]['stop']
            
            for exon in exon_list:
                exon_start = gene_dict[transcript]['exons'][exon]['start']
                exon_stop = gene_dict[transcript]['exons'][exon]['stop']
                
                # If the end point of the CNV is exonic, quit 
                # Covers both start and stop
                if (stop >= exon_start and stop <= exon_stop) or \
                    (start >= exon_start and start <= exon_stop):
                    any_overlap = True
                    # End point is exonic, don't calculate
                    exonic_end = True
                    out_dict[gene][transcript] = {'fail': 'CNV end point is exonic, FS not calculated'}
                    break
                    
                elif outside:
                    if start > exon_stop:
                        last_end = exon_stop
                    elif start < exon_start:
                        any_overlap = True
                        outside = False
                        inside = True
                        exons_affected.append(str(exon))
                        bases_affected += (exon_stop - exon_start)-1
                        # Do the calculation of probe confidence here
                        confidence.append(get_confidence(chrom, last_end, exon_start, start, stop, False)) 
                        last_end = exon_stop
                    else:
                        print "Shouldn't reach this condition!"
                        this = raw_input
                elif inside:
                    if stop > exon_stop:
                        any_overlap = True
                        exons_affected.append(str(exon))
                        bases_affected += (exon_stop - exon_start)-1
                        last_end = exon_stop
                    elif stop < exon_start:
                        inside = False
                        # Do the calculation of probe confidence here   
                        confidence.append(get_confidence(chrom, last_end, exon_start, start, stop, True))                      
                    else:
                        print "Shouldn't reach this condition!"
                        this = raw_input
        if not any_overlap:
            out_dict[gene][transcript] = {'fail': 'No overlap for this transcript'}
        elif not exonic_end:
            modulo = bases_affected % 3
            if modulo != 0:
                frameshift = True
            else:
                frameshift = False
            
            out_dict[gene][transcript] = {'type': cnv_type,
                                          'exons aff': exons_affected,
                                          'bases aff': bases_affected,
                                          'frame change': bases_affected % 3,
                                          'fs': frameshift,
                                          'confidence': confidence,
                                          }
    return out_dict
    
def get_confidence(chrom, start, stop, cnv_start, cnv_stop, start_in):
    """
    This method will use a chromosome and span of bases as arguments to identify 
    all array probes within the region
    This will then use the CNV span in bases to determine whether there are probes 
    inside and outside the region, which can be used to infer confidence
    A message about the confidence will be returned.
    """
    
    # Use the chromosome number to select an appropriate region of the probe dataframe
    chromosome_of_interest = 'chr{}'.format(chrom)
    chrom_in_roi = df[df['chrom'] == chromosome_of_interest]
    
    # Probes selected by span of bases on that c'some
    probes_in_roi = chrom_in_roi[chrom_in_roi['start'] >= start]
    probes_in_roi = probes_in_roi[probes_in_roi['stop'] <= stop]
    # Can this be done as a one step assignment?
    
    if len(probes_in_roi) == 0:
        return 'No probes in this intronic region'
    elif len(probes_in_roi) == 1:
        return 'Only one intronic probe at CNV end'
    elif start_in:
        # Do the calculations
        inside = 0
        outside = 0
        
        # Bear in mind that the probes will be used to determine the calls;
        # CNV should never split a probe site, makes no sense
        
        # Sort values
        probes_in_roi.sort_values('start')
        
        # Take each probe in turn
        for index, probe in probes_in_roi.iterrows():
            if probe['stop'] < cnv_stop:
                inside += 1
            elif probe['start'] > cnv_stop:
                outside += 1
        if inside != 0 and outside != 0:
            return 'Significant confidence CNV is intronic'
        elif (inside == 0 and outside != 0):
            return 'Multiple intronic probes, but all are under CNV'
        elif (inside != 0 and outside == 0):
            return 'Multiple intronic probes, but all are outside CNV'
        else: 
            return 'Big fat error has occurred'
    else:
        # Do the calculations, assuming the start of the region is not in CNV
        inside = 0
        outside = 0
        
        # Sort values
        probes_in_roi.sort_values('start')
        
        # Take each probe in turn
        for index, probe in probes_in_roi.iterrows():
            if probe['stop'] < cnv_start:
                outside += 1
            elif probe['start'] > cnv_start:
                inside += 1
        if inside != 0 and outside != 0:
            return 'Significant confidence CNV is intronic'
        elif (inside == 0 and outside != 0):
            return 'Multiple intronic probes, but all are under CNV'
        elif (inside != 0 and outside == 0):
            return 'Multiple intronic probes, but all are outside CNV'
        else:
            return 'Big fat error has occurred'
            
    
# Establish file names
cases = 'fake_results.csv'
genes = 'ddg2p_with_exons.cpickle'
probes = 'DDG2P_probes.bed'
summary_out = 'Partial_Duplications.txt'
fs_out = 'frame_shift_summary.txt'

# Import the gene locations 
with open(genes, 'r') as handle:
    master_dict = cPickle.load(handle)
    
# Import the probe sites
# This uses the Pandas read_csv module, choosing the first 3 columns of the file and providing them with the names shown
df = pd.read_csv(probes, sep='\t', header=None, names=['chrom', 'start', 'stop'], usecols=[0,1,2])

# Go through the CNVs, one by one
with open(cases, 'r') as handle:
    reader = csv.DictReader(handle)
    with open(fs_out, 'w') as outhandle:
        for row in reader:
        
            # Select values from the results row
            chromosome = row['Chromosome']
            start = int(row['Start'])
            stop = int(row['Stop'])
            copy_num = int(row['Copy number'])
            cnv_type = row['Type']
            
            # Progress indicator
            print 'Sample: %s' % row['Sample']
            
            # Isolate the relevant part of the BED file based on the chromosome
            # These indexes are character values, not numbers (to account for chroms X & Y)
            dict_section = master_dict[chromosome]
            
            # For each row, review all genes on the chromosome
            for gene in dict_section:
            
                # Each transcript is reviewed-
                """
                This method isn't particularly efficient;
                    Scans through every gene and transcript first, to find if there are any overlaps
                    Then repeats the process to find the specific genes which the CNV overlaps
                    
                    This should be optimised, but works so quickly that it may not be worth refactoring
                """
                check = False
                intra = False
                for transcript in [transcript for transcript in dict_section[gene].keys()
                                   if transcript not in ['hi_score', 'hs_score']]:
                    # Use the first and last bases of the whole gene region
                    # This method doesn't account for gene direction, and may need to be reworked

                    gene_start = dict_section[gene][transcript]['txStart']
                    gene_stop = dict_section[gene][transcript]['txEnd']
                    
                    
                    # If there is partial overlap, add the gene list to the set (only one entry per transcript)
                    # 
                    #                   [---------------------EXON--------------------]
                    #           [---Clause 1---]
                    #                                  [---Clause 2---]
                    #                                                           [---Clause 3---]
                    
                    # first condition pair checks for partial overlaps
                    if (start <= gene_start and (stop >= gene_start and stop <= gene_stop)) \
                        or ((start >= gene_start and start <= gene_stop) and stop >= gene_stop):
                        # If a gene is found, add the details of the CSV, gene, and sample to the master set
                        check = True
                    # second block checks for an intragenic duplication, handled separately
                    # Isn't this re-checked later?
                    elif (start >= gene_start and start <= gene_stop) and (stop >= gene_start and stop < gene_stop):
                        check = True
                        intra = True
                        
                if check:
                    output_strings = work_out_frameshift(chromosome, start, stop, dict_section[gene], copy_num, gene)
                    outhandle.write('\nSample number: %s\n' % row['Sample'])
                    outhandle.write('\nCNV: {}:{}-{}\n'.format(chromosome, start, stop))
                    for line in output_strings:
                        outhandle.write(line)
                  

