## DDG2P_pure
 
 
A program to interpret genetic duplications (CNV) with reference to a set of array probes to identify likely frameshifts. This is done to add an additional component to the interpretation of duplications, as a duplication leading to a frameshift or degraded transcript could be interpreted as a loss of a functional allele

To run this analysis, the tool uses several key components:

1.  A BED file containing the exons of any genes relevant to the analysis (currently uses a dump of all hg19 genes)
2.  A BED file of probe sites for the array being used
3.  A CSV containing the CNV results to be analysed
4.  (optional) Gene annotations (e.g. DDG2P supplement files)

The software requirements for this program are defined in req.txt, and were exported from the conda virtual environment the application was developed in.

The program requires a CSV file containing CNV details, a pickle file with all genes of interest (DDG2P in this case) and a BED file containing probe sites. These are currently specified around line 509 of `frameshift_comp_with_probes.py`, though this will be changed to argparse in a future release.

This repository also contains a tool for manually masking regions of a reference genome (FastA). This will be shifted into its own repository soon
