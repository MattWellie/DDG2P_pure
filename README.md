 DDG2P_pure
 ==========
 
A program to interpret genetic duplications (CNV) with reference to a set of array probes to identify likely frameshifts. This is done to add an additional component to the interpretation of duplications, as a duplication leading to a frameshift or degraded transcript could be interpreted as a loss of a funtional allele

To run this analysis, the tool uses several key components:

1.  A BED file containing the exons of any genes relevant to the analysis (currently uses a dump of all hg19 genes)
2.  A BED file of probe sites for the array being used
3.  A CSV containing the CNV results to be analysed
4.  (optional) Gene annotations (e.g. DDG2P supplement files)
