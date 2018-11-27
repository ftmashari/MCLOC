# Protein Clustering Algorithm
This program clusters orthologous proteins. It needs Python 2.7 or higher, MCL, and HMMER 3.

**USAGE: MCLOC.py [-h] [-e EVA] [-s THR] [-c CHNK] p o**

## Positional Arguments
**p**	Directory containing all proteomes in fasta format. Each proteome should be in a separate file.
**o**	Path to the output file

## Optional Arguments
\-h, \-\-help            	show this help message and exit
\-e EVA, \-\-eva EVA     	E-value threshold for phmmer
\-s THR, \-\-thr THR     	Similarity threshold
\-c CHNK, \-\-chnk CHNK  	Chunk size
