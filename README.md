This repository contains a collection of python scripts used to quantify exon skipping.
Scripts were developed for Python 2.7

Individual scripts in this repository:

flatten_gtf.py

A script to generate a flattened exon map in .BED format from a .GTF description file

exon_skip.py

A script that divides the expression of an exon (e.g. in FPKM) by the number of times
that exon is skipped.

FPKM_per_gene.py

A script to summarize FPKM for all exon sin a gene and use the same expression value 
for all exons in that gene. Requires exon with FPKM expression values and a gene map
in .BED format
