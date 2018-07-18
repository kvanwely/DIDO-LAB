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

less_skip_total.py

Takes the group of exons undergoing less skipping as compared to the control situation, 
and calculates the base composition for each position (including upsteram and downstream
regions). 

more_skip_total.py

As for the previous script, but uses the exons that are skipped more as compared to controls.

less_skip_motif.py

Takes the group of exons undergoing lee skipping as compared to the control situation,
and calculates the frequency of K-mers in upstream, exonic, and downstream regions.

more_skip_motif.py

As for the previous script, but uses the exons that are skipped more as compared to controls. 

get_total_rnd2.py

Used to calculate base composition for a random subset of exons (from a .BED input).

get_motif_rnd2.py

Used to calculate the frequency of K-mer motifs in a random subset of exons (from a .BED input). 
