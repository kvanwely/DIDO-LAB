
#this script will take an input gene model (.BED) and FPKM output as generated by RSeQR fragcount 
#each value in the input that overlaps with the gene and strand will be added
#output will be reported as a set of values equal for every exon in a single gene

import sys
from argparse import ArgumentParser, FileType
import os
import string

# get reference genes from bed file and format output genes for printing

def get_fpkm_exons (ref_chrom,fpkm_file,fpkm_exons):
    del fpkm_exons[:]
    fpkmfile = open(fpkm_file,'r')
    for line in fpkmfile:       
        line = line.rstrip('\n')
        if line.startswith('#'):
            continue
        try:
            exon_data = line.split('\t')
        except ValueError:
            continue
            
        if exon_data[0] == ref_chrom:
            tmp_exon = [exon_data[0],exon_data[1],exon_data[2],exon_data[3],exon_data[4],exon_data[5],exon_data[6],exon_data[7],exon_data[8]]
            fpkm_exons.append(tmp_exon)
           
    fpkmfile.close()       
    return fpkm_exons
 
def go_gene1 (fpkm_file,gene_file):

    gene1 = []
    gene2 = []
    exon = []
    gene1_data = []
    test_exon = []
    gene1_list = []
    gene2_list = []
    fpkm_exons = [] 
    prev_chrom = ''
    num_exons = 0
        
    ref_file = open(gene_file)
    print >>sys.stderr, "Getting gene regions ..."
    
    for line in ref_file:
        line = line.rstrip('\n')
        if line.startswith('#'):
            continue
        try:
            gene1_data = line.split('\t')
        except ValueError:
            continue 
        gene1 = [gene1_data[0],gene1_data[1],gene1_data[2],gene1_data[3],gene1_data[4],gene1_data[5]]
        gene1_list.append(gene1)
    
    fpkmfile = open(fpkm_file,'r')
    out_file_tmp = fpkm_file.split(".")
            
    out_file = out_file_tmp[0] + '.G1_FPKM.xls'
    outfile = open(out_file,'w')
    
    for gene1 in gene1_list:
        
        gene1_fpkm = 0.0
        gene1_frag = 0.0
        gene1_fpm = 0.0
        ref_chrom = gene1[0].rstrip("\t")

        if ref_chrom != prev_chrom:
            get_fpkm_exons(ref_chrom,fpkm_file,fpkm_exons)
            prev_chrom = ref_chrom
            
        for exon in fpkm_exons:
            if exon[3] == gene1[3] and exon[5] == gene1[5]:
                if int(exon[1]) >= int(gene1[1]) - 2 and int(exon[2]) <= int(gene1[2]) + 2:
                    gene1_frag = gene1_frag + float(exon[6])
                    gene1_fpm = gene1_fpm + float(exon[7])
                    gene1_fpkm = gene1_fpkm + float(exon[8])
                    num_exons = num_exons + 1
       
        gene2 = [gene1[0],gene1[1],gene1[2],gene1[3],gene1[4],gene1[5],str(gene1_frag),str(gene1_fpm),str(gene1_fpkm)]
        gene2_list.append(gene2)
        print >>outfile, '\t'.join((gene2[0],gene2[1],gene2[2],gene2[3],gene2[4],gene2[5],str(int(gene1_frag)),str(gene1_fpm),str(gene1_fpkm)))
        print >>sys.stderr, "Total exons : " + str(num_exons) + "          \r",
        
    outfile.close ()
        
    print >>sys.stderr, "\n", 

    out_file = out_file_tmp[0] + '.G2_FPKM.xls'
    outfile = open(out_file,'w')  
    print >>outfile, '\t'.join(('#chrom','st','end','accession','mRNA_size','gene1_strand','Frag_count','FPM','FPKM'))
    num_exons = 0
      
    for line in fpkmfile:  
        prev_chrom = ""
        exon_frag = ""
        exon_fpm = ""
        exon_fpkm = ""    
        line = line.rstrip('\n')
        if line.startswith('#'):
            continue
        try:
            exon = line.split('\t')
        except ValueError:
            continue
            
        if exon[0] != prev_chrom:
            del gene1_list[:]
            for gene2 in gene2_list:
                if exon[0] == gene2[0]:
                    gene1_list.append(gene2) 
            prev_chrom = exon[0]
            
        for gene1 in gene1_list:
            if exon[3] == gene1[3] and exon[5] == gene1[5]:
                 if int(exon[1]) >= int(gene1[1]) - 2 and int(exon[2]) <= int(gene1[2]) + 2:  
                     exon_frag = gene1[6]
                     exon_fpm = gene1[7]
                     exon_fpkm = gene1[8]
        
        print >>outfile, '\t'.join((exon[0],exon[1],exon[2],exon[3],exon[4],exon[5],exon_frag,exon_fpm,exon_fpkm))
        
        num_exons = num_exons + 1
        print >>sys.stderr, "Done  exons : " + str(num_exons) + "               \r", 

    outfile.close ()
    print >>sys.stderr, "\n", 
    print >>sys.stderr, "All done !"
    
# argument parser
if __name__ == '__main__':
    parser = ArgumentParser(
        description='get gene1 names for represented exons')
    parser.add_argument('fpkm_file',nargs='?',type=str, help='input exons') 
    parser.add_argument('gene1_file',nargs='?',type=str, help='reference gene1 exons')       
    args = parser.parse_args()
    if not (args.fpkm_file and args.gene1_file):
            parser.print_help()
            exit(1)
    go_gene1 (args.fpkm_file,args.gene1_file)
    
        