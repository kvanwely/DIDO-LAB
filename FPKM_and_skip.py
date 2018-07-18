#a single script to combine FPKM counting and exon skipping
#uses Tophat2 .BAM coverage and .BED jucntions
#requires a flattend gene mode in .BED format as guide

'''--------------------------------------------------------------------------------
calculate exon skipping for each region in a BED file.
--------------------------------------------------------------------------------'''

#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
     print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " needs python2.7!\n"
     sys.exit()

import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
import sets
from time import strftime

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

#import my own modules
from qcmodule import SAM
from qcmodule import bam_cigar
#changes to the paths

#changing history to this module

__author__ = "Karel van Wely"
__credits__ = "Liguo Wang for the FPKM_count script"
__license__ = "GPL"
__version__="0.1 alpha"
__maintainer__ = " "
__email__ = " "
__status__ = "Pre-production"

def fetch_cigar_exon(chrom, st, cigar):
	''' fetch exon regions defined by cigar. st must be zero based return list of tuple of (chrom,st, end) '''
	chrom_st = st
	exon_bound = []
	for c,s in cigar:	#code and size
		if c==0:		#match
			exon_bound.append((chrom, chrom_st,chrom_st + s))
			chrom_st += s
		elif c==1:		#insertion to ref
			continue
		elif c==2:		#deletion to ref
			chrom_st += s
		elif c==3:		#gap or intron
			chrom_st += s
		elif c==4:		#soft clipping. We do NOT include soft clip as part of exon
			chrom_st += s
		else:
			continue
	return exon_bound

def build_range(refgene):
    '''build ranges for exonic region'''
    ranges={}
    refgenefile = open(refgene,'r')
    for line in refgenefile:
        try:
            if line.startswith(('#','track','browser')):continue  
            # Parse fields from gene tabls
            fields = line.split()
            chrom = fields[0].upper()
            exon_starts = int(fields[1]) + 1          # bed files are zero-based
            exon_ends = int(fields[2])
            gene_name = fields[3]
            strand = fields[5].replace(" ","_")
        except:
            print >>sys.stderr," skipped this line: " + line,
            continue
          
        if chrom not in ranges:
            ranges[chrom] = Intersecter()
            ranges[chrom].add_interval( Interval( exon_starts, exon_ends) )
            
    refgenefile.close()
    return ranges

def get_fpkm_exons (ref_chrom,fpkm_file,fpkm_exons):
    ''' open an existing FPKM file and extract data '''
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

def open_exons (ref_file,ref_exons):
    ''' open a BED file  for coordinates '''
    del ref_exons[:]
    reffile = open(ref_file)
    for line in reffile:
        if line.startswith('#'):
            continue
        exon = line.rstrip('\n')
        ref_exons.append(exon)
    reffile.close ()
    return ref_exons

def open_junctions (do_chrom,in_file,test_junctions):
    ''' get junctions for do_chrom from tophat2 junctions.bed file and put values in same order as exons '''
    del test_junctions[:]
    infile = open(in_file)
    for line in infile:       
        line = line.rstrip('\n')
        if line.startswith('#'):
            continue
        try:
            test_chrom, test_left, test_right, test_name, test_dept, test_strand, someval_1, someval_2, colorstuff, someval_3, overhangs, someval_4  = line.split('\t')
        except ValueError:
            continue
        try:
            overhang_left, overhang_right = overhangs.split(',')
        except ValueError:
            continue
        test_strand = test_strand[0]
        if test_chrom != do_chrom:
            continue
        else:
            junction = [test_chrom,test_left,test_right,test_name,test_dept,test_strand,overhang_left,overhang_right]
            test_junctions.append(junction)
    infile.close ()
    
    test_junctions = sorted (test_junctions,key=lambda item: int(item[1])) 
    return test_junctions

def go_gene (fpkm_file,gene_file):
    ''' calculate expression data per entire gene, outputs a gene and exon list ''' 
    gene1 = []
    gene2 = []
    gene3 = []
    exon = []
    gene1_data = []
    test_exon = []
    gene1_list = []
    gene2_list = []
    fpkm_exons = [] 
    prev_chrom = ''
    num_exons = 0
        
    ref_file = open(gene_file)
    print >>sys.stderr, "Getting gene regions from " + gene_file + " ... "
    
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
    
    out_file_tmp = fpkm_file.split(".")            
    out_file = out_file_tmp[0] + '.GENE_1.xls'
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
        print >>sys.stderr, "Counting exons for genes : " + str(num_exons) + "          \r",
        
    outfile.close ()
    print >>sys.stderr, "\n", 

    num_exons = 0
    fpkmfile = open(fpkm_file,'r')
    out_file_tmp = fpkm_file.split(".")
    out_file = out_file_tmp[0] + '.GENE_2.xls'
    outfile = open(out_file,'w')
    
    for line in fpkmfile:  
        prev_chrom = ""
        exon_frag = "-1"
        exon_fpm = "-1"
        exon_fpkm = "-1"    
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
            
        gene1_list = sorted (gene1_list,key=lambda item: int(item[1]))
        
        for gene1 in gene1_list:
            if exon[3] == gene1[3] and exon[5] == gene1[5]:
                 if int(exon[1]) >= int(gene1[1]) - 2 and int(exon[2]) <= int(gene1[2]) + 2:  
                     exon_frag = gene1[6]
                     exon_fpm = gene1[7]
                     exon_fpkm = gene1[8]
                
        print>> outfile, '\t'.join((exon[0],exon[1],exon[2],exon[3],exon[4],exon[5],exon_frag,exon_fpm,exon_fpkm))
        num_exons = num_exons + 1
        print >>sys.stderr, "Adding expr data to list : " + str(num_exons) + "                 \r", 
    
    # clean up 
    
    fpkmfile.close()
    outfile.close()
    del gene1_list[:]
    del gene2_list[:] 
    print >>sys.stderr,"\n"
    return()
    
def skip_counter(distance,min_skip_count,ref_file,in_file):
    ''' tests how many times each exon in ref_file is skipped, skips are loaded from in_file'''
    exon = []
    junction = []
    ref_exons = [exon]
    test_junctions = [junction]
    prev_chrom = 'abc'
    ref_chrom = 'cba'
      
    out_file_tmp = ref_file.split(".")
    out_file = out_file_tmp[0] + '.SKIP.xls'
    out2file = open(out_file,'w')
    print >>sys.stderr, "Loading junctions, please be patient ... \n"
    
    open_exons (ref_file,ref_exons)
    for exon in ref_exons:
        out_line = ''
        skip_count = 0
                
        try:
            exon_data = exon.split('\t')
        except ValueError:
            continue
        
        ref_chrom = exon_data[0] 
        ref_left = exon_data[1]
        ref_right = exon_data[2]
        ref_name = exon_data[3]
        ref_dept = exon_data[4] 
        ref_strand = exon_data[5] 
        frag_count = exon_data[6]
        FPM_count = exon_data[7]  
        FPKM_count = exon_data[8]    
        ref_strand = ref_strand[0]
        ref_leftcoord = int(ref_left)
        ref_rightcoord = int(ref_right)
           
        if ref_chrom != prev_chrom:
            open_junctions (ref_chrom,in_file,test_junctions)
            prev_chrom = ref_chrom
        
        test_left_coord = 0
        test_junc_index = 0
        
        while test_left_coord < ref_rightcoord and len(test_junctions) > test_junc_index:
            junction = test_junctions[test_junc_index]
            test_chrom = junction[0]
            test_leftcoord = int(junction[1]) + int(junction[6])
            test_rightcoord = int(junction[2]) - int(junction[7])
            test_name = junction[3]
            test_dept_val = int(junction[4])
            test_strand = junction[5]
            
            test_junc_index = test_junc_index + 1
              
            if (test_rightcoord - test_leftcoord) > distance:
                continue
              
            if ref_strand == test_strand and test_leftcoord <= ref_leftcoord and test_rightcoord >= ref_rightcoord:
                skip_count = skip_count + test_dept_val
            
            if test_rightcoord < (ref_leftcoord - (5 * distance)):     # bit of safe distance before deleting "old" junctions
                del test_junctions[test_junc_index-1]
                test_junc_index = 0

        FPKM_val = float(FPKM_count)
        skip_val = float(skip_count)
        if FPKM_val > 0.0 and skip_count >= min_skip_count:
            SKIP_ratio = str(skip_val / FPKM_val)
        elif skip_count >= min_skip_count:
            SKIP_ratio = 'inf'
        else:
            SKIP_ratio = "0.0"
           
        out_line = exon + '\t' + str(skip_count) + '\t' + SKIP_ratio
        out2_line = out_line.replace('.',',')
        print >>out2file,out2_line
        out_stuff = test_chrom + "  :  " + ref_name + "                                  \r"
        print >>sys.stderr, out_stuff,
      
    out2file.close()
    print >>sys.stderr, "\n"

def main():
    usage="%prog [options]" + '\n' + __doc__ + "\n"
    parser = OptionParser(usage,version="%prog ")
    parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM format (SAM is not supported). [required]")
    parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
    parser.add_option("-r","--refgene",action="store",type="string",dest="refgene_bed",help="Reference gene model in bed fomat. [required]")
    parser.add_option("-j","--junctions",action="store",type="string",dest="junctions_bed",help="Junctions of skips in bed format.")
    parser.add_option("-d","--exon_distance",action="store",type="int",dest="do_dist",default=50000,help="Maximum skipping distance. default=%default")
    parser.add_option("-t","--threshold",action="store",type="int",dest="min_skips",default=3,help="Minimum skipping threshold. default=%default")
    parser.add_option("-u","--skip-multi-hits",action="store_true",dest="skip_multi",help="How to deal with multiple hit reads. Presence this option renders program to skip multiple hits reads.")
    parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be called \"uniquely mapped\". default=%default")
    parser.add_option("-s","--skipcount",action="store",type="string",dest="skip_FPKMcount",help="Only do the exon skipping analysis, requires FPKM file as input")
    parser.add_option("-g","--genebased",action="store",type="string",dest="gene_based",help="Do gene based caluclations for skipping (requires additional gene index bed)")
    
    (options,args)=parser.parse_args()
            
    if not ((options.output_prefix and options.input_file and options.refgene_bed) or (options.skip_FPKMcount and options.junctions_bed)):
        parser.print_help()
        sys.exit(0)
    
    tmp_out_str = str(options.map_qual)                    
    print >>sys.stderr, "minimum map quality : " + tmp_out_str
    tmp_out_str = str(options.min_skips)                    
    print >>sys.stderr, "minimum skip thresh : " + tmp_out_str
    tmp_out_str = str(options.do_dist)                    
    print >>sys.stderr, "maximum skip distnc : " + tmp_out_str
    
    if options.skip_FPKMcount and not options.gene_based:
        print >>sys.stderr, "skipping counting phase through -s option"
        if not (os.path.exists(options.skip_FPKMcount) and os.path.exists(options.junctions_bed)):
            print >>sys.stderr, "cannot find FPKM or junctions file"
            sys.exit(0)
        else:
            skip_counter(options.do_dist,options.min_skips,options.skip_FPKMcount,options.junctions_bed)
            print >>sys.stderr, "All done !\n'"
            sys.exit(0)
        
    if options.skip_FPKMcount and options.gene_based:
        print >>sys.stderr, "skipping counting phase through -s option"
        if not (os.path.exists(options.gene_based) and os.path.exists(options.skip_FPKMcount) and os.path.exists(options.junctions_bed)):
            print >>sys.stderr, "cannot find gene, junctions, or FPKM file"
            sys.exit(0)
        else:
            go_gene(options.skip_FPKMcount,options.gene_based)
            out_file_tmp1 = options.skip_FPKMcount
            out_file_tmp2 = out_file_tmp1.split(".")            
            out_file_tmp3 = out_file_tmp2[0] + '.GENE_2.xls'
            skip_counter(options.do_dist,options.min_skips,out_file_tmp3,options.junctions_bed)
            print >>sys.stderr, "All done !\n'"
            sys.exit(0)
                           
    for file in (options.input_file, options.refgene_bed):
        if not os.path.exists(file):
            print >>sys.stderr, file + " does NOT exists" + '\n'
            sys.exit(0)
            
    if not os.path.exists(options.input_file + '.bai'):
        print >>sys.stderr, "cannot find index file of input BAM file"
        print >>sys.stderr, options.input_file + '.bai' + " does not exists"
        sys.exit(0)
            
    obj = SAM.ParseBAM(options.input_file) 
    out1_file = options.output_prefix + '.FPKM.xls'
    out1file = open(out1_file,'w')
     
    #++++++++++++++++++++++++++++++++++++counting fragments
    print >>sys.stderr, "Extract exon regions from "+ options.refgene_bed + ' ...'
    gene_ranges = build_range(options.refgene_bed)
    
     #++++++++++++++++++++++++++++++++++++++++++++++++    
    print >>out1file, '\t'.join(('#chrom','st','end','accession','mRNA_size','gene_strand','Frag_count','FPM','FPKM'))
     
    gene_finished = 0
    denominator = 0
    gene = []
    gene_list = []
     
    #calculate raw count for each gene
    for line in open(options.refgene_bed,'r'):
        frag_count_fr = 0
        mRNA_size = 0
        exon_ranges = Intersecter()
        if line.startswith(('#','track','browser')):
            continue
        fields = line.split()
        chrom = fields[0]
        exon_starts = int(fields[1])
        exon_ends = int(fields[2])
        gene_name = fields[3]
        gstrand = fields[5].replace(" ","_")
    
        mRNA_size += (exon_ends - exon_starts)
        exon_ranges.add_interval(Interval(exon_starts,exon_ends))
          
        # extract reads for mapped gene region
        try:
            alignedReads = obj.samfile.fetch(chrom,exon_starts,exon_ends)
        except:
            continue
        for aligned_read in alignedReads:
            flag = 0
            if aligned_read.is_qcfail:continue             # skip low quanlity
            if aligned_read.is_duplicate:continue          # skip duplicate read
            if aligned_read.is_secondary:continue          # skip non primary hit
            if options.skip_multi:
                if aligned_read.mapq < options.map_qual:
                    continue
        
            map_list = fetch_cigar_exon(chrom,aligned_read.pos,aligned_read.cigar)
            for map_coords in map_list:
                frag_st = map_coords[1]
                frag_end = map_coords[2]
                if len(exon_ranges.find(frag_st, frag_st +1 ))  > 0 or len(exon_ranges.find(frag_end, frag_end +1 )) > 0:
                    frag_count_fr += 1
        
        denominator = denominator + frag_count_fr  
        frag_count_str = str(frag_count_fr) 
        mRNA_size_str = str(mRNA_size)       
        gene = [chrom, exon_starts, exon_ends, gene_name, mRNA_size_str, gstrand, frag_count_str]
        gene_list.append(gene)
        gene_finished +=1
        print >>sys.stderr, "%d transcripts finished\r" % (gene_finished),
    
    FPM_fr = 0.0
    FPKM_fr = 0.0
    print >>sys.stderr, "\n"
    for gene in gene_list:
        chrom = gene[0]
        exon_starts = gene[1]
        exon_ends = gene[2]
        gene_name = gene[3]
        mRNA_size = float(gene[4])
        gstrand = gene[5]
        frag_count = float(gene[6])
        FPM_fr = frag_count * 1000000 / denominator
        FPKM_fr = frag_count * 1000000000 / (denominator * mRNA_size)
        print >>out1file, '\t'.join([str(i) for i in (chrom, exon_starts, exon_ends, gene_name, mRNA_size, gstrand, frag_count, FPM_fr,FPKM_fr)])
        out_stuff = "Reporting: " + gene_name + "                                    \r"
        print >>sys.stderr, out_stuff,
      
    out1file.close()
    
    if options.gene_based and options.junctions_bed and not options.skip_FPKMcount:
        if not (os.path.exists(options.gene_based) and os.path.exists(options.junctions_bed)):
            print >>sys.stderr, "cannot find gene or junctions bed file"
            sys.exit(0)
        else:
            go_gene(out1_file,options.gene_based) 
            out_file_tmp = out1_file.split(".")            
            out_file = out_file_tmp[0] + '.GENE_2.xls'
            skip_counter(options.do_dist,options.min_skips,out_file,options.junctions_bed)
            
    if options.junctions_bed and not (options.skip_FPKMcount or options.gene_based):
        if not (os.path.exists(options.junctions_bed)):
            print >>sys.stderr, "cannot find junctions bed file"
            sys.exit(0)
        else:    
            del gene_list[:]
            skip_counter(options.do_dist,options.min_skips,out1_file,options.junctions_bed)
    
    print >>sys.stderr, "All done !\n'"

if __name__ == '__main__':
     main()
