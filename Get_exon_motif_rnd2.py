
#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
     print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " RSeQC needs python2.7!\n"
     sys.exit()

import re
from optparse import OptionParser
import warnings
import collections
from collections import defaultdict as defdic
import math
try:
    from string import maketrans as mkts
except ImportError:
    mkts = str.maketrans
import sets
import string
import random

__license__ = "GPL"
__version__="2.6"
    
# get junctions from tophat2 junctions.bed file and put values in same order as exons. 

def open_fasta (do_chrom,fasta_prefix,fasta_data):
     
    del fasta_data[:]
       
    fasta2_file = fasta_prefix + ".fa." + do_chrom    
    fastafile = open(fasta2_file,"r")
    for line in fastafile:
        # line = line.rstrip('\n')
                
        if line.startswith('>') or line.startswith('chr'):
            continue
        else:
            line = line.rstrip('\n')
            fasta_data.append(line)
            
    fastafile.close ()
    return fasta_data
    
def motif_maker(motif_length,motif_list):

    bases = {'A', 'C', 'G', 'T'}
    
    tmp_motifs_2 = []
    
    del motif_list[:]
    
    for a_base in bases:
        motif_list.append(a_base)
    
    while (len(motif_list[-1])) < motif_length:
         
        for motif_1 in motif_list:
            for a_base in bases:
                motif_2 = motif_1 + a_base
                tmp_motifs_2.append(motif_2)
                
        del motif_list[:] 
        for motif_2 in tmp_motifs_2:       
            motif_list.append(motif_2)
        del tmp_motifs_2[:]        
        
    return motif_list
        
def get_sequence(left_coord,right_coord,strand,fasta_data):
    
    all_sequence = "".join(fasta_data)
    
    fasta_sequence = all_sequence[left_coord:right_coord].upper()
    
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    
    if strand == "-":    
        fasta_sequence = "".join(complement.get(base,base) for base in reversed(fasta_sequence))
        
    return(fasta_sequence)
            
def get_rnd_junctions (num_exons,bed_file,test_junctions):
      
    del test_junctions[:]   

    temp_junctions_1 = []
    
    selected_junctions = []
    
    infile = open(bed_file,"r")
    for line in infile:

        line = line.rstrip('\n')
        
        if line.startswith("chr"):
            try:
                junc_data = line.split('\t')
            except ValueError:
                continue    
            
            curr_chrom = junc_data[0]
            if "ran" in curr_chrom or len(curr_chrom) > 5:
                continue
            curr_left = junc_data[1]
            curr_right = junc_data[2]
            curr_strand = junc_data[5]

            curr_ID = junc_data[3] + "_" + curr_left + "_" + curr_right
            
            linedata_1 = curr_chrom + '@@' + curr_left + '@@' + curr_right + '@@' + curr_ID + '@@' + curr_strand 
                            
            temp_junctions_1.append(linedata_1)
            
    temp_junctions_1.sort()
    
    num_selected = 0
    
    while num_selected < num_exons:
        
        select_junc = random.choice(temp_junctions_1)
        
        try:
            select_index = selected_junctions.index(select_junc)
        except ValueError:
            linedata_2 = select_junc
            test_junctions.append(linedata_2)    
            selected_junctions.append(select_junc)
            num_selected += 1            
            print >>sys.stderr, str(num_selected) + "      \r",
        continue
    
    selected_junctions.sort()
    
    for linedata_3 in selected_junctions:
        linedata_2 = linedata_3.replace("@@","\t") 
        test_junctions.append(linedata_2)
    
    del selected_junctions[:]
    del temp_junctions_1[:]
            
    infile.close ()
    return test_junctions
    
def get_fasta (num_exons,bed_file,fasta_prefix,out_file,numbases):
    
    list_data = []
    junction = []
    juncs_used = []
    fasta_data = []
    test_junctions = []
    test_exons = []
    already_tested = []
    prev_chrom = 'abc' 
    transcript_list = []
    out_series = []
            
    motif_list = []
    motif_maker (numbases,motif_list)    
    
    #open_junctions(m_file,test_junctions)
    
    get_rnd_junctions(num_exons,bed_file,test_junctions)
    
    outfile = open(out_file,'w')
    
    outcounter = 0
    
    for junction in test_junctions:       
        
        try:
            test_chrom, test_left, test_right, test_ID, test_strand = junction.split('\t')
        except ValueError:
            continue
                
        test_leftcoord = int(test_left)
        test_rightcoord = int(test_right)
        test_name = test_ID
        test_strand = test_strand[0]
                
        if test_chrom != prev_chrom:
            open_fasta (test_chrom,fasta_prefix,fasta_data)
            prev_chrom = test_chrom
            del juncs_used[:]

        #selection of 5prime or 3prime ends of exons
        prime_5_coord = test_leftcoord + 3
        prime_3_coord = test_rightcoord - 3

        list_id = str(prime_3_coord)   
        list_data = [test_chrom,list_id]   
         
        try:
            junc_index = juncs_used.index(list_data)
        except ValueError:
            juncs_used.append(list_data)
            out_sequence_1 = get_sequence (prime_5_coord,prime_3_coord,test_strand,fasta_data)
            print >>sys.stderr, test_chrom + "  :  " + list_id + "                      \r",
            #out_line = ">" + test_name + "_" + str(prime_3_coord)
            #print >>outfile, out_line
            #out_line = out_sequence_1
            #print >>outfile, out_line
            if len(out_sequence_1) > 100:
                out_series.append(out_sequence_1)
            
    print >>sys.stderr, "\n"  
    
    for motif_2 in motif_list:
        motif_freq = 0.0
        for out_sequence_1 in out_series:
            motif_num = out_sequence_1.count(motif_2)
            motif_freq = motif_freq + (float(motif_num) / (float(len(out_series) * float(len(out_sequence_1) ))))
            
        out_line = motif_2 + '\t' + str(motif_freq) 
        print >>sys.stderr,  out_line + "                         \r",
        print >>outfile, out_line 

    print >>sys.stderr,  out_line           
    print >>sys.stderr, "Done  "

def main():
    usage = "%prog [options]" + "\n"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-e","--number_exons",action="store",type="int",dest="num_exons",help="number of exons to use")
    parser.add_option("-b","--bed-file",action="store",type="string",dest="bed_file",help="Exon map as bedfile")  
    parser.add_option("-f","--fasta-file",action="store",type="string",dest="fasta_prefix",help="fasta prefix")
    parser.add_option("-o","--out-file",action="store",type="string",dest="out_file",help="output file")   
    parser.add_option("-n","--num_bases",action="store",type="int",dest="numbases",help="dept cutoff") 
    
    (options,args)=parser.parse_args()

    if not (options.num_exons and options.bed_file and options.fasta_prefix and options.out_file and options.numbases) :
        parser.print_help()
        sys.exit(0)
    else:
        get_fasta (options.num_exons,options.bed_file,options.fasta_prefix,options.out_file,options.numbases)

if __name__ == '__main__':
     main()
    
        