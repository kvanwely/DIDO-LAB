
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
    
def get_fasta (num_exons,bed_file,fasta_prefix,out_file):
    
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
    
    #open_junctions(m_file,test_junctions)
    
    get_rnd_junctions(num_exons,bed_file,test_junctions)
    
    outfile = open(out_file,'w')
    
    outcounter = 0
    print >>sys.stderr, "\n"                     
    print >>sys.stderr, "Upstream     "
    for junction in test_junctions:       
        
        try:
            test_chrom, test_left, test_right, test_ID, test_strand = junction.split('\t')
        except ValueError:
            continue
                
        if test_chrom != prev_chrom:
            open_fasta (test_chrom,fasta_prefix,fasta_data)
            prev_chrom = test_chrom
            del juncs_used[:]
            
        #test_chrom = junction[0] 
        test_leftcoord = int(test_left)
        test_rightcoord = int(test_right)
        test_name = test_ID
        test_strand = test_strand[0]   
        
        if test_strand == "-":
            prime_5_coord = test_rightcoord + 205
            prime_3_coord = test_rightcoord + 3
        if test_strand == "+":
            prime_3_coord = test_leftcoord - 205   
            prime_5_coord = test_leftcoord - 3
        
        list_data = [test_chrom,test_left,test_right]    
         
        try:
            junc_index = juncs_used.index(list_data)
        except ValueError:
            juncs_used.append(list_data)
            out_sequence_1 = get_sequence(prime_3_coord,prime_5_coord,test_strand,fasta_data)
            print >>sys.stderr, test_chrom + "  :  " + str(prime_3_coord) + "                      \r",
            #out_line = ">" + test_name + "_" + str(prime_3_coord)
            #print >>outfile, out_line
            #out_line = out_sequence_1
            #print >>outfile, out_line
            if len(out_sequence_1) > 100:
                out_series.append(out_sequence_1)
            
    seq_order = 0
    while seq_order < (len(out_sequence_1)-1) :
        num_a = 0
        num_c = 0
        num_g = 0
        num_t = 0
        for out_sequence_1 in out_series:
            if out_sequence_1[seq_order] == "A":
                num_a = num_a + 1
            if out_sequence_1[seq_order] == "C":
                num_c = num_c + 1
            if out_sequence_1[seq_order] == "G":
                num_g = num_g + 1                
            if out_sequence_1[seq_order] == "T":
                num_t = num_t + 1
            
            #out_sequence_1 = out_sequence_1[1:]
        seq_order = seq_order + 1    
        print >>sys.stderr,  str(seq_order) + "                      \r",
        out_line = str(seq_order) + '\t' + str(num_a) + '\t' + str(num_c) + '\t' + str(num_g) + '\t' + str(num_t)
        print >>outfile, out_line
    
    print >>sys.stderr, test_chrom + "  :  " + str(prime_3_coord) + "                        "    
    print >>sys.stderr, "Exon     "
    del out_series[:]
    
    for junction in test_junctions:       
        
        try:
            test_chrom, test_left, test_right, test_ID, test_strand = junction.split('\t')
        except ValueError:
            continue
            
        if test_chrom != prev_chrom:
            open_fasta (test_chrom,fasta_prefix,fasta_data)
            prev_chrom = test_chrom
            del juncs_used[:]
            
        test_leftcoord = int(test_left) 
        test_rightcoord = int(test_right) 
        test_strand = test_strand[0]  
          
        list_data = [test_chrom,test_left,test_right] 
                        
        test_leftcoord = int(test_left) + 3
        test_rightcoord = int(test_right) - 3    
         
        try:
            junc_index = juncs_used.index(list_data)
        except ValueError:
            juncs_used.append(list_data)
            out_sequence_1 = get_sequence(test_leftcoord,test_rightcoord,test_strand,fasta_data)
            seq_length = str(len(out_sequence_1))
            print >>sys.stderr, test_chrom + "  :  " + seq_length + "                      \r",
            #out_line = ">" + test_name + "_" + str(prime_3_coord)
            #print >>outfile, out_line
            #out_line = out_sequence_1
            #print >>outfile, out_line
            if len(out_sequence_1) > 30:
                out_series.append(out_sequence_1) 
    
    percent_a = 0.0
    percent_c = 0.0
    percent_g = 0.0
    percent_t = 0.0
    
    spacing_val = 0.0
    
    per_a_list = []
    per_c_list = []
    per_g_list = []
    per_t_list = []
    
    for count_1 in range(0,101):
         per_a_list.append(percent_a)
         per_c_list.append(percent_c)
         per_g_list.append(percent_g)
         per_t_list.append(percent_t)
         
    for out_sequence_2 in out_series:
        base = ""
        spacing_val = 101.0 / len(out_sequence_2)
        #print >>sys.stderr,  str(len(out_sequence_2)) + " : " + str(spacing_val)
        total_counter = 0.0
        count_2 = 0
    
        while count_2 < len(out_sequence_2)-1:
            spacing_counter = spacing_val
            base = out_sequence_2[count_2]
            while spacing_counter >= 1.0:
                count_3 = int(total_counter)
                if base == "A":
                    per_a_list[count_3] = per_a_list[count_3] + 1.0
                if base == "C":
                    per_c_list[count_3] = per_c_list[count_3] + 1.0     
                if base == "G":
                    per_g_list[count_3] = per_g_list[count_3] + 1.0  
                if base == "T":
                    per_t_list[count_3] = per_t_list[count_3] + 1.0 
                spacing_counter = spacing_counter - 1.0
                total_counter = total_counter + 1.0
            count_3 = int(total_counter)
            if base == "A":
                per_a_list[count_3] = per_a_list[count_3] + spacing_counter
            if base == "C":
                per_c_list[count_3] = per_c_list[count_3] + spacing_counter    
            if base == "G":
                per_g_list[count_3] = per_g_list[count_3] + spacing_counter  
            if base == "T":
                per_t_list[count_3] = per_t_list[count_3] + spacing_counter
            total_counter = total_counter + spacing_counter   
            count_2 = count_2 + 1
                                                               
    for count_5 in range(0,99):
        out_line = str(count_5+1) + '\t' + str(per_a_list[count_5]) + '\t' + str(per_c_list[count_5]) + '\t' + str(per_g_list[count_5]) + '\t' + str(per_t_list[count_5]) 
        out_line = out_line.replace(".",",")
        print >>outfile, out_line
    
    print >>sys.stderr, test_chrom + "  :  " + seq_length + "                        "                
    print >>sys.stderr, "Downstream  "
    del out_series[:]

    for junction in test_junctions:       
        
        try:
            test_chrom, test_left, test_right, test_ID, test_strand = junction.split('\t')
        except ValueError:
            continue
                
        if test_chrom != prev_chrom:
            open_fasta (test_chrom,fasta_prefix,fasta_data)
            prev_chrom = test_chrom
            del juncs_used[:]
            
        #test_chrom = junction[0] 
        test_leftcoord = int(test_left)
        test_rightcoord = int(test_right)
        test_name = test_ID
        test_strand = test_strand[0]   
        
        if test_strand == "+":
            prime_5_coord = test_rightcoord + 205
            prime_3_coord = test_rightcoord + 3
        if test_strand == "-":
            prime_3_coord = test_leftcoord - 205   
            prime_5_coord = test_leftcoord - 3
        
        list_data = [test_chrom,test_left,test_right]    
         
        try:
            junc_index = juncs_used.index(list_data)
        except ValueError:
            juncs_used.append(list_data)
            out_sequence_1 = get_sequence(prime_3_coord,prime_5_coord,test_strand,fasta_data)
            print >>sys.stderr, test_chrom + "  :  " + str(prime_3_coord) + "                      \r",
            #out_line = ">" + test_name + "_" + str(prime_3_coord)
            #print >>outfile, out_line
            #out_line = out_sequence_1
            #print >>outfile, out_line
            if len(out_sequence_1) > 100:
                out_series.append(out_sequence_1)
            
    seq_order = 0
    while seq_order < (len(out_sequence_1)-1) :
        num_a = 0
        num_c = 0
        num_g = 0
        num_t = 0
        for out_sequence_1 in out_series:
            if out_sequence_1[seq_order] == "A":
                num_a = num_a + 1
            if out_sequence_1[seq_order] == "C":
                num_c = num_c + 1
            if out_sequence_1[seq_order] == "G":
                num_g = num_g + 1                
            if out_sequence_1[seq_order] == "T":
                num_t = num_t + 1
            
            #out_sequence_1 = out_sequence_1[1:]
        seq_order = seq_order + 1    
        print >>sys.stderr,  str(seq_order) + "                      \r",
        out_line = str(seq_order) + '\t' + str(num_a) + '\t' + str(num_c) + '\t' + str(num_g) + '\t' + str(num_t)
        print >>outfile, out_line
    
    print >>sys.stderr, test_chrom + "  :  " + str(prime_3_coord) + "                        "    
    print >>sys.stderr, "Done  "

def main():
    usage = "%prog [options]" + "\n"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-e","--number_exons",action="store",type="int",dest="num_exons",help="number of exons to use")
    parser.add_option("-b","--bed-file",action="store",type="string",dest="bed_file",help="Exon map as bedfile")  
    parser.add_option("-f","--fasta-file",action="store",type="string",dest="fasta_prefix",help="fasta prefix")
    parser.add_option("-o","--out-file",action="store",type="string",dest="out_file",help="output file")   
    
    (options,args)=parser.parse_args()

    if not (options.num_exons and options.bed_file and options.fasta_prefix and options.out_file) :
        parser.print_help()
        sys.exit(0)
    else:
        get_fasta (options.num_exons,options.bed_file,options.fasta_prefix,options.out_file)

if __name__ == '__main__':
     main()
    
        