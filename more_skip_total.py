
#this script will take the output from the skipping and get general base composition for the exons skipped significantly more

#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
     print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " RSeQC needs python2.7!\n"
     sys.exit()

import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
import sets

__license__ = "GPL"
__version__="2.6"

# get junctions from tophat2 junctions.bed file and put values in same order as exons. 

def test_for_tss (test_leftcoord,test_rightcoord,test_strand,transcript_list):
    
    is_tss = -1
    
    for transcript in transcript_list:
        try:
            test3_chrom, test3_left, test3_right, test3_ID, dat3_0, test3_strand = transcript.split("\t")            
        except ValueError:
            continue   
            
        test3_leftcoord = int(test3_left)
        test3_rightcoord = int(test3_right)
        
        if ((test_leftcoord - test3_leftcoord > -10) and (test_leftcoord - test3_leftcoord < 10)):
            #if (test_strand == "+"):
            is_tss = 1
            
        if ((test_rightcoord - test3_rightcoord > -10) and (test_rightcoord - test3_rightcoord < 10)):
            #if (test_strand == "-"):
            is_tss = 1    
        
    return (is_tss)

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
        
def open_junctions (skip_file,test_junctions):
     
    del test_junctions[:]
        
    infile = open(skip_file,"r")
    for line in infile:
        
        line = line.rstrip('\n')
        
        if not line.startswith("ID"):
            test_junctions.append(line)
        
    infile.close ()
    return test_junctions
    
def open_exp_list(chrom,gene_exp):
    
    gene_exp_list = []
    
    gene_ex_file = open(gene_exp,'r')
    for line in gene_ex_file:
        try:
            if line.startswith(('#','track','browser')):
                continue
            line = line.rstrip('\n')
            line_dat = line.split('\t')
            if line_dat[0].upper() == chrom.upper():
                gene_exp_list.append(line)
        except:
            print >>sys.stderr," skipped this line: " + line,
            continue
    
    gene_ex_file.close()
    return(gene_exp_list)
    
def get_fasta (cutoff,skip_file,fasta_prefix,transcripts,out_file):
    
    list_data = []
    junction = []
    juncs_used = []
    fasta_data = []
    test_junctions = []
    out_series = []
    prev_chrom = 'abc' 
    out_sequence_1 = ""
    
    open_junctions(skip_file,test_junctions)
    
    outfile = open(out_file,'w')
    
    print >>sys.stderr, "\n"                     
    print >>sys.stderr, "Upstream     "

    for junction in test_junctions:     
        
        more_less_skip = 0.0
        test_dept_val = 1.0  
        
        try:
            test_chrom, test_left, test_right, test_ID, test_strand, dat_1, dat_2, dat_3, dat_4, dat_5, dat_6, dat_7, dat_8, dat_0, val_FDR, level_diff = junction.split('\t')
        except ValueError:
            continue
        
        test_dept_val = float(val_FDR)
        more_less_skip = float(level_diff)
        
        if (test_dept_val > cutoff) or (more_less_skip < 0):
            continue
        
        if test_chrom != prev_chrom:
            open_fasta (test_chrom,fasta_prefix,fasta_data)
            transcript_list = open_exp_list(test_chrom,transcripts)
            prev_chrom = test_chrom
            del juncs_used[:]
                        
        #test_chrom = junction[0] 
        test_leftcoord = int(test_left)
        test_rightcoord = int(test_right)
        test_name = test_ID
        test_strand = test_strand[0]  
        
        is_tss = test_for_tss(test_leftcoord,test_rightcoord,test_strand,transcript_list)    
        if is_tss < 0:
            continue
        
        #selection of 5prime or 3prime ends of exons
        if test_strand == "-":
            prime_5_coord = test_rightcoord + 205
            prime_3_coord = test_rightcoord + 3
        if test_strand == "+":
            prime_3_coord = test_leftcoord - 205   
            prime_5_coord = test_leftcoord - 3
        
        list_id = str(prime_3_coord)   
        list_data = [test_chrom,list_id]   
         
        try:
            junc_index = juncs_used.index(list_data)
        except ValueError:
            juncs_used.append(list_data)
            out_sequence_1 = get_sequence (prime_3_coord,prime_5_coord,test_strand,fasta_data)
            print >>sys.stderr, test_chrom + "  :  " + list_id + "                      \r",
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
    
    print >>sys.stderr, test_chrom + "  :  " + list_id + "                        "    
    print >>sys.stderr, "Exon     "
    del out_series[:]
        
    for junction in test_junctions:       
        
        more_less_skip = 0.0
        test_dept_val = 1.0  
        
        try:
            test_chrom, test_left, test_right, test_ID, test_strand, dat_1, dat_2, dat_3, dat_4, dat_5, dat_6, dat_7, dat_8, dat_0, val_FDR, level_diff = junction.split('\t')
        except ValueError:
            continue
        
        test_dept_val = float(val_FDR)
        more_less_skip = float(level_diff)
        
        if (test_dept_val > cutoff) or (more_less_skip < 0):
            continue
        
        if test_chrom != prev_chrom:
            open_fasta (test_chrom,fasta_prefix,fasta_data)
            transcript_list = open_exp_list(test_chrom,transcripts)    
            prev_chrom = test_chrom
            del juncs_used[:]
                        
        test_leftcoord = int(test_left) 
        test_rightcoord = int(test_right) 
        test_strand = test_strand[0]  
          
        list_data = [test_chrom,test_left,test_right]   
        
        is_tss = test_for_tss(test_leftcoord,test_rightcoord,test_strand,transcript_list)    
        if is_tss < 0:
            continue
            
        test_leftcoord = int(test_left) + 3
        test_rightcoord = int(test_right) - 3    
         
        try:
            junc_index = juncs_used.index(list_data)
        except ValueError:
            juncs_used.append(list_data)
            out_sequence_1 = get_sequence (test_leftcoord,test_rightcoord,test_strand,fasta_data)
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
        
        more_less_skip = 0.0
        test_dept_val = 1.0  
        
        try:
            test_chrom, test_left, test_right, test_ID, test_strand, dat_1, dat_2, dat_3, dat_4, dat_5, dat_6, dat_7, dat_8, dat_0, val_FDR, level_diff = junction.split('\t')
        except ValueError:
            continue
        
        test_dept_val = float(val_FDR)
        more_less_skip = float(level_diff)
        
        if (test_dept_val > cutoff) or (more_less_skip < 0):
            continue
        
        if test_chrom != prev_chrom:
            open_fasta (test_chrom,fasta_prefix,fasta_data)
            transcript_list = open_exp_list(test_chrom,transcripts)    
            prev_chrom = test_chrom
            del juncs_used[:]
                        
        #test_chrom = junction[0] 
        test_leftcoord = int(test_left)
        test_rightcoord = int(test_right)
        test_name = test_ID
        test_strand = test_strand[0]  
        
        is_tss = test_for_tss(test_leftcoord,test_rightcoord,test_strand,transcript_list)    
        if is_tss < 0:
            continue
            
        #selection of 5prime or 3prime ends of exons
        if test_strand == "-":
            prime_5_coord = test_leftcoord - 3
            prime_3_coord = test_leftcoord - 205
        if test_strand == "+":
            prime_3_coord = test_rightcoord + 3   
            prime_5_coord = test_rightcoord + 205
        
        list_id = str(prime_3_coord)   
        list_data = [test_chrom,list_id]   
         
        try:
            junc_index = juncs_used.index(list_data)
        except ValueError:
            juncs_used.append(list_data)
            out_sequence_1 = get_sequence (prime_3_coord,prime_5_coord,test_strand,fasta_data)
            print >>sys.stderr, test_chrom + "  :  " + list_id + "                      \r",
            #out_line = ">" + test_name + "_" + str(prime_3_coord)
            #print >>outfile, out_line
            #out_line = out_sequence_1
            #print >>outfile, out_line
            if len(out_sequence_1) > 100:
                out_series.append(out_sequence_1)
            
            
    print >>sys.stderr, "\n"  
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

    print >>sys.stderr, test_chrom + "  :  " + list_id + "                      "                   
    print >>sys.stderr, "Done  "

def main():
    usage = "%prog [options]" + "\n"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-p","--p-value",action="store",type="float",dest="cutoff",help="dept cutoff")    
    parser.add_option("-m","--mats-file",action="store",type="string",dest="skip_file",help="mats results file")
    parser.add_option("-f","--fasta-file",action="store",type="string",dest="fasta_prefix",help="fasta prefix")
    parser.add_option("-o","--out-file",action="store",type="string",dest="out_file",help="output file")   
    parser.add_option("-t","--trans",action="store",type="string",dest="transcripts",help="transcript bed file") 
    
    (options,args)=parser.parse_args()

    if not (options.cutoff and options.skip_file and options.fasta_prefix and options.transcripts and options.out_file):
        parser.print_help()
        sys.exit(0)
    else:
        get_fasta (options.cutoff,options.skip_file,options.fasta_prefix,options.transcripts,options.out_file)

if __name__ == '__main__':
     main()
    
        