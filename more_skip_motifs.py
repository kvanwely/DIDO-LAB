
#this script will take the output from the skipping and get motifs for the exons skipped significantly more

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
    
def get_fasta (cutoff,numbases,skip_file,fasta_prefix,transcripts,out_file):
    
    list_data = []
    junction = []
    juncs_used = []
    fasta_data = []
    test_junctions = []
    out_series = []
    prev_chrom = 'abc' 
    out_sequence_1 = ""
    
    motif_list = []
    motif_maker (numbases,motif_list)
    
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
                
    for motif_2 in motif_list:
        motif_freq = 0.0
        for out_sequence_1 in out_series:
            motif_num = out_sequence_1.count(motif_2)
            motif_freq = motif_freq + (float(motif_num) / (float(len(out_series) * float(len(out_sequence_1) ))))
            
        out_line = motif_2 + '\t' + str(motif_freq) 
        print >>sys.stderr,  out_line + "                         \r",
        print >>outfile, out_line
    
    print >>outfile, "@@@\n"
    print >>sys.stderr,  out_line 
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
                        
        test_leftcoord = int(test_left) + 3
        test_rightcoord = int(test_right) - 3
        test_strand = test_strand[0]  
        
        is_tss = test_for_tss(test_leftcoord,test_rightcoord,test_strand,transcript_list)    
        if is_tss < 0:
            continue
          
        list_data = [test_chrom,test_left,test_right]   
         
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

    for motif_2 in motif_list:
        motif_freq = 0.0
        for out_sequence_1 in out_series:
            motif_num = out_sequence_1.count(motif_2)
            motif_freq = motif_freq + (float(motif_num) / (float(len(out_series) * float(len(out_sequence_1) ))))
            
        out_line = motif_2 + '\t' + str(motif_freq) 
        print >>sys.stderr,  out_line + "                         \r",
        print >>outfile, out_line
    
    print >>outfile, "@@@\n"
    print >>sys.stderr,  out_line       
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
    parser.add_option("-p","--p-value",action="store",type="float",dest="cutoff",help="dept cutoff")    
    parser.add_option("-n","--num_bases",action="store",type="int",dest="numbases",help="dept cutoff")  
    parser.add_option("-m","--mats-file",action="store",type="string",dest="skip_file",help="mats results file")
    parser.add_option("-f","--fasta-file",action="store",type="string",dest="fasta_prefix",help="fasta prefix")
    parser.add_option("-o","--out-file",action="store",type="string",dest="out_file",help="output file") 
    parser.add_option("-t","--trans",action="store",type="string",dest="transcripts",help="transcript bed file")   
    
    (options,args)=parser.parse_args()

    if not (options.numbases and options.cutoff and options.skip_file and options.fasta_prefix and options.transcripts and options.out_file):
        parser.print_help()
        sys.exit(0)
    else:
        get_fasta (options.cutoff,options.numbases,options.skip_file,options.fasta_prefix,options.transcripts,options.out_file)

if __name__ == '__main__':
     main()
    
        