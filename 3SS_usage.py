
#this script will take a jucntion comparison file as input, find singificntly altered 3SS, and calculate their parameters
#three control sets and three test sets are expected in columns
#this script uses the Yeo & Burge data vor 3'SS strength calculation, adapted from https://github.com/kepbod/maxentpy
#the associated matrix files are required and can be downloaded from https://github.com/kepbod/maxentpy/tree/master/maxentpy/data
#an exon model in BED format is required to test for valid junctions
#genome data in fasta format are required (using chromosome names as extension) to get splice site sequences

#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
     print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + "this script needs python2.7!\n"
     sys.exit()

import re
import string
from optparse import OptionParser
import warnings
import string
from collections import defaultdict
import collections
import math
try:
    from string import maketrans
except ImportError:
    maketrans = str.maketrans
import sets
from scipy import stats
import numpy as np

__license__ = "GPL"
__version__="2.6"

bgd_5 = {'A': 0.27, 'C': 0.23, 'G': 0.23, 'T': 0.27}
cons1_5 = {'A': 0.004, 'C': 0.0032, 'G': 0.9896, 'T': 0.0032}
cons2_5 = {'A': 0.0034, 'C': 0.0039, 'G': 0.0042, 'T': 0.9884}

bgd_3 = {'A': 0.27, 'C': 0.23, 'G': 0.23, 'T': 0.27}
cons1_3 = {'A': 0.9903, 'C': 0.0032, 'G': 0.0034, 'T': 0.0030}
cons2_3 = {'A': 0.0027, 'C': 0.0037, 'G': 0.9905, 'T': 0.0030}
            
def load_matrix3(matrix_file):
    #matrix_f = dir_path + '/data/score3_matrix.txt'
    matrix = defaultdict(dict)
    with open(matrix_file, 'r') as f:
        for line in f:
            n, m, s = line.split()
            matrix[int(n)][int(m)] = float(s)
    return matrix
    
def hashseq(fa):
    table = maketrans('ACGT', '0123')
    seq = fa.translate(table)
    return sum(int(j) * 4**(len(seq) - i - 1) for i, j in enumerate(seq))
    
def score3(fa, matrix=None):

    t_over_c = 0.0
    
    # check length of fa
    if len(fa) != 23:
        sys.exit('Wrong length of fa!')
    # check matrix
    if not matrix:
        matrix = load_matrix3()
    # for key elements
    key = fa[18:20].upper()
    score = cons1_3[key[0]] * cons2_3[key[1]] / (bgd_3[key[0]] * bgd_3[key[1]])
    # for rest elements
    rest = (fa[:18] + fa[20:]).upper()
    rest_score = 1
    rest_score *= matrix[0][hashseq(rest[:7])]
    rest_score *= matrix[1][hashseq(rest[7:14])]
    rest_score *= matrix[2][hashseq(rest[14:])]
    rest_score *= matrix[3][hashseq(rest[4:11])]
    rest_score *= matrix[4][hashseq(rest[11:18])]
    rest_score /= matrix[5][hashseq(rest[4:7])]
    rest_score /= matrix[6][hashseq(rest[7:11])]
    rest_score /= matrix[7][hashseq(rest[11:14])]
    rest_score /= matrix[8][hashseq(rest[14:18])]
    
    rest_2 = fa[:18].upper()
    t_over_c =  4 * ( float(rest_2.count("T")) +1 ) / ( float(rest_2.count("C")) +1)
    # final score
    return math.log(score * rest_score, 2)

def get_sequence(left_coord,right_coord,strand,fasta_data):
    
    all_sequence = "".join(fasta_data)
    all_sequence = all_sequence.replace(" ","")
    
    fasta_sequence = all_sequence[left_coord:right_coord].upper()
    
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    
    if strand == "-":    
        fasta_sequence = "".join(complement.get(base,base) for base in reversed(fasta_sequence))
     
    #print >>sys.stderr, fasta_sequence    
    return(fasta_sequence)
    
def open_fasta (do_chrom,fasta_prefix,fasta_data):
     
    del fasta_data[:]
       
    fasta2_file = fasta_prefix + ".fa." + do_chrom    
    fastafile = open(fasta2_file,"r")
    for line in fastafile:
        # line = line.rstrip('\n')
                
        if line.startswith('>') or line.startswith('chr'):
            continue
        else:
            line = line.rstrip("\n")
            fasta_data.append(line)
            
    fastafile.close ()
    return fasta_data
    
def open_exon_data (do_chrom,r_file,ref_exons):
     
    del ref_exons[:]
    reffile = open(r_file)
    for line in reffile: 
         
        line = line.rstrip('\n')
        
        tmp_dat = line.split("\t")
        tmp_chrom = tmp_dat[0]
        if line.startswith("chr") and tmp_chrom == do_chrom:
            ref_exons.append(line)
        
    reffile.close ()
    return ref_exons
                               
def open_junctions (m_file,t_chrom):
    
    t_junctions = [] 
    del t_junctions[:]
        
    infile = open(m_file,"r")
    for line in infile:
        
        line = line.rstrip('\n')
        line = line.replace(",",".")
        
        if "DIV" in line:
            continue
            
        if "NA"  in line:
            continue
               
        if line.startswith("chr") and t_chrom == "all":
            t_junctions.append(line)
        else:    
            tmp_dat = line.split("\t")
            tmp_chrom = tmp_dat[0]
            if line.startswith("chr") and tmp_chrom == t_chrom:
                t_junctions.append(line)
               
    infile.close ()
    return t_junctions
        
def get_3SS (cutoff,m_file,r_file,x_file,fasta_prefix,out_prefix):
    
    junction = []
    test_junctions = []    
    super_junctions = []
    fasta_data = []
    
    prev_chrom = 'abc' 
    
    used_3SS = []
    used_5SS = []
    
    ref_exons = []
    exon_line = []
    
    testable_junctions = []
    temp_junctions = []
    counter_junctions = []
                           
    test_junctions = open_junctions(m_file,"all")   
    print >>sys.stderr, "Pass 1: \n" 
    # pass 1 serves to select significantly changed 3prime sites
    
    for tes_junction in test_junctions:      
        
        del temp_junctions[:]
        SS3_mapped = 0

        try:
            test_data = tes_junction.split('\t')
        except ValueError:
            continue
                
        test_chrom = test_data[0]
        test_left = int(test_data[2])
        test_right = int(test_data[3])
        test_strand = test_data[1]
                
        if test_chrom != prev_chrom:
            #open_fasta (test_chrom,fasta_prefix,fasta_data)
            del ref_exons[:]
            ref_exons = open_exon_data(test_chrom,r_file,ref_exons)  
            del super_junctions [:]
            super_junctions = open_junctions(m_file,test_chrom)
            prev_chrom = test_chrom   
            
        if test_strand == "+":
            test_3SS = test_right
        if test_strand == "-":
            test_3SS = test_left
                                   
        for ref_exon in ref_exons:
            exon_line = ref_exon.split("\t")
            exon_strand = exon_line[5]
            exon_strand = exon_strand[0]
            exon_left = int(exon_line[1])
            exon_right = int(exon_line[2])
            
            if exon_strand == "+":
                exon_3SS_1 = exon_left - 2
                exon_3SS_2 = exon_left + 2
            if exon_strand == "-":
                exon_3SS_1 = exon_right - 2
                exon_3SS_2 = exon_right + 2
            
            if exon_strand == test_strand and test_3SS > exon_3SS_1 and test_3SS < exon_3SS_2 :
                SS3_mapped += 1
                                
        if SS3_mapped < 1:
            continue
                        
        # get an indication of usage of this junction 
        usage_WT_1 = float(test_data[4]) 
        usage_WT_2 = float(test_data[5]) 
        usage_WT_3 = float(test_data[6])  
        
        usage_mut_1 = float(test_data[7])
        usage_mut_2 = float(test_data[8]) 
        usage_mut_3 = float(test_data[9])  
                
        if test_strand == "+":
            prime_5_coord = test_right - 20
            prime_3_coord = test_right + 3 
        if test_strand == "-":
            prime_5_coord = test_left - 3
            prime_3_coord = test_left + 20
        
        sup_usage_WT_1 = 0.0
        sup_usage_WT_2 = 0.0
        sup_usage_WT_3 = 0.0
        sup_usage_mut_1 = 0.0
        sup_usage_mut_2 = 0.0
        sup_usage_mut_3 = 0.0
                        
        sup_SS3_mapped = 0
                                           
        for sup_junction in super_junctions :
                        
            try:
                sup_data = sup_junction.split('\t')
            except ValueError:
                continue
            
            sup_strand = sup_data[1]
                                        
            if test_strand != sup_strand:
                continue
            
            sup_left = int(sup_data[2])
            sup_right = int(sup_data[3])
                                            
            if test_strand == "+" and sup_strand == "+" :
                if test_left - sup_left > -3 and test_left - sup_left < 3 : 
                    if not ( test_right - sup_right > -3 and test_right - sup_right < 3 ):
                                               
                        sup_3SS = sup_right
                        
                        for ref_exon in ref_exons:
                            exon_line = ref_exon.split("\t")
                            exon_strand = exon_line[5]
                            exon_strand = exon_strand[0]
                            exon_left = int(exon_line[1])
                            exon_right = int(exon_line[2])
            
                            if exon_strand == "+":
                                exon_3SS_1 = exon_left - 2
                                exon_3SS_2 = exon_left + 2
                            if exon_strand == "-":
                                exon_3SS_1 = exon_right - 2
                                exon_3SS_2 = exon_right + 2
            
                            if exon_strand == sup_strand and sup_3SS > exon_3SS_1 and sup_3SS < exon_3SS_2 :
                                sup_SS3_mapped += 1
                                                                     
                                temp_junctions.append(sup_junction)
                                                            
                                sup_usage_WT_1 = sup_usage_WT_1 + float(sup_data[4]) 
                                sup_usage_WT_2 = sup_usage_WT_2 + float(sup_data[5])
                                sup_usage_WT_3 = sup_usage_WT_3 + float(sup_data[6])
                                
                                sup_usage_mut_1 =  sup_usage_mut_1 + float(sup_data[7])
                                sup_usage_mut_2 =  sup_usage_mut_2 + float(sup_data[8])
                                sup_usage_mut_3 =  sup_usage_mut_3 + float(sup_data[9])
                                                                                                           
            if test_strand == "-" and sup_strand == "-" :
                if test_right - sup_right > -3 and test_right - sup_right < 3 :
                    if not ( test_left - sup_left > -3 and test_left - sup_left < 3 ):
                                                                    
                        sup_3SS = sup_left
                                        
                        for ref_exon in ref_exons:
                            exon_line = ref_exon.split("\t")
                            exon_strand = exon_line[5]
                            exon_strand = exon_strand[0]
                            exon_left = int(exon_line[1])
                            exon_right = int(exon_line[2])
            
                            if exon_strand == "+":
                                exon_3SS_1 = exon_left - 2
                                exon_3SS_2 = exon_left + 2
                            if exon_strand == "-":
                                exon_3SS_1 = exon_right - 2
                                exon_3SS_2 = exon_right + 2
            
                            if exon_strand == sup_strand and sup_3SS > exon_3SS_1 and sup_3SS < exon_3SS_2 :
                                sup_SS3_mapped += 1
                                                 
                                temp_junctions.append(sup_junction)
                                                                                                
                                sup_usage_WT_1 = sup_usage_WT_1 + float(sup_data[4]) 
                                sup_usage_WT_2 = sup_usage_WT_2 + float(sup_data[5])
                                sup_usage_WT_3 = sup_usage_WT_3 + float(sup_data[6])
                                
                                sup_usage_mut_1 =  sup_usage_mut_1 + float(sup_data[7])
                                sup_usage_mut_2 =  sup_usage_mut_2 + float(sup_data[8])
                                sup_usage_mut_3 =  sup_usage_mut_3 + float(sup_data[9])
                                          
        if sup_SS3_mapped > 0: 
            
            if sup_usage_WT_1 > 0.0 :
                ratio_WT_1 = usage_WT_1 / sup_usage_WT_1
            else:
                sup_usage_WT_1 = 0.00001
                ratio_WT_1 = usage_WT_1 / sup_usage_WT_1
                
            if sup_usage_WT_2 > 0.0 :
                ratio_WT_2 = usage_WT_2 / sup_usage_WT_2
            else:
                sup_usage_WT_2 = 0.00001
                ratio_WT_2 = usage_WT_2 / sup_usage_WT_2
                
            if sup_usage_WT_3 > 0.0 :
                ratio_WT_3 = usage_WT_3 / sup_usage_WT_3
            else:
                sup_usage_WT_3 = 0.00001
                ratio_WT_3 = usage_WT_3 / sup_usage_WT_3
        
            if sup_usage_mut_1 > 0.0 :
                ratio_mut_1 = usage_mut_1 / sup_usage_mut_1
            else:
                sup_usage_mut_1 = 0.00001
                ratio_mut_1 = usage_mut_1 / sup_usage_mut_1
                
            if sup_usage_mut_2 > 0.0 :
                ratio_mut_2 = usage_mut_2 / sup_usage_mut_2
            else:
                sup_usage_mut_2 = 0.00001
                ratio_mut_2 = usage_mut_2 / sup_usage_mut_2
                
            if sup_usage_mut_3 > 0.0 :
                ratio_mut_3 = usage_mut_3 / sup_usage_mut_3
            else:
                sup_usage_mut_3 = 0.00001
                ratio_mut_3 = usage_mut_3 / sup_usage_mut_3
               
            usage_p_value = 1.1
              
            if ratio_WT_1 > -1.0 and ratio_WT_2 > -1.0 and ratio_WT_3 > -1.0:
                if ratio_mut_1 > -1.0 and ratio_mut_2 > -1.0 and ratio_mut_3 > -1.0:
                                        
                    usage_fin_WT = ([ratio_WT_1,ratio_WT_2,ratio_WT_3])
                    usage_fin_mut = ([ratio_mut_1,ratio_mut_2,ratio_mut_3])
                    usage_test,usage_p_value = stats.ttest_ind(usage_fin_WT,usage_fin_mut,equal_var=False)
                    
                    print >>sys.stderr, test_chrom + "  :  " + str(test_left) + "   \t" + str(len(testable_junctions)) +  "   \t" + str(len(counter_junctions)) + "               \r",   
                                         
            if usage_p_value > cutoff:
                continue
            else:
                testable_junctions.append(tes_junction)  
                for tmp_junction in temp_junctions:
                    try:
                        counter_index = counter_junctions.index(tmp_junction)
                    except ValueError:                   
                        counter_junctions.append(tmp_junction)
                del temp_junctions[:]   
            
    #compare selected junctions to competiton (pairwise)
    print >>sys.stderr, "\n"  
    print >>sys.stderr, "Pass 2: \n"   
    
    ss_matrix = load_matrix3(x_file)  
    
    outfile_1 = open(out_prefix + ".all.txt",'w') 
                 
    for tes_junction in testable_junctions:      
        
        try:
            test_data = tes_junction.split('\t')
        except ValueError:
            continue
            
        test_chrom = test_data[0]
        test_left = int(test_data[2])
        test_right = int(test_data[3])
        test_strand = test_data[1]
                
        if test_chrom != prev_chrom:
            open_fasta (test_chrom,fasta_prefix,fasta_data) 
            prev_chrom = test_chrom
            
        test_size = test_right - test_left
             
        usage_WT_1 = float(test_data[4]) 
        usage_WT_2 = float(test_data[5]) 
        usage_WT_3 = float(test_data[6])  
        
        usage_mut_1 = float(test_data[7])
        usage_mut_2 = float(test_data[8]) 
        usage_mut_3 = float(test_data[9])  
                
        if test_strand == "+":
            prime_5_coord = test_right - 20
            prime_3_coord = test_right + 3 
        if test_strand == "-":
            prime_5_coord = test_left - 3
            prime_3_coord = test_left + 20
                
        out_sequence_1 = get_sequence (prime_5_coord,prime_3_coord,test_strand,fasta_data)
        out_score_1 = score3(out_sequence_1,ss_matrix)
        
        out_sequence_2 = get_sequence (prime_5_coord-1,prime_3_coord-1,test_strand,fasta_data)
        out_score_2 = score3(out_sequence_2,ss_matrix)
                    
        out_sequence_3 = get_sequence (prime_5_coord+1,prime_3_coord+1,test_strand,fasta_data)
        out_score_3 = score3(out_sequence_3,ss_matrix)
        
        out_sequence_4 = get_sequence (prime_5_coord-2,prime_3_coord-2,test_strand,fasta_data)
        out_score_4 = score3(out_sequence_3,ss_matrix)            

        out_sequence_5 = get_sequence (prime_5_coord+2,prime_3_coord+2,test_strand,fasta_data)
        out_score_5 = score3(out_sequence_3,ss_matrix)
        
        if out_score_2 > out_score_1:
            out_score_1 = out_score_2
            out_sequence_1 = out_sequence_2
            
        if out_score_3 > out_score_1:
            out_score_1 = out_score_3
            out_sequence_1 = out_sequence_3
    
        if out_score_4 > out_score_1:
            out_score_1 = out_score_4
            out_sequence_1 = out_sequence_4

        if out_score_5 > out_score_1:
            out_score_1 = out_score_5
            out_sequence_1 = out_sequence_5
                          
        for sup_junction in counter_junctions :
            
            sup_SS3_mapped = 0
                        
            try:
                sup_data = sup_junction.split('\t')
            except ValueError:
                continue
            
            sup_chrom = sup_data[0]
            sup_strand = sup_data[1]
                     
            if test_strand != sup_strand or sup_chrom != test_chrom:
                continue
            
            sup_left = int(sup_data[2])
            sup_right = int(sup_data[3])
                                            
            if test_strand == "+" and sup_strand == "+" :
                if test_left - sup_left > -3 and test_left - sup_left < 3 : 
                    if not ( test_right - sup_right > -3 and test_right - sup_right < 3 ):
                                                                                                    
                            sup_usage_WT_1 =  float(sup_data[4]) 
                            sup_usage_WT_2 =  float(sup_data[5])
                            sup_usage_WT_3 =  float(sup_data[6])
                                
                            sup_usage_mut_1 =  float(sup_data[7])
                            sup_usage_mut_2 =  float(sup_data[8])
                            sup_usage_mut_3 =  float(sup_data[9])
                                
                            sup_size = sup_right - sup_left
                                  
                            sup_5_coord = sup_right - 20
                            sup_3_coord = sup_right + 3    
                                
                            sup_sequence_1 = get_sequence (sup_5_coord,sup_3_coord,test_strand,fasta_data)
                            sup_score_1 = score3(sup_sequence_1,ss_matrix)
        
                            sup_sequence_2 = get_sequence (sup_5_coord-1,sup_3_coord-1,test_strand,fasta_data)
                            sup_score_2 = score3(sup_sequence_2,ss_matrix)
                
                            sup_sequence_3 = get_sequence (sup_5_coord+1,sup_3_coord+1,test_strand,fasta_data)
                            sup_score_3 = score3(sup_sequence_3,ss_matrix)
        
                            sup_sequence_4 = get_sequence (sup_5_coord-2,sup_3_coord-2,test_strand,fasta_data)
                            sup_score_4 = score3(sup_sequence_3,ss_matrix)            
                        
                            sup_sequence_5 = get_sequence (sup_5_coord+2,sup_3_coord+2,test_strand,fasta_data)
                            sup_score_5 = score3(sup_sequence_3,ss_matrix)
        
                            if sup_score_2 > sup_score_1:
                                sup_score_1 = sup_score_2
                                sup_sequence_1 = sup_sequence_2
            
                            if sup_score_3 > sup_score_1:
                                sup_score_1 = sup_score_3
                                sup_sequence_1 = sup_sequence_3
    
                            if sup_score_4 > sup_score_1:
                                sup_score_1 = sup_score_4
                                sup_sequence_1 = sup_sequence_4

                            if sup_score_5 > sup_score_1:
                                sup_score_1 = sup_score_5
                                sup_sequence_1 = sup_sequence_5
                                
                            sup_SS3_mapped += 1  
                                                                                           
            if test_strand == "-" and sup_strand == "-" :
                if test_right - sup_right > -3 and test_right - sup_right < 3 :
                    if not ( test_left - sup_left > -3 and test_left - sup_left < 3 ):
                                                                                                                    
                            sup_usage_WT_1 = float(sup_data[4]) 
                            sup_usage_WT_2 = float(sup_data[5])
                            sup_usage_WT_3 = float(sup_data[6])
                                
                            sup_usage_mut_1 = float(sup_data[7])
                            sup_usage_mut_2 = float(sup_data[8])
                            sup_usage_mut_3 = float(sup_data[9])
                                
                            sup_size = sup_right - sup_left
                                
                            sup_5_coord = sup_left - 3
                            sup_3_coord = sup_left + 20
                                
                            sup_sequence_1 = get_sequence (sup_5_coord,sup_3_coord,test_strand,fasta_data)
                            sup_score_1 = score3(sup_sequence_1,ss_matrix)
        
                            sup_sequence_2 = get_sequence (sup_5_coord-1,sup_3_coord-1,test_strand,fasta_data)
                            sup_score_2 = score3(sup_sequence_2,ss_matrix)
                    
                            sup_sequence_3 = get_sequence (sup_5_coord+1,sup_3_coord+1,test_strand,fasta_data)
                            sup_score_3 = score3(sup_sequence_3,ss_matrix)
        
                            sup_sequence_4 = get_sequence (sup_5_coord-2,sup_3_coord-2,test_strand,fasta_data)
                            sup_score_4 = score3(sup_sequence_3,ss_matrix)            
                        
                            sup_sequence_5 = get_sequence (sup_5_coord+2,sup_3_coord+2,test_strand,fasta_data)
                            sup_score_5 = score3(sup_sequence_3,ss_matrix)
        
                            if sup_score_2 > sup_score_1:
                                sup_score_1 = sup_score_2
                                sup_sequence_1 = sup_sequence_2
            
                            if sup_score_3 > sup_score_1:
                                sup_score_1 = sup_score_3
                                sup_sequence_1 = sup_sequence_3
    
                            if sup_score_4 > sup_score_1:
                                sup_score_1 = sup_score_4
                                sup_sequence_1 = sup_sequence_4

                            if sup_score_5 > sup_score_1:
                                sup_score_1 = sup_score_5
                                sup_sequence_1 = sup_sequence_5
                                
                            sup_SS3_mapped += 1  
                                  
            if sup_SS3_mapped > 0:
                
                usage_WT = (usage_WT_1 + usage_WT_2 + usage_WT_3) / 3
                usage_mut = (usage_mut_1 + usage_mut_2 + usage_mut_3) / 3
                sup_usage_WT = (sup_usage_WT_1 + sup_usage_WT_2 + sup_usage_WT_3) / 3
                sup_usage_mut = (sup_usage_mut_1 + sup_usage_mut_2 + sup_usage_mut_3) / 3
                                                         
                out_line = test_chrom + "\t" + str(test_left)  + "\t" + str(test_right) + "\t\t" + str(usage_WT) +  "\t" + str(usage_mut) +  "\t" + str(sup_usage_WT) +  "\t" + str(sup_usage_mut) + "\t\t" + str(test_size) + "\t" + str(sup_size) + "\t\t" + str(out_score_1) + "\t" + str(sup_score_1)
                out_line = out_line.replace(".",",") 
                print >>outfile_1, out_line
                print >>sys.stderr, test_chrom + "  :  " + str(test_left) + "                     \r",  
                                                    
    print >>sys.stderr, "\n"                     
    print >>sys.stderr, "Done  "                                 
               
def main():
    usage = "%prog [options]" + "\n"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-p","--p-value",action="store",type="float",dest="cutoff",help="dept cutoff [required]")
    parser.add_option("-m","--m-file",action="store",type="string",dest="m_file",help="junctions file")
    parser.add_option("-r","--r-file",action="store",type="string",dest="r_file",help="ref exons file")
    parser.add_option("-x","--x-file",action="store",type="string",dest="x_file",help="matrix file")
    parser.add_option("-f","--fasta-file",action="store",type="string",dest="fasta_prefix",help="fasta prefix")
    parser.add_option("-o","--out-file",action="store",type="string",dest="out_prefix",help="output prefix")
    
    (options,args)=parser.parse_args()

    if not (options.cutoff and options.m_file and options.r_file and options.x_file and options.fasta_prefix and options.out_prefix):
        parser.print_help()
        sys.exit(0)
    else:
        get_3SS (options.cutoff,options.m_file,options.r_file,options.x_file,options.fasta_prefix,options.out_prefix)

if __name__ == '__main__':
     main()
    
        