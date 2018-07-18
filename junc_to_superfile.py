
#this script will collect junctions from two inputs and output a single flattened file 
#this script is used to combine junctions from multiple datasets into a single superfile
#for multiple samples, run this script stepwise, using the output from prior runs as one input
#input junction smust be ordered by chromosome and left coordinate

#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
     print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " this script needs python2.7!\n"
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

# get junctions from tophat2 junctions.bed file and put values in same order as junc1s. 

        
def open_junctions (junc1_chrom,juncs_file,test3_junctions):
     
    del test3_junctions[:]
        
    infile = open(juncs_file,"r")
    for junc2_line in infile:
        
        junc2_line = junc2_line.rstrip('\n')
        if junc2_line.startswith("chr"):

            try:
                junc2_data = junc2_line.split('\t') 
            except ValueError:
                #print >>sys.stderr,"skipped line"
                continue
        
            junc2_chrom = junc2_data[0]
            if junc2_chrom.upper() == junc1_chrom.upper():
                test3_junctions.append(junc2_line)
            if junc1_chrom.upper() == "ALL":
                test3_junctions.append(junc2_line)
                    
    infile.close ()
    return test3_junctions
    
    
def get_fasta (juncs_file,test_file,out_file):
    
    list_data = []
    junction = []

    fasta_data = []
    test1_junctions = []
    test2_junctions = []
    out_series = []
    prev1_chrom = 'abc' 
    prev2_chrom = 'abc' 
    all_junctions = []
    test_junctions=[]
    out_junctions =[]
    leftover_junctions = []
    used_junctions = []
    
    outfile = open(out_file,'w')
    
    open_junctions("all",juncs_file,test1_junctions)    
        
    for junc34 in test1_junctions:
        all_junctions.append(junc34)
            
    for junc1 in all_junctions:
                
        try:
            junc1_chrom, junc1_left, junc1_right, junc1_name, junc1_depth, junc1_strand, someval1_1, someval1_2, colorstuff1, someval1_3, overhangs1, someval1_4  = junc1.split('\t')
        except ValueError:
            continue
            
        out_junction = junc1
        got_hit = 0
            
        #print >>sys.stderr, junc1_chrom + "      " + junc1_left + "       " + junc1_right + "   \r"
                            
        if len(junc1_chrom) > 5:
            continue 
        
        try:
             overhang1_left, overhang1_right = overhangs1.split(',')
        except ValueError:
            continue
                                            
        junction_1 = [junc1_chrom,junc1_left,junc1_right,junc1_name,junc1_depth,junc1_strand,overhang1_left,overhang1_right]  
        
        junc1_chrom = junction_1[0] 
        junc1_leftcoord = int(junction_1[1]) + int(junction_1[6])
        junc1_rightcoord = int(junction_1[2]) - int(junction_1[7])
        junc1_name = junction_1[3]
        junc1_depth_val = int(junction_1[4])
        junc1_strand = junction_1[5]
        junc1_strand = junc1_strand[0]
        
        if ((junc1_chrom != prev1_chrom) and (prev1_chrom != 'abc')):
            del test2_junctions[-1]
            for test2_junction in test2_junctions:
                leftover_junctions.append(test2_junction)
                            
        if junc1_chrom != prev1_chrom:
            del test2_junctions[:]
            open_junctions(junc1_chrom,test_file,test2_junctions)
            prev1_chrom = junc1_chrom
            
        j_counter = 0
        
        out_overhang_left = overhang1_left
        out_overhang_right = overhang1_right
        out_junc_left = junc1_left
        out_junc_right = junc1_right
        out_someval_1 = someval1_1
        out_someval_2 = someval1_2
        out_someval_3 = someval1_3
        out_someval_4 = someval1_4
        
        while len(test2_junctions) > 0 and j_counter < len(test2_junctions)-1:
                        
            junction_tmp = test2_junctions[j_counter]
            
            try:
                junc2_chrom, junc2_left, junc2_right, junc2_name, junc2_depth, junc2_strand, someval2_1, someval2_2, colorstuff2, someval2_3, overhangs2, someval2_4  = junction_tmp.split('\t')
            except ValueError:
                continue
        
            try:
                 overhang2_left, overhang2_right = overhangs2.split(',')
            except ValueError:
                continue
            
            junction_2 = [junc2_chrom,junc2_left,junc2_right,junc2_name,junc2_depth,junc2_strand,overhang2_left,overhang2_right]  
            
            junc2_chrom = junction_2[0] 
            junc2_leftcoord = int(junction_2[1]) + int(junction_2[6])
            junc2_rightcoord = int(junction_2[2]) - int(junction_2[7])
            junc2_name = junction_2[3]
            junc2_depth_val = int(junction_2[4])
            junc2_strand = junction_2[5]
            junc2_strand = junc2_strand[0]
            
            #overhang2_total = overhang2_left + overhang2_right
                        
            if junc2_rightcoord < (junc1_leftcoord - 500000):
                try:
                    used_index = used_junctions.index(test2_junctions[0])
                except ValueError:
                    try:
                        used_index = leftover_junctions.index(test2_junctions[0])
                    except ValueError:
                        leftover_junctions.append(test2_junctions[0])
                del test2_junctions[0]
                j_counter = 0
                continue
            
            if junc2_leftcoord > (junc1_rightcoord + 500000):
                j_counter = len(test2_junctions) + 10
                continue
                            
            if junc2_chrom.upper() == junc1_chrom.upper() and junc2_strand[0] == junc1_strand[0]:
                if junc1_leftcoord - junc2_leftcoord > -3 and junc1_leftcoord - junc2_leftcoord < 3:
                    if junc1_rightcoord - junc2_rightcoord > -3 and junc1_rightcoord - junc2_rightcoord < 3:
                        total_depth_val = junc1_depth_val + junc2_depth_val
                        
                        if int(overhang2_left) > int(overhang1_left):
                            out_overhang_left = overhang2_left
                            out_junc_left = junc2_left

                        if int(overhang2_right) > int(overhang1_right):
                            out_overhang_right = overhang2_right
                            out_junc_right = junc2_right

                        out_overhangs = out_overhang_left + "," + out_overhang_right
                        out_junction = junc1_chrom + "\t" + out_junc_left + "\t" + out_junc_right + "\t" + junc1_name + "\t" + str(total_depth_val) + "\t" + junc1_strand[0] + "\t" + out_someval_1 + "\t" + out_someval_2 + "\t" + colorstuff1 + "\t" + out_someval_3 + "\t" + out_overhangs + "\t" + out_someval_4
                        
                        print >>outfile, out_junction
                        used_junctions.append(test2_junctions[j_counter])                       
                        del test2_junctions[j_counter]
                        j_counter = len(test2_junctions)
                        got_hit = 1
                                                
            j_counter += 1
            
        if got_hit < 1:
            print >>outfile, out_junction
                
        print >>sys.stderr, junc1_chrom + "      " + str(junc1_left) + "       " + str(len(test2_junctions) + len(used_junctions) + len(leftover_junctions)) + "       " + str(len(leftover_junctions)) + "                 \r",
    
    # for last chromosome
    del test2_junctions[-1]
    for test2_junction in test2_junctions:
        leftover_junctions.append(test2_junction)
        
    # process leftovers     
    print >>sys.stderr, junc1_chrom + "      " + str(junc1_left) + "       " + str(len(test2_junctions) + len(used_junctions) + len(leftover_junctions)) + "       " + str(len(leftover_junctions)) + "                 \r",
    print >>sys.stderr, "\n"      
    
    prev_chrom = "abc"
               
    for junction_leftover in leftover_junctions:
        
        got_hit = 0
        
        try:
            junc1_chrom, junc1_left, junc1_right, junc1_name, junc1_depth, junc1_strand, someval1_1, someval1_2, colorstuff1, someval1_3, overhangs1, someval1_4  = junction_leftover.split('\t')
        except ValueError:
            continue
                    
        junction_1 = [junc1_chrom,junc1_left,junc1_right,junc1_name,junc1_depth,junc1_strand,overhang1_left,overhang1_right]  
        
        junc1_chrom = junction_1[0] 
        junc1_leftcoord = int(junction_1[1]) + int(junction_1[6])
        junc1_rightcoord = int(junction_1[2]) - int(junction_1[7])
        junc1_name = junction_1[3]
        junc1_depth_val = int(junction_1[4])
        junc1_strand = junction_1[5]
        junc1_strand = junc1_strand[0]
        
        if junc1_chrom != prev_chrom:
            del test1_junctions[:]
            for all2_junction in all_junctions:
                try:
                    junc2_chrom, junc2_left, junc2_right, junc2_name, junc2_depth, junc2_strand, someval2_1, someval2_2, colorstuff2, someval2_3, overhangs2, someval2_4  = all2_junction.split('\t')
                except ValueError:
                    continue
                if junc2_chrom.upper() == junc1_chrom.upper():
                    test1_junctions.append(all2_junction)
            prev_chrom = junc1_chrom
                        
        for fin_junction in test1_junctions:
            
            if got_hit > 0:
                continue
                
            try:
                junc2_chrom, junc2_left, junc2_right, junc2_name, junc2_depth, junc2_strand, someval2_1, someval2_2, colorstuff2, someval2_3, overhangs2, someval2_4  = fin_junction.split('\t')
            except ValueError:
                continue
            
            junction_2 = [junc2_chrom,junc2_left,junc2_right,junc2_name,junc2_depth,junc2_strand,overhang2_left,overhang2_right] 
                            
            junc2_chrom = junction_2[0] 
            junc2_leftcoord = int(junction_2[1]) + int(junction_2[6])
            junc2_rightcoord = int(junction_2[2]) - int(junction_2[7])
            junc2_name = junction_2[3]
            junc2_depth_val = int(junction_2[4])
            junc2_strand = junction_2[5]
            junc2_strand = junc2_strand[0]
            
            if junc2_chrom.upper() == junc1_chrom.upper() and junc2_strand[0] == junc1_strand[0]:
                if junc1_leftcoord - junc2_leftcoord > -3 and junc1_leftcoord - junc2_leftcoord < 3:
                    if junc1_rightcoord - junc2_rightcoord > -3 and junc1_rightcoord - junc2_rightcoord < 3:
                        got_hit = 10
                        
        if got_hit < 1:
            print >>sys.stderr, junc1_chrom + "      " + str(junc1_left) + "                     \r",
            print >>outfile, junction_leftover 
                                    
    print >>sys.stderr, "\n",           
    print >>sys.stderr, "Done  "

def main():
    usage = "%prog [options]" + "\n"
    parser = OptionParser(usage,version="%prog " + __version__)   
    parser.add_option("-j","--juncs-file",action="store",type="string",dest="juncs_file",help="input junctions")
    parser.add_option("-t","--test-file",action="store",type="string",dest="test_file",help="test junctions")
    parser.add_option("-o","--out-file",action="store",type="string",dest="out_file",help="output file")   

    (options,args)=parser.parse_args()

    if not (options.juncs_file and options.test_file and options.out_file):
        parser.print_help()
        sys.exit(0)
    else:
        get_fasta (options.juncs_file,options.test_file,options.out_file)

if __name__ == '__main__':
     main()
    
        