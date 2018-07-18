#split an input GTF into several bed files
#useful for making flattened geen models


import sys
from argparse import ArgumentParser, FileType
#from operator import itemgetter

def get_bed_exons (ref_chrom,bed_file,bed_exons):
    del bed_exons[:]
    bedfile = open(bed_file,'r')
    for line in bedfile:       
        line = line.rstrip('\n')
        if line.startswith('#'):
            continue
        try:
            exon_data = line.split('\t')
        except ValueError:
            continue
            
        if exon_data[0] == ref_chrom:
            tmp_exon = [exon_data[0],exon_data[1],exon_data[2],exon_data[3],exon_data[4],exon_data[5]]
            bed_exons.append(tmp_exon)
           
    bedfile.close()       
    return bed_exons
    
def go_gene (bed_file,gene_file):

    gene = []
    exon = []
    gene_data = []
    gene_list = []
    bed_exons = [] 
    prev_chrom = ''
        
    ref_file = open(gene_file)
    print >>sys.stderr, "Getting genes ..."
    
    for line in ref_file:
        line = line.rstrip('\n')
        if line.startswith('#'):
            continue
        try:
            gene_data = line.split('\t')
        except ValueError:
            continue 
        gene = [gene_data[0],gene_data[1],gene_data[2],gene_data[3],gene_data[4],gene_data[5]]
        gene_list.append(gene)
    
    bedfile = open(bed_file,'r')
    out_file_tmp = bed_file.split(".")
    out_file = out_file_tmp[0] + '.spl_genes.bed'
    outfile = open(out_file,'w')
    
    for gene in gene_list:
        
        num_exons = 0
        ref_chrom = gene[0].rstrip("\t")

        if ref_chrom != prev_chrom:
            get_bed_exons(ref_chrom,bed_file,bed_exons)
            prev_chrom = ref_chrom
            
        for exon in bed_exons:
            if exon[3] == gene[3] and exon[5] == gene[5]:
                if int(exon[1]) >= int(gene[1]) - 2 and int(exon[2]) <= int(gene[2]) + 2:
                    num_exons = num_exons + 1
                    
        if num_exons > 1:
            print >>outfile, '\t'.join((gene[0],gene[1],gene[2],gene[3],gene[4],gene[5]))
            print >>sys.stderr, "Spliced gene : " + gene[3] + "          \r",
            
    outfile.close()
    return()

def reduce_exons(tmp2_exons,min_exon_size,sign_str):
    
    exon_2 = []
    exon_6 = []
    exon_7 = []
    tmp6_exons = []
    tmp7_exons = []
            
    new_num_exons = len(tmp2_exons)
    old_num_exons = 0
    
    while new_num_exons > old_num_exons:
                
        for exon_2 in tmp2_exons:
            tmp6_exons.append(exon_2)
        
        old_num_exons = len(tmp2_exons)        
        get_index = 0

        for exon_2 in tmp2_exons:
            
            out_stuff = sign_str + exon_2[3] + "                 \r"
            print >>sys.stderr, out_stuff,
            
            num_equals = 0
            ref_chrom = exon_2[0]
            ref_left = exon_2[1]
            ref_right = exon_2[2]
            ref_name = exon_2[3]
            ref_trans = exon_2[4]
            ref_strand = exon_2[5]
            ref_strand = ref_strand[0]
            
            ex2_left = int(ref_left)
            ex2_right = int(ref_right)
            
            for exon_6 in tmp6_exons:
                test_chrom = exon_6[0]
                test_left = exon_6[1]
                test_right = exon_6[2]
                test_name = exon_6[3]
                test_trans = exon_6[4]
                test_strand = exon_6[5]
                test_strand = test_strand[0]
        
                ex6_left = int(test_left)
                ex6_right = int(test_right)
                    
                if test_strand == ref_strand and ex6_left == ex2_left and ex6_right == ex2_right:
                    if num_equals > 0:
                        del tmp2_exons[get_index]
                        ex2_left = 0
                        ex2_right = 0
                        old_num_exons = 0
                    else:
                        num_equals = num_equals + 1
                    
                if test_strand == ref_strand and ((ex6_left == ex2_left and ex6_right > ex2_right) or (ex6_left < ex2_left and ex6_right == ex2_right)):
                    del tmp2_exons[get_index]
                    ex2_left = 0
                    ex2_right = 0
                    old_num_exons = 0
                    continue
                                   
                if test_strand == ref_strand and test_name == ref_name and ex6_left < ex2_left and ex6_right > ex2_right:
                    del tmp2_exons[get_index]
                    ex2_left = 0
                    ex2_right = 0
                    old_num_exons = 0
                    continue
                         
                if sign_str.find("red") > -1  and test_name == ref_name and ex6_left < ex2_left and ex6_right < ex2_right and ex6_right > (ex2_left - min_exon_size):
                    new_left_str = str(ex6_left)
                    new_right_str = str(ex2_right)
                    
                    new_exon =[ref_chrom,new_left_str,new_right_str,ref_name,ref_trans,ref_strand]
                                                            
                    del tmp2_exons[get_index]
                    ex2_left = 0
                    ex2_right = 0
                    tmp7_exons.append(new_exon)
                    old_num_exons = 0
                    continue
                    
                if sign_str.find("red") > -1 and test_name == ref_name and ex6_right > ex2_right and ex6_left > ex2_left and ex6_left < (ex2_right + min_exon_size):
                    new_right_str = str(ex6_right)
                    new_left_str = str(ex2_left)
                    
                    new_exon =[ref_chrom,new_left_str,new_right_str,ref_name,ref_trans,ref_strand]
                                        
                    del tmp2_exons[get_index]
                    ex2_left = 0
                    ex2_right = 0
                    tmp7_exons.append(new_exon)
                    old_num_exons = 0
                    continue
                
                if sign_str.find("map") > -1  and test_strand == ref_strand and ex6_left < ex2_left and ex6_right < ex2_right and ex6_right > (ex2_left - min_exon_size):
                    new_left_str = str(ex6_left)
                    new_right_str = str(ex2_right)
                    
                    new_exon =[ref_chrom,new_left_str,new_right_str,ref_name,ref_trans,ref_strand]
                                                            
                    del tmp2_exons[get_index]
                    ex2_left = 0
                    ex2_right = 0
                    tmp7_exons.append(new_exon)
                    old_num_exons = 0
                    continue
                    
                if sign_str.find("map") > -1 and test_strand == ref_strand and ex6_right > ex2_right and ex6_left > ex2_left and ex6_left < (ex2_right + min_exon_size):
                    new_right_str = str(ex6_right)
                    new_left_str = str(ex2_left)
                    
                    new_exon =[ref_chrom,new_left_str,new_right_str,ref_name,ref_trans,ref_strand]
                                        
                    del tmp2_exons[get_index]
                    ex2_left = 0
                    ex2_right = 0
                    tmp7_exons.append(new_exon)
                    old_num_exons = 0
                    continue
        
            get_index = get_index + 1
                            
        for exon_7 in tmp7_exons:
            tmp2_exons.append(exon_7)
                    
        new_num_exons = len(tmp2_exons)
        del tmp7_exons[:]
        del tmp6_exons[:]

    tmp2_exons = sorted (tmp2_exons,key=lambda item: int(item[1])) 
    return(tmp2_exons)

def gtf_splitter (do_dist,in_file):
    
    exon_1 = []
    exon_2 = []
    exon_3 = []
    exon_4 = []
    exon_5 = []
    exon = []    
    gene_data = []
    tmp1_exons = []
    tmp2_exons = []
    tmp3_exons = []
    tmp4_exons = []
    tmp5_exons = []
    buff_exons = []
    data_string = ''
    test_name = ''
    test_dept = '255'
    ref_dept = '255'    
    prev_chrom = ''
    prev1_chrom = ''
    prev2_chrom = ''   
    
    num_lines = 0
    min_exon_size = int(do_dist)
    
    infile = open(in_file)
    out_file = in_file + '.exons.bed'
    outfile = open(out_file,'w')
    out2_file = in_file + '.unique.bed'
    out2file = open(out2_file,'w')
    
    for line in infile:
       
        line = line.rstrip('\n')
        
        if line.startswith('#'):
            continue
        
        try:
            exon_1 = line.split('\t')
        except ValueError:
            continue
        
        test_chrom = exon_1[0]
        someval_2 = exon_1[2]
        test_left = str(int(exon_1[3])-1)
        test_right = exon_1[4]
        test_strand = exon_1[6]
        
        gene_data = exon_1[8].split(';')
        
        test_name = ''
        trans_name = ''
        for data_string in gene_data:
            if data_string.find("gene_id") > -1:
                tmp_1 = data_string
                tmp_1 = tmp_1.replace("gene_id","")
                tmp_1 = tmp_1.replace(" ","")
                tmp_1 = tmp_1.lstrip('"')
                test_name = tmp_1.rstrip('"')
            if data_string.find("transcript_id") > -1:       
                tmp_11 = data_string
                tmp_11 = tmp_11.replace("transcript_id","")
                tmp_11 = tmp_11.replace(" ","")
                tmp_11 = tmp_11.lstrip('"')
                trans_name = tmp_11.rstrip('"')
           
        exon_size = int(test_right) - int(test_left)
        test_strand = test_strand[0]
        
        if someval_2.find("exon") > -1 and exon_size > min_exon_size:
            exon_2 = [test_chrom,test_left,test_right,test_name,trans_name,test_strand]  
                       
            tmp1_exons.append(exon_2)
            tmp2_exons.append(exon_2)
                            
            exon_3 = [test_chrom,test_left,test_right,test_name,trans_name,test_strand]
            
            if test_chrom != prev1_chrom:
                del tmp3_exons[:]
                prev1_chrom = test_chrom
              
            try:
                tmp_exon = tmp3_exons.index(exon_3)
            except ValueError:
                tmp3_exons.append(exon_3)
                num_lines = num_lines + 1
                out_line = test_chrom +  "\t" + test_left +  "\t" +  test_right +  "\t" + trans_name +  "\t" + test_dept +  "\t" + test_strand 
                print >>outfile, out_line
                out_stuff = "exons: " + str(num_lines)  +  "     \r"
                print >>sys.stderr, out_stuff,
                
            exon_4 = [test_chrom,test_left,test_right,test_strand]
            exon_5 = [test_chrom,test_left,test_right,test_name,test_dept,test_strand]
            
            if test_chrom != prev2_chrom:
                del tmp4_exons[:]
                prev2_chrom = test_chrom
            
            try:
                tmp_exon = tmp4_exons.index(exon_4)
            except ValueError:
                tmp4_exons.append(exon_4)
                tmp5_exons.append(exon_5)                
                out_line = test_chrom +  "\t" + test_left +  "\t" +  test_right +  "\t" + test_name +  "\t" + test_dept +  "\t" + test_strand 
                print >>out2file, out_line
                     
    out_stuff = "exons: " + str(num_lines)  +  "     \n"
    print >>sys.stderr, out_stuff  
                 
    infile.close()
    outfile.close()
    out2file.close()
    
    tmp1_exons = sorted (tmp1_exons,key=lambda item: int(item[1]))
    tmp1_exons = sorted (tmp1_exons,key=lambda item: item[0])    
    tmp2_exons sorted (tmp2_exons,key=lambda item: int(item[1]))
    tmp2_exons = sorted (tmp2_exons,key=lambda item: item[0])
    tmp5_exons = sorted (tmp5_exons,key=lambda item: int(item[1]))
    tmp5_exons = sorted (tmp5_exons,key=lambda item: item[0])
        
    del tmp3_exons[:]
    del tmp4_exons[:]
    del buff_exons[:]    
    num_lines = 0
    prev_chrom = '' 
        
    out3_file = in_file + '.trans.bed'
    out3file = open(out3_file,'w')
    
    for exon_1 in tmp1_exons:
                
        ref_chrom = exon_1[0]
        ref_left = exon_1[1]
        ref_right = exon_1[2]
        ref_name = exon_1[3]
        ref_trans = exon_1[4]
        ref_strand = exon_1[5]
        ref_strand = ref_strand[0]
            
        ex1_left = int(ref_left)
        ex1_right = int(ref_left)
        
        new_left = ref_left
        tmp_left = int(ref_left)
        new_right = ref_right
        tmp_right = int(ref_right)
            
        if ref_chrom != prev_chrom:
            prev_chrom = ref_chrom
            del buff_exons[:]
            for exon_2 in tmp2_exons:
                test_chrom = exon_2[0]
                if test_chrom == prev_chrom:
                    buff_exons.append(exon_2)
                
        for exon in buff_exons:
            test_chrom = exon[0]
            test_left = exon[1]
            test_right = exon[2]
            test_name = exon[3]
            test_trans = exon[4]
            test_strand = exon[5]
            test_strand = test_strand[0]            
                        
            if test_chrom == ref_chrom and test_trans == ref_trans and ref_name == test_name and test_strand == ref_strand:
                ex2_left = int(test_left)
                ex2_right = int(test_right)
                
                if ex2_left  < tmp_left:
                    tmp_left = ex2_left 
       
                if ex2_right > tmp_right:
                    tmp_right = ex2_right
        
        new_left = str(tmp_left)
        new_right = str(tmp_right)
               
        exon_3 = [ref_chrom,new_left,new_right,ref_name,ref_trans,ref_strand]  
        
        try:
            tmp_exon = tmp3_exons.index(exon_3)
        except ValueError:
            tmp3_exons.append(exon_3)
            num_lines = num_lines + 1
            out_line = ref_chrom +  "\t" + new_left +  "\t" +  new_right +  "\t" + ref_trans +  "\t" + test_dept +  "\t" + ref_strand 
            print >>out3file, out_line
            out_stuff = "trans: " + str(num_lines)  +  "     \r"
            print >>sys.stderr, out_stuff,
    
    out3file.close()        
    out_stuff = "trans: " + str(num_lines)  +  "     \n"
    print >>sys.stderr, out_stuff
    
    num_lines = 0
    
    out_file = in_file + '.genes.bed'
    out2file = open(out_file,'w')
    
    del tmp1_exons[:]
    del tmp2_exons[:]    
    del buff_exons[:]  
    prev_chrom = '' 
    
    for exon_3 in tmp3_exons:
        tmp1_exons.append(exon_3)
        
    for exon_1 in tmp1_exons:
        ref_chrom = exon_1[0]
        if ref_chrom != prev_chrom:
            del tmp2_exons[:]
            prev_chrom = ref_chrom
            for exon_3 in tmp3_exons:
                test_chrom = exon_3[0]
                if test_chrom == ref_chrom:
                    tmp2_exons.append(exon_3)
            reduce_exons(tmp2_exons,min_exon_size,"reducing ... ")
    
            for exon_2 in tmp2_exons:
                buff_exons.append(exon_2)
            
    del tmp3_exons[:]

    for exon in buff_exons:
        tmp3_exons.append(exon)
    
    tmp3_exons = sorted (tmp3_exons,key=lambda item: int(item[1]))
    tmp3_exons = sorted (tmp3_exons,key=lambda item: item[0])
    
    del tmp1_exons[:]
    del tmp2_exons[:]    
    del buff_exons[:]
    prev_chrom = ''
    
    print >>sys.stderr, "\n"
    
    for exon_3 in tmp3_exons:
        
        tmp_left = 100000000000
        tmp_right = 0
        out_line = '' 
        ref_chrom = exon_3[0]
        ref_left = exon_3[1]
        ref_right = exon_3[2]
        ref_name = exon_3[3]
        ref_trans = exon_3[4]
        ref_strand = exon_3[5]
        ref_strand = ref_strand[0]
            
        new_left = ref_left
        tmp_left = int(ref_left)
        new_right = ref_right
        tmp_right = int(ref_right)
        
        if ref_chrom != prev_chrom:
            prev_chrom = ref_chrom
            del buff_exons[:]
            for exon_4 in tmp4_exons:
                test_chrom = exon_4[0]
                if test_chrom == prev_chrom:
                    buff_exons.append(exon_4)   
                                                       
        for exon in buff_exons:
            test_chrom = exon[0]
            test_left = exon[1]
            test_right = exon[2]
            test_name = exon[3]
            test_trans = exon[4]
            test_strand = exon[5]
            test_strand = test_strand[0]
                        
            if test_chrom == ref_chrom and test_strand == ref_strand and ref_name == test_name:
                ex2_left = int(test_left)
                ex2_right = int(test_right)
                  
                if ex2_left < tmp_left and ex2_right > tmp_right and ref_name == test_name:
                    #ref_name = test_name
                    tmp_left = ex2_left
                    tmp_right = ex2_right
                    new_left = str(tmp_left)
                    new_right = str(tmp_right)

                if ex2_left < tmp_left and ex2_right >= tmp_left and ex2_right <= tmp_right:
                    tmp_left = ex2_left
                    new_left = str(tmp_left)
           
                if ex2_right > tmp_right and ex2_left <= tmp_right and ex2_left >= tmp_left:
                    tmp_right = ex2_right
                    new_right = str(tmp_right)
                    
        exon_2 = [ref_chrom,new_left,new_right]
        exon_1 = [ref_chrom,new_left,new_right,ref_name,ref_dept,ref_strand]
                                                    
        try:
            tmp_exon = tmp2_exons.index(exon_2)
        except ValueError:
            tmp2_exons.append(exon_2)
            tmp1_exons.append(exon_1)
            out_line = ref_chrom +  "\t" + new_left +  "\t" +  new_right +  "\t" + ref_name +  "\t" + ref_dept +  "\t" + ref_strand 
            print >>out2file, out_line
            num_lines = num_lines + 1
            out_stuff = "genes: " + str(num_lines)  +  "     \r"
            print >>sys.stderr, out_stuff,
        
    out_stuff = "genes: " + str(num_lines)  +  "     \n"
    print >>sys.stderr, out_stuff
                        
    out2file.close()
        
    out_file = in_file + '.map_exons.bed'
    out2file = open(out_file,'w')
    out_file = in_file + '.spl_exons.bed'
    out3file = open(out_file,'w')
    
    prev_chrom = ''
    
    for exon_1 in tmp1_exons:
        
        del tmp2_exons[:]     
        
        ref_chrom = exon_1[0]
        ref_left = exon_1[1]
        ref_right = exon_1[2]
        ref_name = exon_1[3]
        ref_trans = exon_1[4]
        ref_strand = exon_1[5]
        ref_strand = ref_strand[0]
        
        ex1_left = int(ref_left)-2
        ex1_right = int(ref_right)+2
        
        if ref_chrom != prev_chrom:
            del buff_exons[:]
            prev_chrom = ref_chrom
            for exon_5 in tmp5_exons:
                test_chrom = exon_5[0]
                if test_chrom == ref_chrom:
                    buff_exons.append(exon_5)
                            
        for exon in buff_exons:
            test_chrom = exon[0]
            test_left = exon[1]
            test_right = exon[2]
            test_name = exon[3]
            test_trans = exon[4]
            test_strand = exon[5]
            test_strand = test_strand[0]
            
            buff_left = int(test_left)
            buff_right = int(test_right)
             
            if test_name == ref_name and test_strand == ref_strand and buff_left > ex1_left and buff_right < ex1_right:
                tmp2_exons.append(exon)
 
        reduce_exons(tmp2_exons,min_exon_size,"mapping ... ")
        
        tmp2_exons = sorted (tmp2_exons,key=lambda item: int(item[1]))
        tmp2_exons = sorted (tmp2_exons,key=lambda item: item[0])
        
        for exon_2 in tmp2_exons:        
            out_line = exon_2[0] +  "\t" + exon_2[1] +  "\t" + exon_2[2] +  "\t" + exon_2[3] +  "\t" + exon_2[4] +  "\t" + exon_2[5]
            print >>out2file, out_line
            
        if len(tmp2_exons) > 1:
            for exon_2 in tmp2_exons:      
                out_line = exon_2[0] +  "\t" + exon_2[1] +  "\t" + exon_2[2] +  "\t" + exon_2[3] +  "\t" + exon_2[4] +  "\t" + exon_2[5]
                print >>out3file, out_line
           
    out2file.close()
    out3file.close()
    
    # added for getting spliced genes
    bed_file = in_file + '.spl_exons.bed'
    gene_file = in_file + '.genes.bed'
    go_gene (bed_file,gene_file)
    
    print >>sys.stderr, "\n"                
    print >>sys.stderr, "All done. \n"
    
# argument parser
if __name__ == '__main__':
    parser = ArgumentParser(
        description='split a GTF in idividual BED files')  
    parser.add_argument('do_dist', nargs='?', type=str, help='tolerance')   
    parser.add_argument('in_file',nargs='?',type=str, help='input gtf')     
    args = parser.parse_args()
    if not args.do_dist and args.in_file:
            parser.print_help()
            exit(1)
    gtf_splitter (args.do_dist,args.in_file)
    
        