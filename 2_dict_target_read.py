#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 13:38:09 2020

@author: kun-linho
"""
import sys
from collections import Counter
from datetime import datetime
# each chromosome
# min = 23840680 
# max = 123780325

## use 50 million bp as a cut off
now = datetime. now()
sam_file='/Volumes/Research_Data/Pan_cancer/CMT-SRR7780976-test2.sam'
file_name='2_dict_CMT-SRR7780976'
read_length= 101

## criteria is Reads start site <=Exon end site and exon start site <= reads end site
def withinRegion(reads_position, exom_start, exom_end):
     reads_end = reads_position + read_length-1
     reads_start = reads_position
     if (reads_start <= exom_end and exom_start <= reads_end):
         return 1
     elif (reads_start > exom_end):
         return 0
     elif (reads_end < exom_start):
         return -1


def binarySearch (arr, left, right, reads_position, read_length): 
    
    # Check base case 
    if right >= left: 
        mid = int((left+right)/2)
        exon_start_value = int (arr[mid][0])
        exon_end_value = int (arr[mid][1])
        # If element is present at the middle itself 
        if (withinRegion(reads_position, exon_start_value, exon_end_value ) == 1):
            return mid 
          
        # If element is smaller than mid, then it can only 
        # be present in left subarray 
        elif (withinRegion(reads_position, exon_start_value, exon_end_value ) == -1):
            return binarySearch(arr, left, mid-1, reads_position, read_length) 
  
        # Else the element can only be present in right subarray 
        elif (withinRegion(reads_position, exon_start_value, exon_end_value ) == 0) : 
            return binarySearch(arr, mid+1, right, reads_position, read_length) 
  
    else: 
        # Element is not present in the array 
        return -2


small_interval_dict ={}
large_interval_dict ={}
with open('/Volumes/Research_Data/Pan_cancer/Mapping_source/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list', 'r') as f:
    file = f.read()


CDS = file.split('\n')[:-1]
    
for i in range(len(CDS)):
    chrom = CDS[i].split(':')[0]
    start = int(CDS[i].split(':')[1].split('-')[0])
    end = int(CDS[i].split(':')[1].split('-')[1])
    if (start >= 50000000 and end >= 50000000):
        if chrom not in large_interval_dict.keys():
            large_interval_dict[chrom]=[(start, end)]
        else:
            large_interval_dict[chrom].append((start, end))
    else: 
         if chrom not in small_interval_dict.keys():
            small_interval_dict[chrom]=[(start, end)]
         else:
            small_interval_dict[chrom].append((start, end))

for i in small_interval_dict.keys():
    small_interval_dict[i].sort()
   
for i in large_interval_dict.keys():
     large_interval_dict[i].sort()

unique = 0 #3
#duplicate = 0 #3
#Onemapped = 0 #5,9
#incorrect = 0 #1
#unmapped = 0 #13

transcript_list =[]    
total = 0 
pass_line =0
with open(sam_file,'r') as f1:   
    for line in f1:
        file_lst = line.split('\t')
        if '@' in file_lst[0]:
            pass
        else :
            total += 1
            reads_name = file_lst[0]
            reads_chr = file_lst[2]
            reads_position = int(file_lst[3])
            status = int(file_lst[1])%16
            if status == 3:
                for ele in file_lst:
                    if 'XT:' in ele:
                        status2 = ele.split(':')[2]
                        if status2 == 'U' or status2 == 'M':
                            unique += 1
                            if reads_position < 50000000:
                                exome_loc = small_interval_dict[reads_chr]
                                if (binarySearch(exome_loc,0, len(exome_loc)-1,reads_position,read_length)!=-2):
                                    transcript_list.append(reads_name)
                                    pass_line+=1
                                    
                            else:
                                 exome_loc = large_interval_dict[reads_chr]
                                 if (binarySearch(exome_loc,0, len(exome_loc)-1,reads_position,read_length)!=-2):
                                    transcript_list.append(reads_name)
                                    pass_line+=1
                                

dup = [key for (key, value) in Counter(transcript_list).items() if value > 1 and key]
pairs = total / 2

summary = open('/Users/kun-linho/Desktop/'+file_name+'_CDS_mapping.txt','w')
summary.write('Total_reads'+'\t'+'Total_uniq_reads\t'+'uniq_mapped_rate\t'+'Total_read_pairs\t'+'uniq_CDS_region\t'+'uniq_CDS_region_paris_rates\t'+'CDS_uniq_pairs\n')
summary.write(str(total)+'\t'+str(unique)+'\t'+str(unique/total)+'\t'+str(pairs)+'\t'+str(pass_line)+'\t'+str(pass_line/unique)+'\t'+str(len(dup))+'\n')
summary.close()
end = datetime. now()
print(end - now)
"""
for i in dup:
    duplist.write(i+'\n')
    
duplist.close()


allList = open('/scratch/kh31516/Original_Melanoma/'+file_name+'_all_CDS_Mapping_list.txt','w')
for i in transcript_list:
    allList.write(i+'\n')

allList.close()
"""  