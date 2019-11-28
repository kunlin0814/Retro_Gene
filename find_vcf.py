# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:52:10 2019

@author: abc73_000
"""

with open('C:/Users/abc73_000/Desktop/retro_gene_list_information.txt','r') as f:
    candidate_file = f.read()
    candidate_list =candidate_file.split('\n')[1:-1]

with open('C:/Users/abc73_000/Desktop/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS','r') as f:
    file = f.read()
    vcf_file = file.split('\n')[:-1]

cnt = 0
for i in range(len(vcf_file)):
    chr_number = vcf_file[i].split('\t')[0]
    loc = int(vcf_file[i].split('\t')[1])
    for j  in range(len(candidate_list)):
        other_chr = candidate_list[j].split('\t')[1]
        if chr_number == other_chr :
            first_pos = int(candidate_list[j].split('\t')[2])
            last_pos = int(candidate_list[j].split('\t')[3])
            if (first_pos<loc< last_pos) :
               print(chr_number)
               print(cnt)
           #if (int(candidate_list[0].split('\t')[2])<=loc<=int(candidate_list[0].split('\t')[3])):
           #    cnt+=1
       

print(cnt)




