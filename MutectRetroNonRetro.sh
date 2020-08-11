#!/bin/bash
#PBS -N retro_non_retro-${Cancer_Type}
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=30gb
#PBS -q batch


Retro_gene_source='/work/szlab/kh31516_Lab_Share_script/Retro_gene_source'
Retro_gene_script='/work/szlab/kh31516_Lab_Share_script/Retro_gene_script'
source_folder='/scratch/kh31516/Pan_cancer/TGen/Mutect'
output_folder='/scratch/kh31516/Pan_cancer/TGen/Mutect/Mutation_Rates'
Total_file_name='/scratch/kh31516/Pan_cancer/TGen/source/Tgene_pairs'
Cancer_Type='Tgen_OSA'
module load Anaconda3/2018.12
source activate py35

cd  ${source_folder}
mkdir -p ${output_folder}

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "file_name" "retro_mutation_rate" "retor_PASS" "retro_callable" "non_retro_mutation_rate" \
"non_retro_PASS" "nonr_retro_callable" "combine_mutation_rate" "combine_PASS" "combine_callable" "Cancer_Type" >> ${output_folder}/New_Mut_retro_non_retro_${Cancer_Type}_summary.txt

while read line;
do
    mkdir ${output_folder}/${line}
    cd  ${source_folder}/${line}
    #awk '$7 == "PASS" {print $0}' ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut
    printf "%s\t%s\n" "${line}" "${source_folder}/${line}/${line}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt" > ${output_folder}/${line}/sample_wigs_file.txt
    a=$(cat ${line}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt | grep "^1" | wc -l | bc) # denominator, all callable
    b=$(cat ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut | grep -v "^#" | grep -i "PASS"| wc -l |bc) #numerator, all PASS
    c=$(python $Retro_gene_script/find_vcf.py $Retro_gene_source/new_retro_gene_list_information_CanFam3.1.99gtf.txt ${source_folder}/${line}/${line}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut) # retro numerator
    java -cp $Retro_gene_script/ GetCallableCounts ${output_folder}/${line}/sample_wigs_file.txt ${output_folder}/${line}/retro_callable_summary.txt # give the retro gene callable summary
    d=$(python $Retro_gene_script/calcaulte_retro_callable.py ${output_folder}/${line}/retro_callable_summary.txt) # retro denominator
    e=$((b-c)) # non-retro numerator
    f=$((a-d)) # non-retro denominator
    g=$(echo "$((c))*1000000/$((d))" | bc -l) # retro mutation rate per million
    h=$(echo "$((e))*1000000/$((f))" | bc -l) # non-retro mutation rate per million
    i=$(echo "$((b))*1000000/$((a))" | bc -l)
    printf  "%s\t%4f\t%d\t%d\t%4f\t%d\t%d\t%4f\t%d\t%d\t%s\n" "${line}" "${g}" "${c}" "${d}" "${h}" "${e}" "${f}" \
    "${i}" "${b}" "${a}" "${Cancer_Type}" >> ${output_folder}/New_Mut_retro_non_retro_${Cancer_Type}_summary.txt
    rm -r ${output_folder}/${line}
    #rm -r ${output_folder}/${line}/retro_callable_summary.txt
done < $Total_file_name