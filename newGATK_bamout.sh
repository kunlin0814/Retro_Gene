#PBS -S /bin/bash
#PBS -N LAB3_dam_SAMEA4672968_bamout
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=50gb
#PBS -q batch
#PBS -m ae
#PBS -M kh31516@uga.edu


result='/scratch/kh31516/Pan_cancer/others/PRJEB26319/results/LAB3_dam_SAMEA4672968' #
reference='/scratch/kh31516/Melanoma/Melanoma_source' 
dataN='/scratch/kh31516/Pan_cancer/others/PRJEB26319/data/SAMEA4672967'
dataT='/scratch/kh31516/Pan_cancer/others/PRJEB26319/data/SAMEA4672968'
sraRunTable='/scratch/kh31516/Pan_cancer/others/PRJEB26319/source/PRJEB26319'
annovar_index='/scratch/kh31516/Melanoma_source/annovar'


module load SAMtools/1.9-foss-2016b
module load SRA-Toolkit/2.9.1-centos_linux64
module load BWA/0.7.17-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load GATK/3.8-1-Java-1.8.0_144
#module load MuTect/1.1.7-Java-1.7.0_80
module load Perl/5.26.1-GCCcore-6.4.0

: "
####### Download #######
mkdir $result
mkdir $dataN
mkdir $dataT
cd $dataN
#fastq-dump --split-files --gzip SAMEA4672967
/scratch/kh31516/Pan_cancer/others/PRJEB26319/scripts/download_sra_sample.sh SAMEA4672967 $sraRunTable
cd $dataT
#fastq-dump --split-files --gzip SAMEA4672968
/scratch/kh31516/Pan_cancer/others/PRJEB26319/scripts/download_sra_sample.sh SAMEA4672968 $sraRunTable


####### BWA Mapping #######
cd $result
# Tumor
time bwa aln $reference/canFam3.fa $dataT/SAMEA4672968_1.fastq.gz > $result/SAMEA4672968_1.sai
time bwa aln $reference/canFam3.fa $dataT/SAMEA4672968_2.fastq.gz > $result/SAMEA4672968_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SAMEA4672968_1.sai \
$result/SAMEA4672968_2.sai \
$dataT/SAMEA4672968_1.fastq.gz \
$dataT/SAMEA4672968_2.fastq.gz \
> $result/SAMEA4672968.sam

# Normal
time bwa aln $reference/canFam3.fa $dataN/SAMEA4672967_1.fastq.gz > $result/SAMEA4672967_1.sai
time bwa aln $reference/canFam3.fa $dataN/SAMEA4672967_2.fastq.gz > $result/SAMEA4672967_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SAMEA4672967_1.sai \
$result/SAMEA4672967_2.sai \
$dataN/SAMEA4672967_1.fastq.gz \
$dataN/SAMEA4672967_2.fastq.gz \
> $result/SAMEA4672967.sam




####### Convert sam to bam file #######
samtools view -bS $result/SAMEA4672968.sam > $result/SAMEA4672968.bam
samtools view -bS $result/SAMEA4672967.sam > $result/SAMEA4672967.bam




####### Delete #######
#$dataN/*.fastq.gz $dataT/*.fastq.gz
rm $result/*.sai




####### Picard #######
# get header information
grep '@SQ\|@PG' $result/SAMEA4672967.sam > $result/header_N
grep '@SQ\|@PG' $result/SAMEA4672968.sam > $result/header_T

# exclude unmapped reads based on FLAG
cat $result/SAMEA4672967.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/SAMEA4672967-cleaned.sam
cat $result/header_N $result/SAMEA4672967-cleaned.sam > $result/fooN
mv $result/fooN $result/SAMEA4672967-cleaned.sam
cat $result/SAMEA4672968.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/SAMEA4672968-cleaned.sam
cat $result/header_T $result/SAMEA4672968-cleaned.sam > $result/fooT
mv $result/fooT $result/SAMEA4672968-cleaned.sam

# picard sort
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/SAMEA4672967-cleaned.sam O=$result/SAMEA4672967_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/SAMEA4672968-cleaned.sam O=$result/SAMEA4672968_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

# picard MarkDuplicates
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/SAMEA4672967_rg_added_sorted.bam O=$result/SAMEA4672967_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$result/SAMEA4672967-output.metrics REMOVE_SEQUENCING_DUPLICATES=true
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/SAMEA4672968_rg_added_sorted.bam O=$result/SAMEA4672968_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$result/SAMEA4672968-output.metrics REMOVE_SEQUENCING_DUPLICATES=true





####### Remove unneeded files #######
rm $result/SAMEA4672967.sam $result/SAMEA4672968.sam 
rm $result/header_N $result/header_T
rm $result/SAMEA4672967_rg_added_sorted.bam $result/SAMEA4672968_rg_added_sorted.bam
rm $result/SAMEA4672967-cleaned.sam $result/SAMEA4672968-cleaned.sam





####### GATK Realign #######
# Generating interval file for sort.bam
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa -I $result/SAMEA4672967_rg_added_sorted_dedupped_removed.bam -o $result/SAMEA4672967_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa -I $result/SAMEA4672968_rg_added_sorted_dedupped_removed.bam -o $result/SAMEA4672968_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

# realign
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa -I $result/SAMEA4672967_rg_added_sorted_dedupped_removed.bam -targetIntervals $result/SAMEA4672967_rg_added_sorted_dedupped_removed.bam.intervals -o $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam --allow_potentially_misencoded_quality_scores
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa -I $result/SAMEA4672968_rg_added_sorted_dedupped_removed.bam -targetIntervals $result/SAMEA4672968_rg_added_sorted_dedupped_removed.bam.intervals -o $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam --allow_potentially_misencoded_quality_scores



cd $result
######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 
module load MuTect/1.1.7-Java-1.7.0_80

time java -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence $reference/canFam3.fa \
--dbsnp $reference/dbsnp_9615.vcf \
--intervals $reference/Canis_familiaris.CanFam3.interval_list  \
--input_file:normal $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam \
--input_file:tumor $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam \
--out $result/LAB3_dam_SAMEA4672968_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file $result/LAB3_dam_SAMEA4672968_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf $result/LAB3_dam_SAMEA4672968_rg_added_sorted_dedupped_removed.MuTect.vcf

####### Mutect2 #######
#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T MuTect2 \
#-R $reference/canFam3.fa \
#-I:tumor $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam \
#-I:normal $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam \
#--dbsnp $reference/dbsnp_9615.vcf \
#-L $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
#-o $result/LAB3_dam_SAMEA4672968_rg_added_sorted_dedupped_removed.MuTect.vcf


module load SAMtools/1.9-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load GATK/3.8-1-Java-1.8.0_144
module load Perl/5.26.1-GCCcore-6.4.0

cd $result

####### Annovar #######
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $result/LAB3_dam_SAMEA4672968_rg_added_sorted_dedupped_removed.MuTect.vcf > $result/LAB3_dam_SAMEA4672968_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $reference/annovar/convert2annovar.pl -format vcf4old $result/LAB3_dam_SAMEA4672968_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > $result/LAB3_dam_SAMEA4672968_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $reference/annovar/annotate_variation.pl --buildver canFam3 $result/LAB3_dam_SAMEA4672968_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $reference






####### remove unneeded files #######
rm $result/*_sorted_dedupped_removed.bam
#rm $result/LAB3_dam_SAMEA4672968_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS $result/LAB3_dam_SAMEA4672968_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
#rm $result/*.bai
"

module load SAMtools/1.9-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load GATK/3.8-1-Java-1.8.0_144
module load Perl/5.26.1-GCCcore-6.4.0

data='/scratch/kh31516/trio-family-temp'
result='/scratch/kh31516/trio-family-temp'
cd $result

################ GATK with coverage condition ###################
#scp $data/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam
#samtools index $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam
#samtools index $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam
# Variant calling
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference/canFam3.fa -I $data/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam.vcf -bamout $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bamout.bam
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference/canFam3.fa -I $data/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam.vcf -bamout $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bamout.bam

# Variant filtering
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $reference/canFam3.fa -V $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam.vcf -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $reference/canFam3.fa -V $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam.vcf -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf

cd $result
################ Annovar ###################
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf > $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS
awk '$7 == "PASS" {print $0}' $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf > $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS
# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS > $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS-avinput
perl $annovar_index/convert2annovar.pl -format vcf4old $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS > $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS-avinput $annovar_index
perl $annovar_index/annotate_variation.pl --buildver canFam3 $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS-avinput $annovar_index

# add gene names to annovar output
python $script/Add_GeneName_N_Signature.py $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS-avinput.exonic_variant_function $reference/Canis_familiaris.CanFam3.1.96.gtf_geneNamePair.txt
python $script/Add_GeneName_N_Signature.py $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS-avinput.exonic_variant_function $reference/Canis_familiaris.CanFam3.1.96.gtf_geneNamePair.txt

############### Delete files ##############
#rm $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS $result/SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS-avinput
#rm $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS $result/SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam.filter.vcf-PASS-avinput


## GATK_callable ##
#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
# -R $reference/canFam3.fa \
# -I SAMEA4672967_rg_added_sorted_dedupped_removed.realigned.bamout.bam \
# -summary SAMEA4672967_table.txt \
# -o SAMEA4672967_callable_status.bed 

#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
# -R $reference/canFam3.fa \
# -I SAMEA4672968_rg_added_sorted_dedupped_removed.realigned.bamout.bam \
# -summary SAMEA4672968_table.txt \
# -o SAMEA4672968_callable_status.bed 