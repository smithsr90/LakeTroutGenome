
###check the best library to deep sequencing

module load mugqic/java/openjdk-jdk1.8.0_72
module load mugqic/picard/2.9.0
module load mugqic/python/2.7.14
module load mugqic/genpipes/3.1.4
module load mugqic/hicup/v0.7.0

#creating the digester file

hicup_digester --genome Lake_Trout --re1 N^GATC,Arima:G^ANTC,Arima LT/LT_pb-polished.contigset.fasta


#create the dictionary file


#replace | in fasta file with _
sed 's/|/_/g' LT_pb-polished.contigset.fasta.ori > LT_pb-polished.contigset.fasta

java -jar $PICARD_HOME/picard.jar CreateSequenceDictionary R=LT/LT_pb-polished.contigset.fasta O=LT_assembly.dict

bowtie2 index

module load mugqic/bowtie2/2.3.1
bowtie2-build LT/LT_pb-polished.contigset.fasta LT_bowtie2_index

#create a sim link for genome fasta file (file name should be match to the bowtie2 indexes)
cd bowtie2_index
ln -s ../LT/LT_pb-polished.contigset.fasta LT_bowtie2_index.fa


#run genpipes

genpipes/pipelines/hicseq/hicseq.py -c genpipes/pipelines/hicseq/hicseq.base.ini custom.ini -r readset_LT_janick_hicseq.txt -s 1-19 -e Arima > hicseqScript.txt
 bash hicseqScript.txt
 
 #checked the file Analysis_Summary_Report_4.html and chose the library to deep sequencing
 
 
 ##Arima mapping pipeline for mapping Hi-C reads with draft contig genome
 ####################################################
 #####scaffolding with deep sequenced library
 
 #! /bin/bash
#PBS -q sw
#PBS -l nodes=1:ppn=12
#PBS -l walltime=2:00:00:00
#PBS -o Arima.out
#PBS -e Arima.err
#PBS -N Arima_pipeline
#PBS -V

##############################################
# ARIMA GENOMICS MAPPING PIPELINE 02/08/2019 #
##############################################

#Below find the commands used to map HiC data.

#Replace the variables at the top with the correct paths for the locations of files/programs on your system.

#This bash script will map one paired end HiC dataset (read1 & read2 fastqs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Commands #
##########################################

 module load mugqic/samtools/1.9
 module load mugqic/bwa/0.7.17
 module load mugqic/java/openjdk-jdk1.8.0_72
 module load mugqic/picard/2.17.3

SRA='LT_muscle_ArimaLucigen_BER806A4-1_S1_L004_'
LABEL='LT_muscle_ArimaLucigen'
BWA='bwa'
SAMTOOLS='samtools'
IN_DIR='scafolding/raw_file'
REF='LT/LT_pb-polished.contigset.fasta'
FAIDX='$REF.fai'
PREFIX='LT_pb-polished.contigset.fasta'
RAW_DIR='scafolding/Arima_pipeline/bam_out'
FILT_DIR='scafolding/Arima_pipeline/filtered'
FILTER='scafolding/Arima_pipeline/mapping_pipeline/filter_five_end.pl'
COMBINER='scafolding/Arima_pipeline/mapping_pipeline/two_read_bam_combiner.pl'
STATS='scafolding/Arima_pipeline/mapping_pipeline/get_stats.pl'
PICARD=$PICARD_HOME'/picard.jar'
TMP_DIR='scafolding/Arima_pipeline/temp'
PAIR_DIR='scafolding/Arima_pipeline/paired'
REP_DIR='scafolding/Arima_pipeline/duplicated'
REP_LABEL=$LABEL\_rep1
MERGE_DIR='scafolding/Arima_pipeline/merged'
MAPQ_FILTER=10
CPU=12

echo "### Step 0: Check output directories exist & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
#$BWA index -a bwtsw -p $PREFIX $REF

echo "### Step 1.A: FASTQ to BAM (1st)"
#$BWA mem -t $CPU $PREFIX $IN_DIR/$SRA\R1_001.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
#$BWA mem -t $CPU $PREFIX $IN_DIR/$SRA\R2_001.fastq.gz  | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
#$SAMTOOLS view -h $RAW_DIR/$SRA\_1.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA\_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
#$SAMTOOLS view -h $RAW_DIR/$SRA\_2.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA\_2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
#perl $COMBINER $FILT_DIR/$SRA\_1.bam $FILT_DIR/$SRA\_2.bam $SAMTOOLS $MAPQ_FILTER | $SAMTOOLS view -bS -t $FAIDX - | $SAMTOOLS sort -@ $CPU -o $TMP_DIR/$SRA.bam -

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/$SRA.bam OUTPUT=$PAIR_DIR/$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none

###############################################################################################################################################################
###                                           How to Accommodate Technical Replicates                                                                       ###
### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                                                    ###
### Technical replicates (eg. one library split across multiple lanes) should be merged before running the MarkDuplicates command.                          ###
### If this step is run, the names and locations of input files to subsequent steps will need to be modified in order for subsequent steps to run correctly.###
### The code below is an example of how to merge technical replicates.                                                                                      ###
###############################################################################################################################################################
#	REP_NUM=X #number of the technical replicate set e.g. 1
#	REP_LABEL=$LABEL\_rep$REP_NUM
#	INPUTS_TECH_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as technical replicates
#   example bash array - INPUTS_TECH_REPS=('INPUT=A.L1.bam' 'INPUT=A.L2.bam' 'INPUT=A.L3.bam')
#	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_TECH_REPS OUTPUT=$TMP_DIR/$REP_LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

$SAMTOOLS index $REP_DIR/$REP_LABEL.bam

perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"


 ###QC stats
cd scafolding/Arima_pipeline

module load mugqic/samtools/1.9
samtools sort -n duplicated_eric/LT_muscle_ArimaLucigen_rep1.bam  > duplicated_eric/LT_muscle_ArimaLucigen_rep1.sort.bam

#download https://github.com/phasegenomics/hic_qc
python hic_qc/hic_qc/hic_qc.py -b duplicated_eric/LT_muscle_ArimaLucigen_rep1.sort.bam -r -o qc_stats_eric/LT_muscle -n 5000000



#############################################################################################
#############################################################################################
##SALSA scaffolding

#SALSA

cd scafolding/SALSA
module load mugqic/bedtools/2.27.0
module load mugqic/samtools/1.9


module load mugqic/bedtools/2.27.0
bamToBed -i ../duplicated/LT_muscle_ArimaLucigen_rep1.bam > LT_muscle_ArimaLucigen_rep1.bam.bed

sort -k 4 LT_muscle_ArimaLucigen_rep1.bam.bed > LT_muscle_ArimaLucigen_rep1.sort.bed

module load mugqic/samtools/1.9
samtools faidx LT_pb-polished.contigset.fasta


python SALSA/SALSA-2.2/run_pipeline.py -a LT_pb-polished.contigset.fasta -l LT_pb-polished.contigset.fasta.fai -b LT_muscle_ArimaLucigen_rep1.sort.bed -e GATC,GANTC -o scaffolds -m yes

####################################
#new polished ref genome

module load mugqic/bedtools/2.27.0
bamToBed -i ../duplicated_eric/LT_muscle_ArimaLucigen_rep1.bam > LT_muscle_ArimaLucigen_rep1.eric.bam.bed

sort -k 4 LT_muscle_ArimaLucigen_rep1.eric.bam.bed > LT_muscle_ArimaLucigen_rep1.sort.eric.bed

module load mugqic/samtools/1.9
samtools faidx lake_trout_after_pilon.fasta


python SALSA/SALSA-2.2/run_pipeline.py -a lake_trout_after_pilon.fasta -l lake_trout_after_pilon.fasta.fai -b LT_muscle_ArimaLucigen_rep1.sort.eric.bed -e GATC,GANTC -o scaffolds -m yes





###quast quality


module purge
module load mugqic_dev/quast/5.0.2

quast.py scaffolds/scaffolds_FINAL.fasta --output-dir metrics/ --threads 1
quast.py LT_pb-polished.contigset.fasta --output-dir metrics/ --threads 12
quast.py scaffolds_eric/scaffolds_FINAL.fasta --output-dir metrics_eric/ --threads 3

quast.py scaffolds_eric/scaffolds_FINAL.fasta scaffolds/scaffolds_FINAL.fasta LT_pb-polished.contigset.fasta lake_trout_after_pilon.fasta --output-dir metrics_all/ --threads 3

quast.py *.fasta --output-dir metrics_all/ --threads 3

quast.py scaffolds_eric_i3/scaffolds_FINAL.fasta lake_trout_after_pilon.fasta --output-dir metrics_manuscript/ --threads 2

quast.py *.fasta scaffolds*/*FINAL.fasta --output-dir metrics_all_new/ --threads 1





#Generate a hic file for visualization

python ${SCRIPT_PATH}/alignments2txt.py -b ${SALSA_OUT_DIR}/alignment_iteration_1.bed  -a ${SALSA_OUT_DIR}/scaffolds_FINAL.agp -l ${SALSA_OUT_DIR}/scaffold_length_iteration_1 > ${SALSA_OUT_DIR}/alignments.txt

awk '{if ($2 > $6) {print $1"\t"$6"\t"$7"\t"$8"\t"$5"\t"$2"\t"$3"\t"$4} else {print}}' ${SALSA_OUT_DIR}/alignments.txt | sort -k2,2d -k6,6d --parallel=16 | awk 'NF'  > ${SALSA_OUT_DIR}/alignments_sorted.txt

java -jar ${JUICER_JAR} pre ${SALSA_OUT_DIR}/alignments_sorted.txt ${SALSA_OUT_DIR}/salsa_scaffolds.hic ${SALSA_OUT_DIR}/chromosome_sizes.tsv


module load mugqic/samtools/1.9
samtools sort -n duplicated/LT_muscle_ArimaLucigen_rep1.bam > duplicated/LT_muscle_ArimaLucigen_rep1.sort.bam
module load mugqic/python/2.7.14
python hic_qc/hic_qc.py -b duplicated/*.sort.bam -r -o qc_stats/LT_muscle
```

