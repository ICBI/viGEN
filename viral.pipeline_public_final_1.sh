#!/bin/bash

# Alignment using Bowtie2
functionBowtie2() {
	echo "Start - Align using Bowtie2 local algin" >> $LOG
	$BOWTIE2_SW"bowtie2" -x $BOWTIE2_REF --local -1 $I_FOLDER$ID"_1.fq" -2 $I_FOLDER$ID"_2.fq" --threads $N_THREAD --al $O_FOLDER$ID/$ID".bowtie2.align" --un $O_FOLDER$ID/$ID".Bowtie2.unaligned" -S $O_FOLDER$ID/$ID".bowtie2.sam" --sensitive > $O_FOLDER$ID/$ERR2 2> $O_FOLDER$ID/$ERR3
	echo "End - Align using Bowtie2 local algin" >> $LOG
}

#Convert SAM file to BAM file
functionConvertSamtoBam() {
    echo "Start - convert sam to bam file using Samtools" >> $LOG
	$SAMTOOLS"samtools" view -bS $O_FOLDER$ID/$ID".bowtie2.sam" > $O_FOLDER$ID/$ID".bowtie2.bam"
    echo "End - convert sam to bam file using Samtools" >> $LOG
}

#Sort BAM coordinate wise
functionSortBam() {
    echo "Start - Samtools sort bam file" >> $LOG
	$SAMTOOLS"samtools" sort $O_FOLDER$ID/$ID".bowtie2.bam" -o $O_FOLDER$ID/$ID".bowtie2.sorted.bam" -@ $N_THREAD
    echo "End - Samtools sort bam file" >> $LOG
} 

#Index BAM file
functionIndexBam() {
    echo "Start - Samtools index bam file" >> $LOG
	$SAMTOOLS"samtools" index $O_FOLDER$ID/$ID".bowtie2.sorted.bam"
    echo "End - Samtools index bam file" >> $LOG
}

# Samtools idx -  genome level counts
functionIdxStatsBam() {
    echo "Start - Samtools Idx" >> $LOG
	$SAMTOOLS"samtools" idxstats $O_FOLDER$ID/$ID".bowtie2.sorted.bam" > $O_FOLDER$ID/$ID".bowtie2.idxstats.txt"
    echo "End - Samtools Idx" >> $LOG
 }

########################################
# Values to set by the user
N_THREAD=4
I_FOLDER="/Users/ls483/Documents/SRA.GEO/input/unmapped_fastq/"
O_FOLDER="/Users/ls483/Documents/SRA.GEO/output_bowtie2/"

LOG=$O_FOLDER"/pipeline_log.txt"
ERR2="system.err.bowtie2.txt"
ERR3="log.bowtie2.txt"
SAMTOOLS="/Users/ls483/Documents/software/samtools/bin/"
BOWTIE2_SW="/Users/ls483/Documents/software/bowtie/bowtie2-2.2.9/"
BOWTIE2_REF="/Users/ls483/Documents/SRA.GEO/viral.reference/virus.bowtie2.refB/virus.bowtie2.refB"
		
for x in $(cat samplenames.txt)
do
    echo $x
    ID=${x%*}

    echo "-----------------------------" >> $LOG
    echo $ID >> $LOG

    START=$(date +%s) 
    
    mkdir $O_FOLDER$ID

    functionBowtie2

    functionConvertSamtoBam

    functionSortBam

    functionIndexBam

    functionIdxStatsBam

    functionFlagstatBam
    
	END=$(date +%s)
	DIFF=$(( ($END - $START)/60 )) 
	echo "It took $DIFF minutes" >> $LOG
done 
