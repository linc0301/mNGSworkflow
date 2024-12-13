#!/bin/bash
# Description: this is the bowtie2 optimization of mNGS workflow
# Version: 1.0
# Author: LinC
WORKDIR=~/workspace/mNGS/data/bowtie2optimize
REFDIR=~/workspace/mNGS/data/reference

i=3

#First round of bowtie
 # pairend alignment --very-sensitive-local 
    bowtie2 --very-sensitive --no-unal -x ${REFDIR}/human/bt2/human_host \
        --threads 24 \
        -1 ~/workspace/mNGS/data/simdata_zymo/sim.clean.R1_${i}.fastq.gz \
        -2 ~/workspace/mNGS/data/simdata_zymo/sim.clean.R2_${i}.fastq.gz \
        --un-conc-gz ${WORKDIR}/zymo_rm_host_${i} \
        -S ${WORKDIR}/zymo_rm_host_${i}.sam
        echo "iteration $i unaligned reads (which should contain non-host sequences) concatenated into a single file and gzipped."
    
    # SAM to BAM 
    samtools view -@ 10 -S -b -o ${WORKDIR}/zymo_rm_host_${i}.bam ${WORKDIR}/zymo_rm_host_${i}.sam
    
    mv ${WORKDIR}/zymo_rm_host_${i}.1 ${WORKDIR}/zymo_rm_host_${i}.R1.fastq.gz
    mv ${WORKDIR}/zymo_rm_host_${i}.2 ${WORKDIR}/zymo_rm_host_${i}.R2.fastq.gz



    bowtie2 --very-sensitive --no-unal -x ${REFDIR}/rRNA/bt2/rRNA \
    --threads 10 \
    -1 ${WORKDIR}/zymo_rm_host_${i}.R1.fastq.gz \
    -2 ${WORKDIR}/zymo_rm_host_${i}.R2.fastq.gz \
    --un-conc-gz ${WORKDIR}/zymo_rm_rRNA_${i} \
    -S ${WORKDIR}/zymo_rm_rRNA_${i}.sam
    echo "iteration $i unaligned reads (which should contain non-rRNA sequences) concatenated into a single file and gzipped."
    
    samtools view -@ 10 -S -b -o ${WORKDIR}/zymo_rm_rRNA_${i}.bam ${WORKDIR}/zymo_rm_rRNA_${i}.sam
    
    mv ${WORKDIR}/zymo_rm_rRNA_${i}.1 ${WORKDIR}/zymo_rm_rRNA_${i}.R1.fastq.gz
    mv ${WORKDIR}/zymo_rm_rRNA_${i}.2 ${WORKDIR}/zymo_rm_rRNA_${i}.R2.fastq.gz

#Second Round 

 bowtie2 --no-unal -x ${REFDIR}/human/bt2/human_host \
    --threads 10 \
    -1 ${WORKDIR}/zymo_rm_rRNA_${i}.R1.fastq.gz \
    -2 ${WORKDIR}/zymo_rm_rRNA_${i}.R2.fastq.gz \
    --un-conc-gz ${WORKDIR}/zymo_rm_host_2_${i} \
    -S ${WORKDIR}/zymo_rm_host_2_${i}.sam
    echo "iteration $i unaligned reads (which should contain non-human sequences) concatenated into a single file and gzipped."



# samtools + bowtie2 


bowtie2 -x ${REFDIR}/human/bt2/human_host \
    --threads 24 \
    -1 ~/workspace/mNGS/data/simdata_zymo/sim.clean.R1_${i}.fastq.gz \
    -2 ~/workspace/mNGS/data/simdata_zymo/sim.clean.R2_${i}.fastq.gz \
    -S ${WORKDIR}/zymo_rm_host_bam_${i}.sam
samtools view -@ 10 -S -b -o ${WORKDIR}/zymo_rm_host_bam_${i}.bam ${WORKDIR}/zymo_rm_host_bam_${i}.sam
      


samtools sort -@ 8 zymo_rm_host_bam_${i}.bam -o zymo_rm_host_bam_${i}.sorted.bam
samtools index zymo_rm_host_bam_${i}.sorted.bam
samtools view -b -f 4 \
    zymo_rm_host_bam_${i}.sorted.bam \
    > zymo_rm_host_bam_${i}.unmapped.bam
    
samtools sort -n zymo_rm_host_bam_${i}.unmapped.bam \
    -O BAM \
    -o zymo_rm_host_bam_${i}.unmapped.sorted.bam
    
samtools fastq \
    -@ 8 \
    zymo_rm_host_bam_${i}.unmapped.sorted.bam \
    -1 zymo_rm_host_bam_${i}.R1.fastq.gz \
    -2 zymo_rm_host_bam_${i}.R2.fastq.gz \
    -n