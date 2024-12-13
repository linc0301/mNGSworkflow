#!/bin/bash
# Description: this is the bowtie2 optimized host removal script
# Version: 1.0
# Author: LinC

show_help() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -h, --help                Show this help message and exit"
    echo "  -R, --refdir <dir>        Specify the reference directory (REFDIR)"
    echo "  -W, --workdir <dir>       Specify the working directory (WORKDIR)"
    echo "  -D, --db <db>             Specify the Kraken2 database directory (DB)"
    echo "  -i, --input <prefix>      Specify the input file prefix (inFile)"
    echo "  -o, --output <prefix>     Specify the output file prefix (outFile)"
    echo "  -p, --paired_end          Specify if the input is paired-end reads"
    exit 0
}

error_handler() {
    echo "An error occurred. Exiting..."
    exit 1
}

trap 'error_handler' ERR

REFDIR=~/workspace/mNGS/data/reference
WORKDIR=~/workspace/mNGS/data
DB=k2_standard_08gb
inFile=""
outFile=""
paired_end=false  # Default to single-end

while getopts "hR:W:D:i:o:p" opt; do
  case $opt in
    h)  show_help               ;;
    R)  REFDIR=$OPTARG          ;;
    W)  WORKDIR=$OPTARG         ;;
    D)  DB=$OPTARG              ;;
    i)  inFile=$OPTARG          ;;
    o)  outFile=$OPTARG         ;;
    p)  paired_end=true        ;;
    *)  show_help               ;;
  esac
done

if [ -z "$REFDIR" ] || [ -z "$WORKDIR" ] || [ -z "$DB" ] || [ -z "$inFile" ] || [ -z "$outFile" ]; then
    echo "Error: Missing arguments."
    show_help
fi

echo "Input file is set to $inFile"



if [ "$paired_end" = false ]; then
    echo "Processing single-end (SE) reads."
    
    # Bowtie2 alignment for single-end reads
    bowtie2 --very-sensitive --no-unal -x ${REFDIR}/human/bt2/human_host \
        --threads 24 \
        -U ${inFile}.fastq.gz \
        --un ${WORKDIR}/bowtie_rm_unaligned.fastq.gz \
        -S ${WORKDIR}/bowtie_rm.sam

    echo "${inFile} unaligned reads (which should contain non-host sequences) concatenated into a single file and gzipped by bowtie2."

    # Bowtie2 alignment against rRNA for unaligned reads from previous step
    bowtie2 --very-sensitive --no-unal -x ${REFDIR}/rRNA/bt2/rRNA \
        --threads 10 \
        -U ${WORKDIR}/bowtie_rm_unaligned.fastq.gz \
        --un ${WORKDIR}/bowtie_rm_rRNA_unaligned.fastq.gz \
        -S ${WORKDIR}/bowtie_rm_rRNA.sam

    echo "${inFile} unaligned reads (which should contain non-rRNA sequences) concatenated into a single file and gzipped by bowtie2."

    # Samtools removal
    bowtie2 -x ${REFDIR}/human/bt2/human_host \
        --threads 24 \
        -U ${WORKDIR}/bowtie_rm_rRNA_unaligned.fastq.gz \
        -S ${WORKDIR}/sam_rm_host.sam

    samtools view -@ 10 -S -b -o ${WORKDIR}/sam_rm_host.bam ${WORKDIR}/sam_rm_host.sam

    samtools sort -@ 8 ${WORKDIR}/sam_rm_host.bam -o ${WORKDIR}/sam_rm_host.sorted.bam
    samtools index ${WORKDIR}/sam_rm_host.sorted.bam

    # Extract unmapped reads
    samtools view -b -f 4 ${WORKDIR}/sam_rm_host.sorted.bam > ${WORKDIR}/unmapped.bam

    # Sort unmapped reads
    samtools sort -n unmapped.bam -O BAM -o ${WORKDIR}/unmapped.sorted.bam

    # Convert sorted unmapped BAM to FASTQ
    samtools fastq \
    -@ 8 \
    ${WORKDIR}/unmapped.sorted.bam > ${WORKDIR}/${outFile}.fastq.gz



else
    echo "Processing paired-end (PE) reads."
    
    # pairend alignment --very-sensitive-local 
    bowtie2 --very-sensitive --no-unal -x ${REFDIR}/human/bt2/human_host \
        --threads 24 \
        -1 ${inFile}.R1.fastq.gz \
        -2 ${inFile}.R2.fastq.gz \
        --un-conc-gz ${WORKDIR}/bowtie_rm \
        -S ${WORKDIR}/bowtie_rm.sam
    
    echo "${inFile} unaligned reads (which should contain non-host sequences) concatenated into a single file and gzipped by bowtie2."

    bowtie2 --very-sensitive --no-unal -x ${REFDIR}/rRNA/bt2/rRNA \
        --threads 10 \
        -1 ${WORKDIR}/bowtie_rm.1 \
        -2 ${WORKDIR}/bowtie_rm.2 \
        --un-conc-gz ${WORKDIR}/bowtie_rm_rRNA \
        -S ${WORKDIR}/bowtie_rm_rRNA.sam
 echo "${inFile} unaligned reads (which should contain non-rRNA sequences) concatenated into a single file and gzipped by bowtie2."
    
#Samtools removal
    bowtie2 -x ${REFDIR}/human/bt2/human_host \
        --threads 24 \
        -1 ${WORKDIR}/bowtie_rm_rRNA.1 \
        -2 ${WORKDIR}/bowtie_rm_rRNA.2 \
        -S ${WORKDIR}/sam_rm_host.sam
    samtools view -@ 10 -S -b -o ${WORKDIR}/sam_rm_host.bam ${WORKDIR}/sam_rm_host.sam
      


    samtools sort -@ 8 ${WORKDIR}/sam_rm_host.bam -o ${WORKDIR}/sam_rm_host.sorted.bam
    samtools index ${WORKDIR}/sam_rm_host.sorted.bam
    samtools view -b -f 4 \
        ${WORKDIR}/sam_rm_host.sorted.bam \
        > ${WORKDIR}/unmapped.bam
    
    samtools sort -n unmapped.bam \
        -O BAM \
        -o ${WORKDIR}/unmapped.sorted.bam
    
    samtools fastq \
        -@ 8 \
        unmapped.sorted.bam \
        -1 ${WORKDIR}/${outFile}.R1.fastq.gz \
        -2 ${WORKDIR}/${outFile}.R2.fastq.gz \
        -n

fi

rm -f bowtie_rm* sam_rm* unmapped*
