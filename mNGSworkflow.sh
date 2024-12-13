#!/bin/bash
# Description: this is the workflow of in silico simulation for mNGS sequencing data
# Version: 1.0
# Author: LinC

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $*" >> "$LOGFILE"
}

error_handling_function() {
    log "An error occurred at line $1."
    cleanup
    exit 1
}

cleanup() {
    log "Cleaning up resources."
    # Additional cleanup code here
}


WORKDIR=~/workspace/mNGS/data
REFDIR=${WORKDIR}/reference
LOGFILE=${WORKDIR}/log/mNGSworkflow.log


cd $WORKDIR

log "Script Started"

# data fetch

##############################
## hg19
# wget -c ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
# tar -zxvf chromFa.tar.gz
# cat chr*.fa > hg19.fa
##############################


## hg38.p14 which matches the RNA fna file
humanref=hg38.p14.fa.gz
if [ ! -f $humanref ]; then
    wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p14/hg38.p14.fa.gz
fi


## zymo
###############################
#prefetch SRR27174691
#fastq-dump --gzip --split-files SRR27174691.sra
###############################

## in silico seq simulation using human genome
## generate 3 groups of sim reads
## microbial: human = 10%

for i in {1..3}
do
    echo "Iteration $i"

    #  wgsim
    rand_num=$(date +%N)
    seed=$rand_num
    wgsim -d 400 -N 18000000 -1 100 -2 100 -S $seed ${humanref} human.R1.fq human.R2.fq

    # bacteria wgsim
    for bacteria in `cat ${WORKDIR}/sra/bacteria.list`; do
    wgsim -e 0.01 -d 400 -s 40 -1 50 -2 50 -N 200000 -S $seed ${WORKDIR}/sra/${bacteria}_*fna \
    ${WORKDIR}/simdata_bacteria/${bacteria}.sim.R1.fq \
    ${WORKDIR}/simdata_bacteria/${bacteria}.sim.R2.fq
    done

    # zymo 
    zymo=${WORKDIR}/sra/SRR27174691/SRR27174691

    # seqtk
    seqtk sample -s $seed ${zymo}_1.fastq.gz 2000000 > zymo.sample.R1.fq
    seqtk sample -s $seed ${zymo}_2.fastq.gz 2000000 > zymo.sample.R2.fq

    # merge
    cat human.R1.fq zymo.sample.R1.fq > merge_${i}.R1.fastq
    cat human.R2.fq zymo.sample.R2.fq > merge_${i}.R2.fastq

    # set QC
    python code/setQC.py -i merge.R1_${i}.fastq -o sim_${i}.R1.fastq
    python code/setQC.py -i merge.R2_${i}.fastq -o sim_${i}.R2.fastq

    # clean up

    rm human.R1.fq human.R2.fq zymo.sample.R1.fq zymo.sample.R2.fq merge_${i}.R1.fastq merge_${i}.R2.fastq

done

# Quality Ctrl
## fastp 

for i in {1..3}
do
    echo "Processing iteration $i"

    input_R1=${WORKDIR}/sim_${i}.R1.fastq
    input_R2=${WORKDIR}/sim_${i}.R2.fastq
    output_R1_gzip=${WORKDIR}/sim.clean_${i}.R1.fastq.gz
    output_R2_gzip=${WORKDIR}/sim.clean_${i}.R2.fastq.gz
    json_report=${WORKDIR}/sim_${i}.json
    html_report=${WORKDIR}/sim_${i}.html


    fastp -i $input_R1 \
          -I $input_R2 \
          -o $output_R1_gzip \
          -O $output_R2_gzip \
          -w 10 \
          --disable_quality_filtering \
          --cut_right \
          --cut_window_size 4 \
          --cut_mean_quality 20 \
          --cut_front \
          --cut_front_window_size 1 \
          --cut_front_mean_quality 20 \
          --cut_tail \
          --cut_tail_window_size 1 \
          --cut_tail_mean_quality 20 \
          --length_required 30 \
          --dont_eval_duplication \
          -j $json_report \
          -h $html_report
    # check if fastp worked
    if [ $? -eq 0 ]; then
        echo "fastp successfully processed iteration $i"
    else
        echo "fastp failed to process iteration $i"
    fi
done


##SE fastp###

fastp -i $input_SE \
      -o $output_SE_gzip \
      -w 10 \
      --disable_quality_filtering \
      --cut_right \
      --cut_window_size 4 \
      --cut_mean_quality 20 \
      --cut_front \
      --cut_front_window_size 1 \
      --cut_front_mean_quality 20 \
      --cut_tail \
      --cut_tail_window_size 1 \
      --cut_tail_mean_quality 20 \
      --length_required 30 \
      --dont_eval_duplication \
      -j $json_report \
      -h $html_report


# remove the host genome
# conda activate bowtieenv
bowtie2-build --threads 30 ${WORKDIR}/${humanref} ${REFDIR}/human/bt2/human_host >${REFDIR}/human/bt2/bt2-build.log 2>&1 &

for i in {1..3}; do
    # pairend alignment --very-sensitive-local 
    bowtie2 --very-sensitive --no-unal -x ${REFDIR}/human/bt2/human_host \
        --threads 10 \
        -1 ${WORKDIR}/sim.clean.R1_${i}.fastq.gz \
        -2 ${WORKDIR}/sim.clean.R2_${i}.fastq.gz \
        --un-conc-gz ${WORKDIR}/zymo_rm_host_${i} \
        -S ${WORKDIR}/zymo_rm_host_${i}.sam
        echo "iteration $i unaligned reads (which should contain non-host sequences) concatenated into a single file and gzipped."
    
    # SAM to BAM 
    samtools view -@ 10 -S -b -o ${WORKDIR}/zymo_rm_host_${i}.bam ${WORKDIR}/zymo_rm_host_${i}.sam
    
    mv ${WORKDIR}/zymo_rm_host_${i}.1 ${WORKDIR}/zymo_rm_host_${i}.R1.fastq.gz
    mv ${WORKDIR}/zymo_rm_host_${i}.2 ${WORKDIR}/zymo_rm_host_${i}.R2.fastq.gz
done

#################################
#unzipped_humanref=${humanref%.gz}
#gzip -d $humanref > $unzipped_humanref
# snap-aligner index ${WORKDIR}/${unzipped_humanref} ${REFDIR}/human/snap/ -locationSize 5 1>${REFDIR}/human/snap/snap-build.log 2>&1 &

# snap-aligner paired ${REFDIR}/human/snap ${WORKDIR}/sim.clean.R1_${i}.fastq.gz ${WORKDIR}/sim.clean.R2_${i}.fastq.gz -o ${WORKDIR}/snap.sam
#################################


# RNA removal

zgrep '^>' GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna | grep "gbkey=rRNA" | awk '{print $1}' | sed 's/>//g' > id.list
seqtk subseq GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz id.list > rRNA.fa

### cluster
cd-hit-est -i rRNA.fa -o rRNA.new.fa -c 1 -n 10 -d 0 -M 20000 -T 10
mv rRNA.new.fa rRNA.fa
rm rRNA.new.fa.*

### 
bowtie2-build --threads 30 ~/workspace/mNGS/data/rRNA.fa ${REFDIR}/rRNA/bt2/rRNA
for i in {1..3}; do
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

done
## conda deactivate

## one step host rm for SE50
for i in `ls simdata_bacteria/GCF*R1.fastq.gz`
do
    echo $i
    bash $WORKDIR/code/host_rm.sh -i ${i%.fastq.gz} -o ${i%.fastq.gz}.host_rm

done


repFile=${WORKDIR}/simdata_bacteria/report
DB=k2_standard_db 
for i in `ls simdata_bacteria/GCF*host_rm*`
do
    i=$(basename $i)

    conda activate base
    kraken2 --db ${DB} --threads 8 ${WORKDIR}/simdata_bacteria/${i} \
    --output ${repFile}/${i%.host*}_output.kraken --report ${repFile}/${i%.host*}.kreport
    conda activate bowtieenv
    bracken -i  ${repFile}/${i%.host*}.kreport -d $DB -o  ${repFile}/${i%.host*}.kreport.bracken -w  ${repFile}/${i%.host*}.bracken.report -t 8
    conda deactivate

done







repFile=${WORKDIR}/simdata_bacteria/report
DB=kraken_std_db
    conda activate base
    kraken2 --db ${DB} --threads 8 ${WORKDIR}/simdata_bacteria/${i} \
    --output ${repFile}/${i%.host*}_output.kraken --report ${repFile}/${i%.host*}.kreport










### Kraken2
#################################
##error with rsync so a mannual building process has been taken
##Edited the raw script 
#################################

##conda activate bowtieenv
##kraken2-build --standard --threads 4 --db ${WORKDIR}/kraken2_db
##conda deactivate


DB=k2_standard_08gb
#################################
##run kraken2&bracken in different env
#################################

### Kraken2 analysis pure zymo###
#conda activate base
#kraken2 --db ${DB} --threads 8 --paired ${WORKDIR}/${zymo}_1.fastq.gz ${WORKDIR}/${zymo}_2.fastq.gz \
#   --output ${repFile}/zymo_output.kraken --report ${repFile}/zymo.kreport
#conda activate bowtieenv  
#bracken -i  ${repFile}/zymo.kreport -d $DB -o  ${repFile}/zymo.kreport.bracken -w  ${repFile}/zymo.bracken.report -t 8
#conda deactivate


repFile=report
for i in  {1..3}; do
    conda activate base
    kraken2 --db ${DB} --threads 24 --paired ${WORKDIR}/zymo_rm_rRNA_${i}.R1.fastq.gz zymo_rm_rRNA_${i}.R2.fastq.gz --output ${repFile}/zymo_${i}_output.kraken --report ${repFile}/zymo_${i}.kreport
    conda activate bowtieenv
    bracken -i  ${repFile}/zymo_${i}.kreport -d $DB -o  ${repFile}/zymo_${i}.kreport.bracken -w  ${repFile}/zymo_${i}.bracken.report -t 8
    conda deactivate

done

