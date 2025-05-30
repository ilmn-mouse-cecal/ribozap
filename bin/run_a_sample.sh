#!/bin/bash
set -euo pipefail

# Variable Set Up
i=${1}
sortMeRnaPath=resources/
thisPath=Analysis/Sample${i}
sampleName=Sample${i}

### Step1: Prepare ###
# Copy Fastq files, change name
mkdir -p ${thisPath}
cp data/Sample${i}/*_R1.fastq ${thisPath}/${sampleName}_R1.fastq
cp data/Sample${i}/*_R2.fastq ${thisPath}/${sampleName}_R2.fastq

bin/merge-paired-reads.sh ${thisPath}/${sampleName}_R1.fastq ${thisPath}/${sampleName}_R2.fastq ${thisPath}/${sampleName}_Merged.fastq

### Step2: SortMeRna ###
sortmerna \
  --workdir '/app/' \
  --ref ${sortMeRnaPath}/rRNA_databases/silva-arc-23s-id98.fasta \
  --ref ${sortMeRnaPath}/rRNA_databases/silva-bac-23s-id98.fasta \
  --ref ${sortMeRnaPath}/rRNA_databases/silva-bac-16s-id90.fasta \
  --ref ${sortMeRnaPath}/rRNA_databases/rfam-5.8s-database-id98.fasta \
  --ref ${sortMeRnaPath}/rRNA_databases/silva-euk-18s-id95.fasta \
  --ref ${sortMeRnaPath}/rRNA_databases/rfam-5s-database-id98.fasta \
  --ref ${sortMeRnaPath}/rRNA_databases/silva-arc-16s-id95.fasta \
  --ref ${sortMeRnaPath}/rRNA_databases/silva-euk-28s-id98.fasta \
  --reads ${thisPath}/${sampleName}_Merged.fastq \
  --aligned ${thisPath}/${sampleName}_SortMeRna \
  --sam \
  --SQ \
  --num_alignments 0


samtools view -Sb ${thisPath}/${sampleName}_SortMeRna.sam > ${thisPath}/${sampleName}_SortMeRna.bam

samtools sort \
  ${thisPath}/${sampleName}_SortMeRna.bam \
  -o ${thisPath}/${sampleName}_SortMeRna.sorted.bam

samtools index ${thisPath}/${sampleName}_SortMeRna.sorted.bam

### Step3: Get Coverage ###
# get coverage
genomeCoverageBed -bga -ibam ${thisPath}/${sampleName}_SortMeRna.sorted.bam -g Genome/allFasta.fasta > ${thisPath}/${sampleName}_genomeCoverage.bed

# Pick Up High Coverage Regions
bin/identify_blocks.py -s ${sampleName} -d ${thisPath} --high 500 --mid 100
bin/merge_closeBy_blocks.py -s ${sampleName} -d ${thisPath}
tail -n +2 ${thisPath}/${sampleName}_highCoverageBlocksGapMerged.tsv | sort -k 4 -nr > ${thisPath}/${sampleName}_highCoverageBlocksGapMerged_CovSorted.bed

# Adding Proportion of the coverage
readCount=`wc -l ${thisPath}/${sampleName}_Merged.fastq | awk '{print $1/4;}'`
bin/add_readPercent.py -s ${sampleName} -d ${thisPath} -n ${readCount}
