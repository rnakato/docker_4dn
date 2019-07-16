#!/bin/bash

fastq1=$1
fastq2=$2
index=$3
prefix=$4
nThreads=$5

if [[ $outdir != '.' ]]
then
  cd $outdir
fi

# unzip fastq files
if [[ $fastq1 =~ \.gz$ ]]
then
  cp $fastq1 fastq1.gz
  gunzip fastq1.gz
else
  cp $fastq1 fastq1
fi
  fastq1=fastq1

if [[ $fastq2 =~ \.gz$ ]]
then
  cp $fastq2 fastq2.gz
  gunzip fastq2.gz
else
  cp $fastq2 fastq2
fi
  fastq2=fastq2


# run bwa
bwa mem -t $nThreads -SP5M $index $fastq1 $fastq2 | samtools view -Shb - > $prefix.bam
