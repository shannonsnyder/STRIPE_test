#!/bin/bash

BWAdir=/home/ssnyde11/scratch/Dpul_Dobt_fastqs/demultiplexed_matched

TagdustDir=/home/ssnyde11/scratch/Dpul_Dobt_fastqs/demultiplexed_matched

module load fastqc/0.11.7

#cd $BWAdir

#echo "Starting fastqc job"

#for fq in unknown_*.fastq

# do

#fastqc $fq

#done

cd $TagdustDir

for fq in *_trno_tagdusted_READ?.fq

do

fastqc $fq

done

echo "Job complete"
