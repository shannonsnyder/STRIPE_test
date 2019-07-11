#!/bin/bash

#SBATCH -n 2                        # number of cores
#SBATCH -t 0-12:00                  # wall time (D-HH:MM)
##SBATCH -A ssnyde11                # Account hours will be pulled from (commented out with double # in front)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification

module load java/latest
module load gcc/7x
module load r/3.5.2
module load bwakit/0.7.12
module load samtools/1.8

echo "Launching job"

cd /home/ssnyde11/scratch/STRIPE_demult_test/scripts #update to reflect our directory structure
./run_tagdust_BWA_alignment.sh > err

echo "Job complete"
