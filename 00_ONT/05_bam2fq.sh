#!/bin/bash -l
 
#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J  bamtofq
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/bamtofq_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/bamtofq_%j.out
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=10

# Specify the parent directory containing the subdirectories with BAM files
parentDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/barcode_classified_bam"

outDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/FASTQ"

module load samtools

# Loop through each subdirectory
for dir in ${parentDir}/FC[1-4]; do
    dirName=$(basename ${dir})
    # Navigate into the subdirectory
    cd ${dir}
    # Convert BAM files to FASTQ and save in the same directory
    for file in *.bam; do
        samtools fastq ${file} > ${outDir}/${dirName}/${file%.bam}.fastq
    done
    # Navigate back to the parent directory
    cd ${parentDir}
done

module unload samtools
