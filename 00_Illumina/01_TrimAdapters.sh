#!/bin/bash -l

#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J TrimAdapters
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/TrimAdapters_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/TrimAdapters_%j.out
#SBATCH -t 7-00:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=4

# activate your environment
trim_galore="/home/fs01/moa4020/miniforge3/envs/trim-galore/bin"

fastqDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/Illumina/fastq"
outDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/Illumina/trim-fastq"

$trim_galore/trim_galore --gzip \
  --cores 4 \
  --phred33 \
  --illumina \
  --output_dir $outDir \
  $fastqDir/*.fastq.gz