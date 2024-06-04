#!/bin/bash -l

#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J STAR_genomeGenerate
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/STAR_genomeGenerate_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/STAR_genomeGenerate_%j.out
#SBATCH -t 7-00:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=16
# activate your environment
STAR="/home/fs01/moa4020/miniforge3/envs/STAR/bin"

genomeDir="/athena/cayuga_0003/scratch/moa4020/refDir/GRCm39"
GTF="/athena/cayuga_0003/scratch/moa4020/refDir/GENCODE_M33/gencode.vM33.chr_patch_hapl_scaff.annotation.gtf"

$STAR/STAR --runMode genomeGenerate \
     --runThreadN 16 \
     --genomeDir ${genomeDir} \
     --genomeFastaFiles ${genomeDir}/GRCm39.genome.fa \
     --sjdbGTFfile ${GTF} \
     --sjdbOverhang 74
