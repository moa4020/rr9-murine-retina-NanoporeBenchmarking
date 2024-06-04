#!/bin/bash -l

#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J Flagstat
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/samtoolsFlagstat_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/samtoolsFlagstat_%j.out
#SBATCH -t 7-00:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1

module load samtools

inDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/Illumina/STARalignment"
outDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/Illumina/STARalignment/flagstat"

# Iterate over the bam files in the directory
for file in ${inDir}/*.Aligned.sortedByCoord.out.bam; do
  filename=$(basename ${file})
  samtools flagstat ${file} > ${outDir}/${filename}.flagstat.txt
  samtools stats ${file} > ${outDir}/${filename}.stats.txt
  echo "$filename done"
done

module unload samtools