#!/bin/bash -l
 
#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J  minimap2_ONTcDNA
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/minimap2_ONTcDNA_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/minimap2_ONTcDNA_%j.out
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=16

refFASTA="/athena/cayuga_0003/scratch/moa4020/refDir/GRCm39/GRCm39.genome.fa"

refBED="/athena/cayuga_0003/scratch/moa4020/refDir/GENCODE_M33/gencode.vM33.chr_patch_hapl_scaff.annotation.bed"

inDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/FASTQ"

outDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/aligned_sam"

alignedbamDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/aligned_bam"

module load minimap2

for flowcell in ${inDir}/FC[1-4]; do
    dirName=$(basename ${flowcell})
    for file in ${flowcell}/*.fastq; do
        basename=$(basename ${file} .fastq)
        minimap2 -t 16 -ax splice --junc-bed  ${refBED} ${refFASTA} ${file} > ${outDir}/${dirName}/${basename}.sam
    done
done

module unload minimap2

module load samtools

for flowcell in ${outDir}/FC[1-4]; do
    dirName=$(basename ${flowcell})
    for samfile in ${flowcell}/*.sam; do
        filename=$(basename ${samfile} .sam)
        samtools view  -o ${alignedbamDir}/${dirName}/${filename}.bam ${samfile}
        samtools sort -o ${alignedbamDir}/${dirName}/${filename}_sorted.bam ${alignedbamDir}/${dirName}/${filename}.bam
        samtools index ${alignedbamDir}/${dirName}/${filename}_sorted.bam
    done
done

module unload samtools