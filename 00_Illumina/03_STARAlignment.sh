#!/bin/bash -l

#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J STAR_alignReads
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/STAR_alignReads_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/STAR_alignReads_%j.out
#SBATCH -t 7-00:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=16

# activate your environment
STAR="/home/fs01/moa4020/miniforge3/envs/STAR/bin"

IndexDir="/athena/cayuga_0003/scratch/moa4020/refDir/GRCm39"

fastqDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/Illumina/trim-fastq"

outDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/Illumina/STARalignment"

for trimmed_fastq in ${fastqDir}/*.fq.gz; do
    echo "Processing sample: ${trimmed_fastq}"
    $STAR/STAR --runMode alignReads \
        --runThreadN 16 \
        --genomeDir $IndexDir \
        --readFilesIn $trimmed_fastq \
        --readFilesCommand zcat \
        --outFileNamePrefix $outDir/"$(basename ${trimmed_fastq%.fq.gz})." \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes All \
        --twopassMode Basic \
        --quantMode TranscriptomeSAM GeneCounts
done

module load samtools

samtools index -M $outDir/*.Aligned.sortedByCoord.out.bam

module unload samtools