#!/bin/bash -l

#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J rsem_counts
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/rsem_counts_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/rsem_counts_%j.out
#SBATCH -t 7-00:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=16

GTF="/athena/cayuga_0003/scratch/moa4020/refDir/GENCODE_M33/gencode.vM33.chr_patch_hapl_scaff.annotation.gtf"
refDir="/athena/cayuga_0003/scratch/moa4020/refDir/GRCm39"
rsem="/home/fs01/moa4020/miniforge3/envs/rsem/bin"
alignmentDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/Illumina/STARalignment"
outDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/Illumina/RSEM"

$rsem/rsem-prepare-reference --gtf $GTF \
    $refDir/GRCm39.genome.fa  \
    $refDir/RSEM-ref-GRCm39

for bamFile in $alignmentDir/*_trimmed.Aligned.toTranscriptome.out.bam; do
    prefix=$(basename $bamFile %_trimmed.Aligned.toTranscriptome.out.bam)
    $rsem/rsem-calculate-expression --num-threads 16 \
        --alignments $bamFile \
        --seed 1245 \
        --seed-length 20 \
        --estimate-rspd \
        --no-bam-output \
        $refDir/RSEM-ref-GRCm39 \
        $outDir/$prefix
done