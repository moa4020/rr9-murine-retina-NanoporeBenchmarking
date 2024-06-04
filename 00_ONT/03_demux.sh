#!/bin/bash -l

#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J dorado_basecall-demux
#SBATCH -p scu-gpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/dorado_basecall-demux_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/dorado_basecall-demux_%j.out
#SBATCH -t 72:00:00
#SBATCH --mem-per-gpu=40G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --gpus-per-task=2

doradoBin="/athena/cayuga_0003/scratch/moa4020/dorado-0.5.0-linux-x64/bin"
inputDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/pod5"
unclassifiedDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/unclassified_bam"
demuxDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/barcode_classified_bam"

for folder in ${inputDir}/FC[1-4]; do
    # Save folder name
    folderName=$(basename ${folder})

    # Basecall
    ${doradoBin}/dorado basecaller ${doradoBin}/dna_r9.4.1_e8_sup@v3.6 ${folder} --kit-name EXP-PBC001 > ${unclassifiedDir}/${folderName}/${folderName}.bam

    # Sort the aligned BAM - supposed to increase the number of demux-ed reads
    # Uncomment and modify the following lines if samtools is required
    # module load samtools
    # samtools sort -o ${unclassifiedDir}/${folderName}/${folderName}_sorted.bam ${unclassifiedDir}/${folderName}/${folderName}.bam
    # module unload samtools

    # Demux
    ${doradoBin}/dorado demux --kit-name EXP-PBC001 --output-dir ${demuxDir}/${folderName}/ ${unclassifiedDir}/${folderName}/${folderName}.bam
done