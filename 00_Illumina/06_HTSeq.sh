#!/bin/bash -l

#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J DEXSeq_countReads
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/DEXSeq_countReads_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/DEXSeq_countReads_%j.out
#SBATCH -t 7-00:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=8

cd /athena/cayuga_0003/scratch/moa4020/rr9_retina/Illumina/STARalignment

samtools  view -h *.Aligned.sortedByCoord.out.bam > *.Aligned.sortedByCoord.out.sam

Prepare annotations
python3 ./dexseq_prepare_annotation.py /athena/cayuga_0003/scratch/moa4020/refDir/GENCODE_M33/gencode.vM33.chr_patch_hapl_scaff.annotation.gtf gencode.vM33.chr_patch_hapl_scaff.annotation.DEXSeq.gff

scriptDir="/home/fs01/moa4020/miniforge3/envs/HTSeq/bin"

cd ${scriptDir}

mamba activate HTSeq

workingDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/Illumina/STARalignment"

for samfile in "${workingDir}"/*trimmed.Aligned.sortedByCoord.out.sam; do
    if [ -e "$samfile" ]; then
        baseName=$(basename "$samfile" .sam)
        outputTxt="${baseName}.txt"
        python3 "${scriptDir}/dexseq_count.py" "${workingDir}/gencode.vM33.chr_patch_hapl_scaff.annotation.DEXSeq.gff" "${samfile}" ${workingDir}/"${outputTxt}"
    else
        echo "No matching files found for pattern: $samfile"
    fi
done

mamba deactivate
