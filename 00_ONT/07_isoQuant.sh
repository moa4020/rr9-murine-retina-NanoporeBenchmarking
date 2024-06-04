#!/bin/bash -l

#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J isoQuant
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/isoQuant_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/isoQuant_%j.out
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=16

module load isoquant

alignedBamDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/aligned_bam"

refFASTA="/athena/cayuga_0003/scratch/moa4020/refDir/GRCm39/GRCm39.genome.fa"

refGTF="/athena/cayuga_0003/scratch/moa4020/refDir/GENCODE_M33/gencode.vM33.chr_patch_hapl_scaff.annotation.gtf"

outDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/featureCounts_isoQuant"

for flowCell in ${alignedBamDir}/FC[1-4]; do
    flowcell=$(basename ${flowCell})
    isoquant.py --force \
        --threads 16 \
        --data_type nanopore \
        --bam_list ${flowCell}/rr9-retina_isoQuant_${flowcell}.txt \
        --reference ${refFASTA} \
        --genedb ${refGTF} \
        --complete_genedb \
        --transcript_quantification all \
        --gene_quantification all \
        --count_exons \
        --output ${outDir}/${basename}
done

module unload isoquant