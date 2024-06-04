#!/bin/bash -l

#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J  extract-fast5_convert-to-pod5
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/extract-fast5_convert-to-pod5_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/extract-fast5_convert-to-pod5_%j.out
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=12

rawData="/athena/cayuga_0003/scratch/moa4020/rr9_retina/rawdata"
outDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/pod5"

for folder in ${rawData}/Retina_FlowCell[1-4]/*/*/fast5_pass; do
    cd ${folder}
    pod5 convert fast5  *.fast5 --output . --one-to-one .
    
    mv ${rawData}/Retina_FlowCell1/*/*/fast5_pass/*.pod5 ${outDir}/FC1
    mv ${rawData}/Retina_FlowCell2/*/*/fast5_pass/*.pod5 ${outDir}/FC2
    mv ${rawData}/Retina_FlowCell3/*/*/fast5_pass/*.pod5 ${outDir}/FC3
    mv ${rawData}/Retina_FlowCell4/*/*/fast5_pass/*.pod5 ${outDir}/FC4
done

sbatch /athena/cayuga_0003/scratch/moa4020/scripts/dorado_basecall-demux.sh
