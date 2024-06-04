#!/bin/bash -l
 
#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J  barcode-rename
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/barcode-rename_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/barcode-rename_%j.out
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1

parentDir="/athena/cayuga_0003/scratch/moa4020/rr9_retina/barcode_classified_bam"

mv ${parentDir}/FC1/EXP-PBC001_barcode01.bam ${parentDir}/FC1/FC1-F9.bam
mv ${parentDir}/FC1/EXP-PBC001_barcode02.bam ${parentDir}/FC1/FC1-F11.bam
mv ${parentDir}/FC1/EXP-PBC001_barcode03.bam ${parentDir}/FC1/FC1-F15.bam
mv ${parentDir}/FC1/EXP-PBC001_barcode04.bam ${parentDir}/FC1/FC1-F16.bam
mv ${parentDir}/FC1/EXP-PBC001_barcode05.bam ${parentDir}/FC1/FC1-F17.bam
mv ${parentDir}/FC1/EXP-PBC001_barcode06.bam ${parentDir}/FC1/FC1-F18.bam
mv ${parentDir}/FC1/EXP-PBC001_barcode07.bam ${parentDir}/FC1/FC1-F19.bam
mv ${parentDir}/FC1/EXP-PBC001_barcode08.bam ${parentDir}/FC1/FC1-F20.bam
mv ${parentDir}/FC1/EXP-PBC001_barcode09.bam ${parentDir}/FC1/FC1-GC9.bam
mv ${parentDir}/FC1/EXP-PBC001_barcode10.bam ${parentDir}/FC1/FC1-GC11.bam
mv ${parentDir}/FC1/EXP-PBC001_barcode11.bam ${parentDir}/FC1/FC1-GC15.bam
mv ${parentDir}/FC1/EXP-PBC001_barcode12.bam ${parentDir}/FC1/FC1-GC16.bam

mv ${parentDir}/FC2/EXP-PBC001_barcode01.bam ${parentDir}/FC2/FC2-GC9.bam
mv ${parentDir}/FC2/EXP-PBC001_barcode02.bam ${parentDir}/FC2/FC2-GC11.bam
mv ${parentDir}/FC2/EXP-PBC001_barcode03.bam ${parentDir}/FC2/FC2-GC15.bam
mv ${parentDir}/FC2/EXP-PBC001_barcode04.bam ${parentDir}/FC2/FC2-GC16.bam
mv ${parentDir}/FC2/EXP-PBC001_barcode05.bam ${parentDir}/FC2/FC2-GC17.bam
mv ${parentDir}/FC2/EXP-PBC001_barcode06.bam ${parentDir}/FC2/FC2-GC18.bam
mv ${parentDir}/FC2/EXP-PBC001_barcode07.bam ${parentDir}/FC2/FC2-GC19.bam
mv ${parentDir}/FC2/EXP-PBC001_barcode08.bam ${parentDir}/FC2/FC2-GC20.bam
mv ${parentDir}/FC2/EXP-PBC001_barcode09.bam ${parentDir}/FC2/FC2-F9.bam
mv ${parentDir}/FC2/EXP-PBC001_barcode10.bam ${parentDir}/FC2/FC2-F11.bam
mv ${parentDir}/FC2/EXP-PBC001_barcode11.bam ${parentDir}/FC2/FC2-F15.bam
mv ${parentDir}/FC2/EXP-PBC001_barcode12.bam ${parentDir}/FC2/FC2-F16.bam

mv ${parentDir}/FC3/EXP-PBC001_barcode01.bam ${parentDir}/FC3/FC3-F17.bam
mv ${parentDir}/FC3/EXP-PBC001_barcode02.bam ${parentDir}/FC3/FC3-F18.bam
mv ${parentDir}/FC3/EXP-PBC001_barcode03.bam ${parentDir}/FC3/FC3-F19.bam
mv ${parentDir}/FC3/EXP-PBC001_barcode04.bam ${parentDir}/FC3/FC3-F20.bam
mv ${parentDir}/FC3/EXP-PBC001_barcode05.bam ${parentDir}/FC3/FC3-GC9.bam
mv ${parentDir}/FC3/EXP-PBC001_barcode06.bam ${parentDir}/FC3/FC3-GC11.bam
mv ${parentDir}/FC3/EXP-PBC001_barcode07.bam ${parentDir}/FC3/FC3-GC15.bam
mv ${parentDir}/FC3/EXP-PBC001_barcode08.bam ${parentDir}/FC3/FC3-GC16.bam
mv ${parentDir}/FC3/EXP-PBC001_barcode09.bam ${parentDir}/FC3/FC3-GC17.bam
mv ${parentDir}/FC3/EXP-PBC001_barcode10.bam ${parentDir}/FC3/FC3-GC18.bam
mv ${parentDir}/FC3/EXP-PBC001_barcode11.bam ${parentDir}/FC3/FC3-GC19.bam
mv ${parentDir}/FC3/EXP-PBC001_barcode12.bam ${parentDir}/FC3/FC3-GC20.bam

mv ${parentDir}/FC4/EXP-PBC001_barcode01.bam ${parentDir}/FC4/FC4-GC17.bam
mv ${parentDir}/FC4/EXP-PBC001_barcode02.bam ${parentDir}/FC4/FC4-GC18.bam
mv ${parentDir}/FC4/EXP-PBC001_barcode03.bam ${parentDir}/FC4/FC4-GC19.bam
mv ${parentDir}/FC4/EXP-PBC001_barcode04.bam ${parentDir}/FC4/FC4-GC20.bam
mv ${parentDir}/FC4/EXP-PBC001_barcode05.bam ${parentDir}/FC4/FC4-F9.bam
mv ${parentDir}/FC4/EXP-PBC001_barcode06.bam ${parentDir}/FC4/FC4-F11.bam
mv ${parentDir}/FC4/EXP-PBC001_barcode07.bam ${parentDir}/FC4/FC4-F15.bam
mv ${parentDir}/FC4/EXP-PBC001_barcode08.bam ${parentDir}/FC4/FC4-F16.bam
mv ${parentDir}/FC4/EXP-PBC001_barcode09.bam ${parentDir}/FC4/FC4-F17.bam
mv ${parentDir}/FC4/EXP-PBC001_barcode10.bam ${parentDir}/FC4/FC4-F18.bam
mv ${parentDir}/FC4/EXP-PBC001_barcode11.bam ${parentDir}/FC4/FC4-F19.bam
mv ${parentDir}/FC4/EXP-PBC001_barcode12.bam ${parentDir}/FC4/FC4-F20.bam
