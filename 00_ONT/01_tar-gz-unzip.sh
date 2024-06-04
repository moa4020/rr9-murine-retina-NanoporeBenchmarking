#!/bin/bash -l

#SBATCH --mail-user=moa4020@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -J  tar-gz-unzip
#SBATCH -p scu-cpu
#SBATCH -e /athena/cayuga_0003/scratch/moa4020/scripts/err/tar-gz-unzip_%j.err
#SBATCH -o /athena/cayuga_0003/scratch/moa4020/scripts/out/tar-gz-unzip_%j.out
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=16

# Set the directory containing the tar.gz files
tar_files_directory="/athena/cayuga_0003/scratch/moa4020/rr9_retina/rawdata/"

# Change to the directory
cd "$tar_files_directory" || exit

# Loop through each tar.gz file
for tar_file in *.tar.gz; do
    # Extract the contents of the tar.gz file
    tar -xzvf "$tar_file"
done

# Optional: Remove the tar.gz files after extraction
# rm *.tar.gz

mv Retina_RR9_FC4 Retina_FlowCell4
mv RR9_Retina_FC3 Retina_FlowCell3

rm */*.pdf
