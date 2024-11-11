#!/bin/bash
#SBATCH --job-name=split_geno
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=32G
#SBATCH --error=logs/%x_%j.err
#SBATCH --output=logs/%x_%j.out

module load gcc/12.1.0
module load plink/2.0

size=100000

for chr in {1..22}
do
  plink2 --bfile /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/data/${size}/ukb_eur \
    --chr ${chr} \
    --make-bed \
    --out /scratch/fnyasimi1/ukb/ukb_eur_${size}.chr${chr}
done
