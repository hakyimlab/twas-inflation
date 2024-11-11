#!/bin/bash
#SBATCH --job-name=fusion_assoc
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --error=logs/%x_%j.err
#SBATCH --output=logs/%x_%j.out

source ~/.bashrc
module load gcc/12.1.0
module load R/4.2.1 

h2=0.85
N=100000
sim=9

for chr in {1..22}
do
  Rscript /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/fusion_twas/FUSION.assoc_test.R \
  --sumstats /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/sim_gwas/n${N}_sim${sim}-${h2}.sumstats \
  --weights /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/fusion_twas/WEIGHTS/GTEx.Whole_Blood.pos \
  --weights_dir /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/fusion_twas/WEIGHTS/ \
  --ref_ld_chr /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/fusion_twas/LDREF/1000G.EUR. \
  --chr ${chr} \
  --out /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/sim_gwas/out/n${N}_sim${sim}-${h2}.${chr}.dat
done


