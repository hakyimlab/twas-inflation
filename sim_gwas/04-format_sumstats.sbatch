#!/bin/bash
#SBATCH --job-name=gather_gwas
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --error=logs/%x_%j.err
#SBATCH --output=logs/%x_%j.out

source ~/.bashrc
module load gcc/12.1.0
module load R/4.2.1 

sim=9
h2=0.85
size=100000

Rscript /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/sim_gwas/data_formatter.R \
  --indir /scratch/fnyasimi1/ukb_idp/idp_gwas_4th \
  --region ${size} \
  --idp_trait n${size}_sim${sim} \
  --output /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/sim_gwas/n${size}_sim${sim}-${h2}.sumstats
