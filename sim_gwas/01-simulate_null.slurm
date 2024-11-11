#! /bin/bash

#SBATCH --job-name=simulate_null
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=1200G
#SBATCH --partition=tier3q
#SBATCH --exclusive


###this is not used SBATCH --partition=tier3q

# SPECIFY LOGGING BEHAVIOR

#SBATCH --output=logs/%x.%j.log
#SBATCH --error=logs/%x.%j.err
source ~/.bashrc
module load gcc/12.1.0
module load R/4.2.1

size=100000
h=0.85 

# Simulate data
Rscript /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/sim_gwas/sim_null_traits.R \
  --plink_genotype /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/data/${size}/ukb_eur \
  --output_file /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/sim_gwas/sim-ukb_pheno_n${size}_h${h}.txt \
  --nsamp ${size}


