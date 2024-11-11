#!/bin/bash

#SBATCH --job-name=fusion_null
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=6:00:00
#SBATCH --partition=tier3q
#SBATCH --mem=1200G

# SPECIFY LOGGING BEHAVIOR

#SBATCH --output=logs/%x.%j.log
#SBATCH --error=logs/%x.%j.err

source ~/.bashrc
module load gcc/12.1.0
module load R/4.2.1
conda activate imlab

db=/gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/fusion_twas/WEIGHTS/GTEx.Whole_Blood.db
name=fusion-$(basename ${db} .db)
size=10000
h2Y=0.5
echo Working on ${name} with size ${size}

# predict
python /home/fnyasimi1/MetaXcan/software/Predict.py \
  --model_db_path ${db} \
  --model_db_snp_key rsid \
  --vcf_genotypes /scratch/fnyasimi1/temp_vcf/${size}/ukb_imp_chr*_v3.vcf.gz \
  --vcf_mode genotyped \
  --prediction_output /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/predictions/${name}-${size}-hapmap.txt \
  --prediction_summary_output /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/predictions/${name}-${size}-hapmap_summary.txt \
  --skip_palindromic \
  --verbosity 9 \
  --throw

# get chi for null null
Rscript /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/null_null_sims.R \
   --predicted_traits /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/predictions/${name}-${size}-hapmap.txt \
   --nsim 1000 --h2 ${h2Y} \
   --output_file /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/nullchi2/full-${name}-${size}-${h2Y}-null-null.txt


# get chi for polygenic null
Rscript /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/polygenic_null_sims.R \
   --predicted_traits /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/predictions/${name}-${size}-hapmap.txt \
   --plink_genotype /scratch/fnyasimi1/temp_vcf/${size}/ukb_imp_chr \
   --nsim 1000 --h2 ${h2Y} \
   --output_file /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/nullchi2/full-${name}-${size}-${h2Y}-polygenic-null.txt

