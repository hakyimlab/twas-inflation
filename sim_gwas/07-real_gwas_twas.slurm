#!/bin/bash
#SBATCH --job-name=real_gwas
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --error=logs/%x_%j.err
#SBATCH --output=logs/%x_%j.out

source ~/.bashrc
module load gcc/12.1.0
module load R/4.2.1
conda activate imlab 


trait=imputed_pgc.scz2
#trait=imputed_CKDGen_Chronic_Kidney_Disease
#trait=imputed_DIAGRAM_T2D_TRANS_ETHNIC
gwas=/gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/sim_gwas/${trait}.sumstats
name=$(basename ${gwas} .gz)
name=$(basename ${name} .sumstats)

# fusion TWAS
#for chr in {1..22}
#do
#  Rscript /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/fusion_twas/FUSION.assoc_test.R \
#  --sumstats ${gwas} \
#  --weights /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/fusion_twas/WEIGHTS/GTEx.Whole_Blood.pos \
#  --weights_dir /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/fusion_twas/WEIGHTS/ \
#  --ref_ld_chr /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/fusion_twas/LDREF/1000G.EUR. \
#  --chr ${chr} \
#  --out /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/sim_gwas/out/${name}.${chr}.dat
#done


# Spredixcan
db=/gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/results/lasso/metsim-invnorm-softimpute.db
#db=/gpfs/data/im-lab/nas40t2/Data/PredictDB/GTEx_v8/models_v1/eqtl/elastic_net_models/en_Whole_Blood.db
covs=${db%.db}.txt.gz
name=metsim-${name}

python /home/fnyasimi1/MetaXcan/software/SPrediXcan.py \
  --gwas_file ${gwas} \
  --snp_column SNP \
  --effect_allele_column A1 \
  --non_effect_allele_column A2 \
  --zscore_column Z \
  --model_db_path ${db} \
  --covariance ${covs} \
  --keep_non_rsid \
  --additional_output \
  --model_db_snp_key rsid \
  --throw \
  --output_file /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/sim_gwas/out/${name}-spred.txt

