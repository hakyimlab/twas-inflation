#!/bin/bash

#SBATCH --job-name=tissue_pred
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=120:00:00
#SBATCH --partition=tier2q
#SBATCH --mem=480G
#SBATCH --exclusive

# SPECIFY LOGGING BEHAVIOR

#SBATCH --output=logs/%x.%j.log
#SBATCH --error=logs/%x.%j.err

source ~/.bashrc
module load gcc/12.1.0
module load R/4.2.1
conda activate imlab

db=/gpfs/data/im-lab/nas40t2/Data/PredictDB/GTEx_v8/models_v1/eqtl/elastic_net_models/en_Whole_Blood.db
working_dir=/gpfs/data/im-lab/nas40t2/festus # update accordingly

name=$(basename ${db} .db)
size=100000
h2Y=0.9
echo Working on ${name} with size ${size}


# predict using vcf to lift over ~98%
python /home/fnyasimi1/MetaXcan/software/Predict.py \
  --model_db_path ${db} \
  --model_db_snp_key varID \
  --vcf_genotypes /scratch/fnyasimi1/temp_vcf/${size}/ukb_imp_chr*_v3.vcf.gz \
  --vcf_mode genotyped \
  --liftover /gpfs/data/im-lab/nas40t2/festus/metabolomics/liftover_chains/hg19ToHg38.over.chain.gz \
  --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
  --prediction_output ${working_dir}/predictions/${name}-${size}-hapmap.txt \
  --prediction_summary_output ${working_dir}/predictions/${name}-${size}-hapmap_summary.txt \
  --verbosity 9 \
  --throw

# get chi2 for null models
Rscript /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/null_null_sims.R \
   --predicted_traits  ${working_dir}/predictions/${name}-${size}-hapmap.txt \
   --nsim 1000 --h2 ${h2Y} \
   --output_file ${working_dir}/nullchi2/${name}-${size}-${h2Y}-null-null.txt

# get chi for polygenic
Rscript /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/polygenic_null_sims.R \
   --predicted_traits ${working_dir}/predictions/${name}-${size}-hapmap.txt \
   --plink_genotype /scratch/fnyasimi1/temp_vcf/${size}/ukb_imp_chr \
   --nsim 1000 --h2 ${h2Y} \
   --output_file ${working_dir}/nullchi2/${name}-${size}-${h2Y}-polygenic-null.txt

# estimate phi and insert to model
Rscript /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/est_phi.R \
   --chi_table ${working_dir}/nullchi2/${name}-${size}-${h2Y}-polygenic-null.txt \
   --sample_size ${size} --h2 ${h2Y} \
   --input_db ${db} \
   --output_db ${working_dir}/dbs/${name}.db

# submission
#for db in $(ls /gpfs/data/im-lab/nas40t2/Data/PredictDB/GTEx_v8/models_v1/eqtl/elastic_net_models/en_*.db); do dn=$(basename ${db} .db); sbatch --job-name=${dn} --export=db=${db} full_gene_pred.slurm; done


