module load gcc/12.1.0 

source ~/.bashrc
conda activate imlab

h2=0.85
N=100000
sim=9

model=metsim
db=/gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/results/lasso/metsim-invnorm-softimpute.db
#db=/gpfs/data/im-lab/nas40t2/Data/PredictDB/GTEx_v8/models_v1/eqtl/elastic_net_models/en_Whole_Blood.db
covs=${db%.db}.txt.gz
gwas=/gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/sim_gwas/n${N}_sim${sim}-${h2}.sumstats

#dbname=$(basename ${db} .db)
#outname=$(basename ${gwas} .txt.gz)_${dbname}.csv

echo Running tissue ${gwas}

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
  --output_file /gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim/ukb/paper-simulations/sim_gwas/out/n${N}_sim${sim}-${h2}-${model}-spred.txt
