#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --error=logs/run_2nd_%x_%j.err
#SBATCH --output=logs/run_2nd_%x_%j.out


source ~/.bash_profile
source ~/.bashrc

#SNMK=/gpfs/data/im-lab/nas40t2/yanyul/softwares/miniconda2/envs/mixqtl/bin/snakemake
SNMK=/gpfs/data/im-lab/nas40t2/festus/metabolomics/guardian/MetaboXcan/brainxcan/snakemake

#conda activate tensorqtl
conda activate /gpfs/data/im-lab/nas40t2/yanyul/softwares/miniconda2/envs/tensorqtl

cd ${SLURM_SUBMIT_DIR}
pipe=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/gw_qtl/run.snmk


$SNMK -s $pipe --configfile ${SLURM_SUBMIT_DIR}/configs/config.$NAME.yaml -p \
  --config chr_num=$CHR \
  > ${SLURM_SUBMIT_DIR}/logs/run_2nd_${NAME}_${CHR}.log 2>&1

