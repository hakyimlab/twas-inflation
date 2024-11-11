# twas-inflation
This repository contains the code used to show inflation in TWAS and other xWAS studies. The general workflow of preprocesing data and estimating phi for different models.

## Prerequisite
* Ensure you have a proper working [PrediXcan tool](https://github.com/hakyimlab/MetaXcan/tree/master), [Fusion tool](http://gusevlab.org/projects/fusion/) and R.

* Download the GTEx models from [here](https://predictdb.org/post/2021/07/21/gtex-v8-models-on-eqtl-and-sqtl/) and Fusion models from [here](http://gusevlab.org/projects/fusion/#gtex-v8-multi-tissue-expression)


## Genotype preprocessing
We utilize the unrelayed indidivuals from the UK Biobank to do the simulations and estimate phi. First we preprocess the genotype to select individuals with low missingness in and select the common SNPs that overlap with GTEx and Hapmap3 to have a manageable size files for analysis. We sample the genotype with a varying sample size for simulations. An example below for selecting individuals with a sample size of 100k, I use a vcf file for easy maching of id between GTEx models and UK Biobank genotype

```{bash}
for i in $(seq 1 22); do

  ## make a vcf
  plink2 \
  --bgen /gpfs/data/ukb/genotypes/v3/ukb_imp_chr${i}_v3.bgen ref-first \
  --sample /gpfs/data/genotype/v3/ukb_imp_v3.sample \
  --recode vcf id-paste=iid \
  --threads 32 \
  --memory 32000 \
  --extract /gpfs/data/paper-simulations/UKB2GTEx_snps.txt \
  --keep /gpfs/data/paper-simulations/samples_100k.txt \
  --out /scratch/temp_vcf/100000/ukb_imp_chr${i}_v3

  gzip /scratch/temp_vcf/100000/ukb_imp_chr${i}_v3.vcf

done
```

## Predict gene expression
For GTEx models we use the PrediXcan tool to impute gene expression. For fusion models there is no direct way of imputing expression so we first convert the Fusion weights into a PrediXcan format to impute gene expression using the code in `fusion_to_predixcan.R`. Once converted then we can run PrediXcan to impute gene expression to be used for simulations


## Under the null simulation
First we run simulations with a null phenotype without polygenic component using the `null_null_sims.R` to obtain the Z^2 values when the phenotype has no polygenic component. Then we run a second simulation using a polygenic phenotype to obtain the Z^2 using the `polygenic_null_sims.R`. In the second simulation we run a varition of the h2 across the different sample sizes we use for simulation. The job submission scripts are in the submission folder.

## summary statistics TWAS
We run TWAS using summary statistics for both GTEx models and predixcan models using their respective tools to obtain TWAS results using null phenotypes. The workflow and code is available in the `sim_gwas` folder. The first step is we simulate a null polygenic phenotype then run a gwas using the genotype we used for simulaton to obtain the summary statistics which we use to run the TWAS.


## Show inflation and paper figures
In the figures folder we have the rmd markdown showing simulations using toy dataset and real genotype to show how TWAS studies are inflated and how to correct the them. We show results using GTEx models and Fusion models.
