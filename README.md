# twas-inflation
This repository contains the code used to show inflation in TWAS and other xWAS studies. The general workflow of preprocesing data and estimating phi for different models.

## Prerequisite
* Ensure you have a proper working [MetaXcan software](https://github.com/hakyimlab/MetaXcan/tree/master), [Fusion tool](http://gusevlab.org/projects/fusion/) and R.

* Download the GTEx models from [here](https://predictdb.org/post/2021/07/21/gtex-v8-models-on-eqtl-and-sqtl/) and Fusion models from [here](http://gusevlab.org/projects/fusion/#gtex-v8-multi-tissue-expression)


## Genotype preprocessing
We utilize the unrelated individuals from the UK Biobank to do the simulations and estimate phi. First, we preprocess the genotype to select individuals with low missingness and select the common SNPs that overlap with GTEx and Hapmap3 to have a manageable size of files for analysis. We sample the genotype with a varying sample size for simulations. An example below for selecting individuals with a sample size of 100k, I use a vcf file for easy matching of id between GTEx models and UK Biobank genotype

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
For GTEx models we use the PrediXcan tool to impute gene expression. For Fusion models there is no direct way of imputing expression so we first convert the Fusion weights into a PrediXcan format to impute gene expression using the code in `fusion_to_predixcan.R`. Once converted then we can run PrediXcan to impute gene expression to be used for simulations


## Under the null simulation
First, we run simulations with a null phenotype without polygenic component using the `null_null_sims.R` to obtain the Z^2 values when the phenotype has no polygenic component. Then we run a second simulation using a polygenic phenotype to obtain the Z^2 using the `polygenic_null_sims.R`. In the second simulation, we run various values of the h2 across the different sample sizes we use for simulation. The job submission scripts are in the submission folder.

## summary statistics TWAS
We run TWAS using summary statistics for both GTEx models and Fusion models using their respective tools to obtain TWAS results using null phenotypes. The workflow and code are available in the `sim_gwas` folder. The first step is to simulate a null polygenic phenotype and then run a GWAS using the genotype we used for simulation to obtain the summary statistics for running the TWAS.


## Show inflation and paper figures
In the figures folder, we have the rmd markdown showing simulations using toy dataset and real genotype to show how TWAS studies are inflated and how to correct them. We show results using GTEx models and Fusion models.

## Application to Real TWAS
We updated our [MetaXcan software v0.8.0](https://github.com/hakyimlab/MetaXcan/releases/tag/v0.8.0) to perform variance control for TWAS studies when using summary statistics. 
First you will need download the latest release of MetaXcan, then download the [latest tissue models](https://uchicago.box.com/s/w0nzszuvuwcsznvo8x4c15o4hujqrwm7) with the phi used for varinace control.

Next determine the sample size (N) and the heritability (h2) of your GWAS summary statistics, these two parameters together with phi are required for variance control. 
To estimate the h2 of your GWAS summary statistics you can use [LDSC](https://github.com/bulik/ldsc) as described in their method. 
Once you have all the required parameters you can run SPrediXcan with two new additional parameters as below;

```{bash}
python /MetaXcan/software/SPrediXcan.py \
  --gwas_file imputed_pgc.scz2.txt.gz \
  --gwas_N 3546445 \
  --gwas_h2 0.34 \
  --snp_column SNP \
  --effect_allele_column A1 \
  --non_effect_allele_column A2 \
  --zscore_column Z \
  --model_db_path en_Whole_Blood.db \
  --covariance en_Whole_Blood.txt.gz \
  --model_db_snp_key rsid \
  --keep_non_rsid \
  --additional_output \
  --throw \
  --output_file imputed_pgc.scz2_SprediXcan_results.csv
```

The output will contain the calibrated values for pvalue and zscore are in the `pvalue` and `zscore` columns respectively. 
This is to allow for backward compatibility of the software and also not to break downstream pipelines which depend on the SPrediXcan results.
The uncalibrated pvalues and zscores are at the right end of the table in `uncalibrated_pvalue` and `uncalibrated_zscore` columns respectively.

**Note:** If you don't provide both the `--gwas_N` and `--gwas_h2` the tool will give you uncalibrated_pvalue and uncalibrated_zscore in the pvalue and zscore column. 
For calibration to work both the parameters should be provided and used with latest models with phi parameter.  

