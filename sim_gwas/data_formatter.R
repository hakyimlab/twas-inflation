#! /usr/bin/env Rscript

# load packages
library(optparse)
library(data.table)
library(tidyverse)
library(arrow)

# create options
option_list <- list(
  make_option(c("--indir"), type="character", default=NULL,
              help="Summary statistics file",
              metavar="character"),
  make_option(c("--output"), type="character", default=NULL,
              help="Summary stats file formated for PRSCS",
              metavar="character"),
  make_option(c("--region"), type="character", default=NULL,
              help="brain region you are processing (dmri,t1)",
              metavar="character"),
  make_option(c("--idp_trait"), type="character", default=NULL,
              help="the idp trait you are processing",
              metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

in.dir = opt$indir
outfile = opt$output
reg = opt$region
idp = opt$idp_trait
snp.dir = '/gpfs/data/im-lab/nas40t2/festus/metabolomics/metsim'
# collate from all chromosomes
t_ss = data.frame()
for(ch in 1:22){
  
  #pq_file = glue::glue("{in.dir}/residual.{reg}.chr{ch}/{idp}.parquet") #original
  #bim_file = glue::glue("{in.dir}/snp_bim/chr{ch}.bim")
  pq_file = glue::glue("{in.dir}/trans_qtl.{reg}.chr{ch}/{idp}.parquet")
  bim_file = glue::glue("{snp.dir}/data/train/indv_chr/metsim_train-hapmap3.chr{ch}.bim")
  # load the data
  dta = read_parquet(pq_file)
  bim = fread(bim_file) %>% 
    setnames(.,colnames(.),c("chr","rsid","cM","pos","alt","ref"))
  
  dta = bim %>% inner_join(dta, by = c("rsid" = "variant_id"))
  t_ss = t_ss %>% bind_rows(dta)
  
  print(paste0("chr",ch, "  ",nrow(dta)))
}

# format it for prscs
ff <- t_ss %>%
  #dplyr::filter(pval < (5 * 10^-4)) %>%
  dplyr::rename(SNP=rsid,A1=alt,A2=ref,BETA=b,P=pval,SE=b_se) %>%
  mutate(Z = BETA/SE) %>%
  mutate(SNP = ifelse(SNP == "", glue::glue('{chrom}:{pos}'),SNP)) %>%
  dplyr::select(SNP,CHR=chrom,POS=pos,A1,A2,Z,BETA,SE,P)

fwrite(ff, file = outfile, sep = "\t")


