#! /usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(arrow)
library(genio)
library(optparse)
library(Matrix)
library(Rfast)
#set.seed(2023)

option_list <- list(
  make_option(c("--nsim"), type="character", default=1000,
              help="The number of simulations to run",
              metavar="numeric"),
  make_option(c("--h2"), type="character", default=0.5,
              help="The heritability of null simulation",
              metavar="numeric"),
  make_option(c("--batch_size"), type="character", default=500,
              help="Size of each batch for simulation",
              metavar="numeric"),
  make_option(c("--predicted_traits"), type="character", default=NULL,
              help="The predicted molecular trait",
              metavar="character"),
  make_option(c("--plink_genotype"), type="character", default=NULL,
              help="The prefix for the genotype used in prediction in plink format",
              metavar="character"),
  make_option(c("--output_file"), type="character", default=NULL,
              help="The output file with average chi2",
              metavar="character"),
  make_option(c("--regress_pcs"), type="character", default=NULL,
              help="Option to regress out the first 5 principal components (TRUE or FALSE)",
              metavar="logical")
  )

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)

nsim <- as.numeric(args$nsim)
predicted_traits <- args$predicted_traits
plink_file_path <- args$plink_genotype
output_file <- args$output_file
regress_pcs <- args$regress_pcs
h <- as.numeric(args$h2)
batch_size <- as.numeric(args$batch_size)

#####
# cor to chi2
## calculate p-value from correlation
cor2zscore = function(cc,nn) 
{
  zz = atanh(cc) * sqrt(nn-3)
}

cor2pval = function(cc,nn) 
{
  zz=cor2zscore(cc,nn)
  pnorm(-abs(zz))*2
}

cor2chi2 = function(cc,n)
{
  cor2zscore(cc,n)^2
}

# genotype loading function
load_geno <- function(plink_file_path) {
  # load plink genotype
  plk_dta <- read_plink(plink_file_path)
  X <- t(plk_dta$X)
  # count the NAs in the genotype
  cc = sum(is.na(X))
  print(glue::glue("Total NAs in genotype: {cc}, imputed by the colmean"))
  # mutate the NAs with col mean
  #X[is.na(X)] <- 0
  X = as.data.table(X)
  for (j in seq_along(X)) {
     set(X, which(is.na(X[[j]])), j, mean(X[[j]], na.rm = TRUE))
  }
  X <- as.matrix(X)
  return(X)
}

# simulating null y function
gen_y_h2 <- function(nsim,plink_pre,h2Y) {
  # split job into chrs
  for (chr in seq(1:22)) {
    # load genotype
    plink_file_path <- glue::glue("{plink_pre}{chr}_v3")
    print(glue::glue("Working on {plink_file_path}"))
    X = load_geno(plink_file_path)
    msnp = ncol(X)
    nsamp = nrow(X)
    betamat = matrix(rnorm(msnp*nsim),msnp, nsim)
    temp_gYmat = X %*% betamat
    
    if (chr == 1) { gYmat = temp_gYmat } else { gYmat = gYmat + temp_gYmat }
  }
  
  # random noise
  epsimat = matrix(rnorm(nsamp*nsim),nsamp, nsim)
  epsimat = scale(epsimat) * sqrt(1 - h2Y)
  # scale with h2
  gYmat = scale(gYmat) * sqrt(h2Y)
  Ymat = gYmat + epsimat
  return(Ymat)
}


# function to calculate chi2
gen_null_chi2_fast <- function(plink_file_path, expr_mat, h2_Y, nsim = 1000) {
  Ymat_null = gen_y_h2(nsim,plink_file_path,h2_Y)
  if (!is.null(regress_pcs)) {
     print("Regressing PCs")
     Ymat_null = regress_pca(Ymat_null,g_pcs)
  }
  nsamp = nrow(expr_mat)
  calculate_chi(expr_mat, Ymat_null,nsamp)
}
  
calculate_chi <- function(expr_mat, Y_mat,nsamp) {
  # scale both X and Y
  X_scale = expr_mat %>% scale()
  Y_mat_h2 = scale(Y_mat)
  
  # transpose and multiply
  cor_mat = (t(X_scale) %*% Y_mat_h2)/nsamp

  # calculate the chi2
  chimat = apply(cor_mat, 2, cor2chi2, n = nsamp)
  
  # transpose sims on rows genes on  columns
  chimat <- chimat %>% data.frame() %>% t() 
  
  return(chimat)
}

# regressing out pcs
regress_pca <- function(Ymat,geno_pcs, npcs = 5) {
 
  # Load pcs
  pca_res <- geno_pcs %>% 
    unite(ID, FID, IID, sep = "_") %>% 
    column_to_rownames(var = "ID")
  
  df_pcs <- pca_res[,c(1:npcs)]
  
  # quick data formating
  pcs = colnames(df_pcs)
  pcs <- do.call(paste, c(as.list(pcs),sep = " + "))
  colnames(Ymat) = paste0("Y", c(1:ncol(Ymat)))
  phenos = colnames(Ymat)

  # Generate empty data frames
  res_df <- data.frame(matrix(,nrow=nrow(Ymat),ncol=ncol(Ymat)))
  rownames(res_df) <- rownames(Ymat)
  df_w_pc = merge(Ymat, df_pcs, by = 'row.names') %>% 
    column_to_rownames(var = "Row.names")

  for (i in 1:ncol(Ymat)){
    comp = phenos[i]
    #print(comp)
    test <- lm(paste0(comp, " ~ ", pcs),data = df_w_pc)
    preds <- data.frame(predict(test))
    resids <- data.frame(resid(test))
    
    colnames(preds) <- comp
    colnames(resids) <- comp
    res_df[i] <- resids
  }
  return(res_df)
}


#######
if(grepl("\\.parquet$", predicted_traits)) {
  pred_traits <- read_parquet(predicted_traits)
} else {
  pred_traits <- fread(predicted_traits) %>% 
    filter(!is.na(IID))
}


if ("indiv" %in% names(pred_traits)) {
  pred_traits <- pred_traits %>%
    dplyr::filter(indiv %in% select_samp$IID) %>%
    column_to_rownames(var = "indiv") %>%
    as.matrix()
} else {
  pred_traits <- pred_traits %>% select(-FID) %>%
    dplyr::filter(IID %in% select_samp$IID) %>% 
    column_to_rownames(var = "IID") %>%
    as.matrix()
}

# loading pcs
if (!is.null(regress_pcs)) {
  g_pcs <- fread(regress_pcs) %>% 
    setnames(., names(.), c("FID","IID", paste0("PC", c(1:20))))
}

# run assoc
n_batches=nsim/batch_size
file_name <- basename(output_file)

for (batch in 1:n_batches) {
  print(glue::glue("Running batch {batch}/{n_batches}"))
  # run assoc
  traits_nulls <- gen_null_chi2_fast(plink_file_path,pred_traits,h,batch_size)
  
  # save the nulls
  fwrite(traits_nulls,
         file = glue::glue("/scratch/fnyasimi1/tmp_data/batch{batch}_{file_name}"),
         sep = "\t")
}

# clean up
rm(X,pred_traits)

# gather
for (batch in 1:n_batches) {
  file_path <- glue::glue("/scratch/fnyasimi1/tmp_data/batch{batch}_{file_name}")
  tmp_df <- fread(file_path)

  if (batch == 1) {
    all_data <- tmp_df 
  } else {
    all_data <- all_data %>% 
      bind_rows(tmp_df)
  }
  file.remove(file_path)
}

# save the nulls
fwrite(all_data,
       file = output_file,
       sep = "\t")


