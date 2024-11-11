#! /usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(arrow)
library(BEDMatrix)
library(optparse)
library(Matrix)
library(Rfast)
set.seed(2023)

option_list <- list(
  make_option(c("--nsim"), type="character", default=1000,
              help="The number of simulations to run",
              metavar="numeric"),
  make_option(c("--h2"), type="character", default=0.5,
              help="The heritability of null simulation",
              metavar="numeric"),
  make_option(c("--predicted_traits"), type="character", default=NULL,
              help="The predicted molecular trait",
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
output_file <- args$output_file
regress_pcs <- args$regress_pcs
h <- as.numeric(args$h2)

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


# UNDER INVESTIGATION
#cor2chi2 <- function(r,n) {
#  # r: correlation coefficient
#  # n: sample size
#  chi2 <- (n - 2) * r^2 / (1 - r^2)
#  return(chi2)
#}

# simulate the y
gen_y_h2 <- function(nsim,pred_expr,h2Y) {

  nsamp = nrow(pred_expr)

  betamat = matrix(rt(nsamp*nsim, df = 2),nsamp, nsim)
  epsimat = matrix(rt(nsamp*nsim, df = 2),nsamp, nsim)
  epsimat = scale(epsimat) * sqrt(1 - h2Y)
  
  Ymat = betamat + epsimat

  return(Ymat)
}

# function to calculate chi2
gen_null_chi2_fast <- function(expr_mat, h2_Y, nsim = 1000) {
  Ymat_null = gen_y_h2(nsim,expr_mat,h2_Y)
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
  
  # make the matrices sparse
  #X_scale = Matrix(X_scale, sparse = TRUE)
  #Y_mat_h2 = Matrix(Y_mat_h2, sparse = TRUE)

  # transpose and multiply
  cor_mat = (t(X_scale) %*% Y_mat_h2)/nsamp
  #cor_mat = mat.mult(t(X_scale), Y_mat_h2)/nsamp

  # make the matrix dense again
  #cor_mat = as.matrix(cor_mat)
  
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
    column_to_rownames(var = "indiv") %>%
    as.matrix()
} else {
  pred_traits <- pred_traits %>% select(-FID) %>% 
    column_to_rownames(var = "IID") %>%
    as.matrix()
}

# loading pcs
if (!is.null(regress_pcs)) {
  g_pcs <- fread(regress_pcs) %>% 
    setnames(., names(.), c("FID","IID", paste0("PC", c(1:20))))
}

# sim the null pheno
# run assoc
traits_nulls <- gen_null_chi2_fast(pred_traits,h,nsim)
  
# save the nulls
fwrite(traits_nulls, 
       file = output_file,
       sep = "\t")


