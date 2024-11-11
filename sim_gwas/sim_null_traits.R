#! /usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(arrow)
library(BEDMatrix)
library(optparse)

option_list <- list(
  make_option(c("--nsamp"), type="character", default=1000,
              help="The sample size you are processing",
              metavar="numeric"),
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

nsamp <- as.numeric(args$nsamp)
plink_file_path <- args$plink_genotype
output_file <- args$output_file
regress_pcs <- args$regress_pcs

# function
# simulate the y
gen_y_h2 <- function(nsim,X,h2Y) {

  nsamp = nrow(X)
  msnp = ncol(X)

  betamat = matrix(rnorm(msnp*nsim),msnp, nsim)
  epsimat = matrix(rnorm(nsamp*nsim),nsamp, nsim)
  epsimat = scale(epsimat) * sqrt(1 - h2Y)
  gYmat = X %*% betamat
  gYmat = scale(gYmat) * sqrt(h2Y)
  Ymat = gYmat + epsimat

  return(Ymat)
}

# sim the null pheno
X <- BEDMatrix(plink_file_path)
X <- X %>% as.matrix()
colnames(X) <- gsub("\\_.*", "",colnames(X))
# fill missing genotype
#sum(is.na(X))
X[is.na(X)] <- 0

# start simulation
traits_nulls <- data.frame(indiv = rownames(X))
h2_ranges <- seq(0,1,0.05)

for (h in h2_ranges[1:10]) {
  # handle extremes
  if (h == 0) h=0.01; if (h == 1) h=0.99
  h_out <- glue::glue("n{nsamp}_h{h}")
  print(glue::glue("Simulating {h_out}"))
  h = 0.9 # set this for specific value
  tmp_null <- gen_y_h2(1,X,h) %>% data.frame() %>%
    setnames(., names(.), c("null_y")) %>%
    rownames_to_column(var = "indiv")

  traits_nulls <- traits_nulls %>%
    left_join(tmp_null, by = "indiv") %>%
    rename(!! h_out := null_y)

}

# format name
traits_nulls <- traits_nulls %>%
  separate(indiv, c("FID","IID"), sep = "_") %>%
  select(-FID)

# save the simulated y
fwrite(traits_nulls,
       file = output_file,
       sep = "\t")

