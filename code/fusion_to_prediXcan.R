#! /usr/bin/env Rscript

library(RSQLite)
library(data.table)
library(tidyverse)

in.dir <- "/gpds/data" # update the paths accordingly

# parser function modified for gene name
make_df <- function(file) {
  load(file)  
  weights <- data.frame(wgt.matrix) 
  snps <- data.frame(snps) 
  rownames(weights) <- c() 
  gene <- str_extract(basename(file), "ENSG[0-9]+\\.[0-9]+")
  prefix <- paste0("GTEx\\.Whole_Blood\\.", gene, "\\.")
  suffix <- "\\.wgt\\.RDat"
  genename <- sub(prefix, "", basename(file))
  genename <- sub(suffix, "", genename)
  
  perf <- data.frame(cv.performance) %>% select(enet) %>% 
    rownames_to_column(var = "metric")
  
  weights$gene <- gene; weights$genename <- genename 
  weights$rsid <- snps$V2
  weights$chromosome <- snps$V1 
  weights$position <- snps$V4
  weights$ref_allele <- snps$V6
  weights$eff_allele <- snps$V5
  weights$pval <- perf %>% filter(metric == "pval") %>%  pull(enet)
  weights$R2 <- perf %>% filter(metric == "rsq") %>%  pull(enet)
  
  weights %>% filter(enet != 0) %>% #rownames_to_column(var = "rsid") %>% 
    select(gene, genename,rsid, chromosome, position, ref_allele, eff_allele,
           enet, pval, R2)
}

# all genes
files <- list.files(path = glue::glue("{in.dir}/ukb/paper-simulations/fusion_twas/WEIGHTS/GTEx.Whole_Blood/"), full.names = T, pattern = "\\.RDat")
# out file
pre_weights = glue::glue("{in.dir}/ukb/paper-simulations/fusion_twas/pre_weights.txt")

fess <- make_df(files[1])

write.table(make_df(files[1]), 
            pre_weights, sep = "\t", quote = FALSE, row.names = FALSE)
# loop through each chr
for(i in 2:length(files)) {
  write.table(make_df(files[i]), 
              pre_weights, append = TRUE, sep = "\t", quote = FALSE, 
              col.names = FALSE, row.names = FALSE)
}


# make the predixcan format database
pre_weights <- fread(glue::glue("{in.dir}/ukb/paper-simulations/fusion_twas/pre_weights.txt"))

weights <- pre_weights %>% 
  mutate(varID = glue::glue("chr{chromosome}_{position}_{ref_allele}_{eff_allele}_b37")) %>% 
  select(gene, rsid, ref_allele, eff_allele, weight = enet, varID)

extras <- pre_weights %>% 
  mutate(pred.perf.qval = NA) %>% 
  select(gene, genename, pred.perf.R2 = R2, 
         pred.perf.pval = pval, pred.perf.qval) %>% 
  distinct() %>% 
  inner_join(pre_weights %>% group_by(gene) %>% 
               summarise(n.snps.in.model = n()), by = "gene")

model_db = glue::glue("{in.dir}/ukb/paper-simulations/fusion_twas/WEIGHTS/GTEx.Whole_Blood.db")
conn <- dbConnect(RSQLite::SQLite(), model_db)
dbWriteTable(conn, "weights", weights)
dbWriteTable(conn, "extra", extras)
dbDisconnect(conn)

fwrite(extras %>% select(gene,genename),
       file = glue::glue("{in.dir}/ukb/paper-simulations/fusion_twas/gene_map_blood.txt"))

