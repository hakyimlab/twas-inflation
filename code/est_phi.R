#! /usr/bin/env R

library(tidyverse)
library(data.table)
library(RSQLite)
library(optparse)

option_list <- list(
  make_option(c("--sample_size"), type="character", default=1000,
              help="The sample size of null simulation",
              metavar="numeric"),
  make_option(c("--h2"), type="character", default=0.5,
              help="The heritability of null simulation",
              metavar="numeric"),
  make_option(c("--chi_table"), type="character", default=NULL,
              help="The estimated chi across simulations",
              metavar="character"),
  make_option(c("--input_db"), type="character", default=NULL,
              help="Original predixcan db to be updated",
              metavar="character"),
  make_option(c("--output_db"), type="character", default=NULL,
              help="The output db with phi",
              metavar="character")
  )

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)

phi_nsamp <- as.numeric(args$sample_size)
chi_table <- args$chi_table
in.db <- args$input_db
out.db <- args$output_db
phi_h2 <- as.numeric(args$h2)


# functions
se_calc <- function(x) sd(x)/sqrt(length(x) - 1)


# estimate phi
chi_df <- fread(chi_table)
gene_phi <- chi_df %>%
  summarise(across(everything(), list(
    avg.chi = ~mean(.),
    se = ~se_calc(.)
  ))) %>%
  pivot_longer(everything(),
               names_to = c("gene", "stat"), names_sep = "_") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(phi = (avg.chi - 1 + 2*se) / (phi_nsamp * phi_h2)) %>%
  mutate(phi = ifelse(phi < 0, 0, phi)) %>%
  mutate(se = ifelse(phi == 0, 0, se))


# load old db for updates
sqlite <- dbDriver("SQLite")
db = dbConnect(sqlite,in.db)
extras = dbGetQuery(db,"select * from extra")
weights = dbGetQuery(db,"select * from weights")
dbDisconnect(db)


# update the extras table
extras <- extras %>%
  inner_join(gene_phi %>% select(gene,phi), by = "gene")


# write out to a new file
driver <- dbDriver('SQLite')
conn <- dbConnect(drv = driver, out.db)
dbWriteTable(conn, 'extra', extras, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX metab_model_summary ON extra (gene)")

dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX weights_snpid ON weights (varID)")
dbExecute(conn, "CREATE INDEX weights_metabolite ON weights (gene)")
dbExecute(conn, "CREATE INDEX weights_snpid_metab ON weights (varID, gene)")

dbDisconnect(conn)


