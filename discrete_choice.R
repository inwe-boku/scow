if(!require(logitr)){
    install.packages('logitr')
    library(logitr)
}
if(!require(tidyverse)){
    install.packages('tidyverse')
    library(tidyverse)
}
if(!require(broom)){
    install.packages('broom')
    library(broom)
}
if(!require(parallel)){
  install.packages('parallel')
  library(parallel)
}

args <- commandArgs(trailingOnly = TRUE)

last <- length(args)
work_dir <- args[1]
datafile <- args[2]
num_runs <- as.integer(args[3])
specif <- args[4]
last_integer <- as.integer(args[5])
integer_model_vars <- c(args[6:last_integer])

setwd(work_dir)

geodat_noe <- read_csv(sprintf('data/processed/dc_data_%s.csv', datafile))
geodat_noe <- mutate(geodat_noe, zoning = replace(zoning, zoning < 1, 0))
geodat_noe <- geodat_noe %>% mutate(across(any_of(integer_model_vars), as.integer))
if (all(c("tree_cover_density", "broadleaved") %in% colnames(geodat_noe))) {
    geodat_noe$"tree_cover_density:broadleaved" <- geodat_noe$tree_cover_density * geodat_noe$broadleaved
}
if (all(c("tree_cover_density", "coniferous") %in% colnames(geodat_noe))) {
    geodat_noe$"tree_cover_density:coniferous" <- geodat_noe$tree_cover_density * geodat_noe$coniferous
}
if (all(c("protected_areas", "important_bird_areas") %in% colnames(geodat_noe))) {
    geodat_noe$"protected_areas:important_bird_areas" <- geodat_noe$protected_areas * geodat_noe$important_bird_areas
}

base_mod <- c(args[6:last])

# drop variables which are not used to save memory
expl_variables_base <- unlist(str_split(base_mod, "\\*"))
geodat_noe <- geodat_noe[,c(expl_variables_base, "zoning")]

# Choice estimation
zone_in <- geodat_noe %>%
  filter(zoning > 0) %>%
  mutate(zoning = as.integer(zoning)) %>%
  mutate(alt = "A") %>%
  mutate(choiceid = row_number())

zone_out <- geodat_noe %>%
  filter(zoning < 1) %>%
  mutate(zoning = as.integer(zoning)) %>%
  mutate(alt = "B")

# data contains columns:
# "id" ... determines the individual,
# "alt" ... determines the alternatives included in the choice set of each observation,
# "choice" ... 0/1 indicating the outcome / chosen alternative,
# obsID ... identifies each unique choice observation

do_estimate <- function(j){
    X <- zone_out %>%
      slice_sample(n = sum(geodat_noe$zoning), replace = FALSE) %>%
      mutate(choiceid = row_number()) %>%
      bind_rows(zone_in) %>% arrange(choiceid)
    mlgt <- logitr(
      data = X,
      outcome = "zoning",
      obsID = "choiceid",
      pars = base_mod
    )
}

# Parallel estimation
# It is recommended to use no more than 50% of cores
NMB_CORES = detectCores() / 2
# configure number of logit model runs
NMB_RUNS <- 1:num_runs

cl <- makeCluster(NMB_CORES)
clusterExport(cl, c("geodat_noe", "zone_out", "zone_in"))
clusterEvalQ(cl, library("logitr", quietly = T))
clusterEvalQ(cl, library("tidyverse"))
clusterEvalQ(cl, library("broom"))

models <- list()
clusterExport(cl, c("base_mod"))
models <- parLapply(cl, NMB_RUNS, fun = do_estimate)

# save list of models
# saveRDS(models, file = sprintf("data/results/models_sig_%s_%s.rds", datafile, specif))
coefs <- lapply(models, tidy, simplify = F)
coefs <- dplyr::bind_rows(coefs, .id = "mod")
modstats <- lapply(models, glance) %>% dplyr::bind_rows(.id = "mod")
write_csv(coefs, sprintf("data/results/spatialdc_coefs_%s_%s.csv", datafile, gsub(":", "_", specif)))
write_csv(modstats, sprintf("data/results/spatialdc_loglik_%s_%s.csv", datafile, gsub(":", "_", specif)))

stopCluster(cl)
