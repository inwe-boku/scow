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

DEBUG <- FALSE

# Function to replace colons with underscores and track original mappings
replace_colons <- function(strings) {
  modified <- gsub(":", "_", strings)
  setNames(modified, strings)  # Create a mapping from original to modified
}

# Function to restore colons using the original-to-modified mapping
restore_colons <- function(values, mapping) {
  sapply(values, function(v) {
    if (v %in% mapping) names(mapping)[mapping == v] else v
  })
}

# Function to compute interaction terms
compute_interaction <- function(term_name, data) {
  term_name <- as.character(term_name)  # Ensure term_name is a character
  vars_needed <- strsplit(term_name, ":")[[1]]
  if (all(vars_needed %in% names(data))) {
    data[[term_name]] <- data[[vars_needed[1]]] * data[[vars_needed[2]]]
  }
  return(data)  # Explicitly return the data frame
}

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

if (!DEBUG){
  args <- commandArgs(trailingOnly = TRUE)

  last <- length(args)
  work_dir <- args[1]
  datafile <- args[2]
  num_runs <- as.integer(args[3])
  specif <- args[4]
  last_integer <- as.integer(args[5])
  integer_model_vars <- c(args[6:last_integer])
  base_mod <- c(args[6:last])
}

if (DEBUG){
  last <- 25
  work_dir <- '/Users/nwesec/repos/scow'
  datafile <- '2014'
  num_runs <- 25
  specif <- 'debug'
  last_integer <- 11
  integer_model_vars <- c(
    "airports_buff",
    "alps_convention",
    "broadleaved",
    #"coniferous",
    "important_bird_areas",
    "protected_areas:important_bird_areas",
    "pastures",
    "preservation",
    "protected_areas",
    "restricted_military_areas",
    "water_bodies"
  )
  base_mod <- append(integer_model_vars, c(
    "distance_buildings_in_greenland",
    "distance_existing_turbines",
    "distance_greenland_zonings",
    "distance_other_building_land",
    "distance_power_lines",
    "distance_residential_buildings",
    "distance_roads",
    "elevation",
    "min_lcoe",
    "overnight_stays",
    "slope",
    "tree_cover_density",
    "tree_cover_density:broadleaved",
    "tree_cover_density:coniferous"))
}

setwd(work_dir)

geodat_noe <- read_csv(sprintf('data/processed/dc_data_%s.csv', datafile))
geodat_noe <- mutate(geodat_noe, zoning = replace(zoning, zoning < 1, 0))
geodat_noe <- geodat_noe %>% mutate(across(any_of(integer_model_vars), as.integer))

interaction_terms <- base_mod[grepl(":", base_mod)]
colons_mapping <- replace_colons(base_mod)
# Identify interaction terms (those with ":") in provided_vars
base_mod <- unname(colons_mapping)

# Correct Reduce usage
geodat_noe <- Reduce(
  function(data, term) compute_interaction(term, data),
  interaction_terms,
  init = geodat_noe
)

colnames(geodat_noe) <- gsub(":", "_", colnames(geodat_noe))
# drop variables which are not used to save memory
# expl_variables_base <- unlist(str_split(base_mod, "\\*"))
geodat_noe <- geodat_noe[,c(base_mod, "zoning")]

# replace : with _ in variable names
#base_mod <- gsub(":", "_", base_mod)

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

# Parallel estimation
# It is recommended to use no more than 50% of cores
NMB_CORES = detectCores() / 2
if (is.na(NMB_CORES)){
    NMB_CORES <- 8
}
# configure number of logit model runs
NMB_RUNS <- 1:num_runs

cl <- parallel::makeCluster(as.integer(NMB_CORES), type = 'FORK')
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

# Restore colons in column names before export
coefs$term <- restore_colons(coefs$term, colons_mapping)

write_csv(coefs, sprintf("data/results/spatialdc_coefs_%s_%s.csv", datafile, gsub(":", "_", specif)))
write_csv(modstats, sprintf("data/results/spatialdc_loglik_%s_%s.csv", datafile, gsub(":", "_", specif)))

stopCluster(cl)
