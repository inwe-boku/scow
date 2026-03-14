bootstrap_stage2_project_library <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    return(invisible(NULL))
  }

  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
  lockfile <- file.path(script_dir, "renv.lock")
  if (!file.exists(lockfile)) {
    return(invisible(NULL))
  }

  if (!requireNamespace("renv", quietly = TRUE)) {
    stop("renv is required to load the tracked stage-2 package lockfile")
  }

  renv::load(project = script_dir, quiet = TRUE)
}

bootstrap_stage2_project_library()

library(logitr)
library(tidyverse)
library(broom)
library(parallel)

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
      slice_sample(n = n_choices, replace = FALSE) %>%
      mutate(choiceid = row_number()) %>%
      bind_rows(zone_in) %>% arrange(choiceid)
    logitr(
      data = X,
      outcome = "zoning",
      obsID = "choiceid",
      pars = base_mod
    )
}

main <- function() {
  if (!DEBUG){
    args <- commandArgs(trailingOnly = TRUE)

    last <- length(args)
    work_dir <- args[1]
    datafile <- args[2]
    num_runs <- as.integer(args[3])
    specif <- args[4]
    integer_count <- as.integer(args[5])
    seed <- as.integer(args[6])
    base_mod <- c(args[7:last])
    if (is.na(num_runs) || num_runs < 1L) {
      stop("num_runs must be a positive integer")
    }
    if (is.na(integer_count) || integer_count < 0L) {
      stop("integer_count must be a non-negative integer")
    }
    if (is.na(seed)) {
      stop("seed must be a valid integer")
    }
    if (integer_count > length(base_mod)) {
      stop("integer_count exceeds the number of provided model variables")
    }
    integer_model_vars <- if (integer_count == 0L) c() else base_mod[seq_len(integer_count)]
  }

  if (DEBUG){
    last <- 25
    work_dir <- '/Users/nwesec/code/scow'
    datafile <- '2014'
    num_runs <- 25L
    specif <- 'debug'
    seed <- 42L
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
    mutate(choiceid = row_number())

  zone_out <- geodat_noe %>%
    filter(zoning < 1) %>%
    mutate(zoning = as.integer(zoning))

  n_choices <- nrow(zone_in)
  if (nrow(zone_out) < n_choices) {
    stop(
      sprintf(
        "Cannot draw %d unzoned alternatives without replacement from only %d available rows.",
        n_choices,
        nrow(zone_out)
      )
    )
  }

  # data contains columns:
  # "choice" ... 0/1 indicating the outcome / chosen alternative,
  # obsID ... identifies each unique choice observation

  # Parallel estimation
  # It is recommended to use no more than 50% of cores.
  # For baseline generation, cap the worker count at 8.
  detected_cores <- detectCores()
  if (is.na(detected_cores)) {
    NMB_CORES <- 8L
  } else {
    NMB_CORES <- as.integer(min(8L, max(1L, floor(detected_cores / 2))))
  }
  # configure number of logit model runs
  NMB_RUNS <- seq_len(num_runs)

  cluster_debug <- DEBUG || identical(Sys.getenv("SCOW_R_DEBUG_CLUSTER"), "1")
  if (cluster_debug) {
    message("Cluster debug: Rscript=", Sys.which("Rscript"))
    message("Cluster debug: PATH=", Sys.getenv("PATH"))
    message(
      "Cluster debug: env=",
      paste(
        names(Sys.getenv(c("R_PARALLEL_PORT", "R_PARALLEL_BACKEND", "RSHCMD"))),
        Sys.getenv(c("R_PARALLEL_PORT", "R_PARALLEL_BACKEND", "RSHCMD")),
        sep = "=",
        collapse = " | "
      )
    )
    message("Cluster debug: NMB_CORES=", NMB_CORES)
    message("Cluster debug: master=127.0.0.1 port=default setup_strategy=sequential")
  }
  rscript_path <- Sys.which("Rscript")
  if (rscript_path == "") {
    stop("Rscript not found on PATH")
  }

  cl <- parallel::makePSOCKcluster(
    NMB_CORES,
    master = "127.0.0.1",
    setup_strategy = "sequential",
    outfile = "",
    rscript = rscript_path
  )
  on.exit(stopCluster(cl), add = TRUE)
  parallel::clusterSetRNGStream(cl, seed)
  clusterEvalQ(cl, library("logitr", quietly = TRUE))
  clusterEvalQ(cl, library("dplyr", quietly = TRUE, warn.conflicts = FALSE))
  clusterExport(cl, c("zone_out", "zone_in", "base_mod", "n_choices"), envir = environment())

  models <- parLapply(cl, NMB_RUNS, fun = do_estimate)

  coefs <- lapply(models, tidy, simplify = FALSE)
  coefs <- dplyr::bind_rows(coefs, .id = "mod")
  modstats <- lapply(models, glance) %>% dplyr::bind_rows(.id = "mod")

  # Restore colons in column names before export
  coefs$term <- restore_colons(coefs$term, colons_mapping)

  dir.create("data/results", recursive = TRUE, showWarnings = FALSE)
  write_csv(coefs, sprintf("data/results/spatialdc_coefs_%s_%s.csv", datafile, specif))
  write_csv(modstats, sprintf("data/results/spatialdc_loglik_%s_%s.csv", datafile, specif))
}

main()
