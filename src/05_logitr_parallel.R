if(!require(gt)){
    install.packages('gt')
    library(gt)  # for generating latex tables
}
if(!require(gtsummary)){
    install.packages('gtsummary')
    library(gtsummary)  # for generating latex tables
}
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

setwd('c:/git_repos/scow')
sc_name <- c("cost_2016")  # "cost_2014",

for (sc in sc_name){

    geodat_noe <- read_csv(sprintf('data/preprocessed/geodat_Niederoesterreich_touch_%s.csv', sc))
    geodat_noe <- mutate(geodat_noe, zoning = replace(zoning, zoning < 1, 0))

    geodat_noe <- geodat_noe %>%
        mutate(across(c("zoning",
                  "existing_turbines",
                  "prt_areas",
                  "brd_areas",
                  "airports_buff",
                  "roads",
                  "waters",
                  "bauland",
                  "greenzoning",
                  "greenbuildings",
                  "broadleaved",
                  "coniferous",
                  "mil_areas",
                  "pastures",
                  "crops",
                  "wohnwidmung",
                  "preserve",
                  "lines",
                  "alpconv"), as.integer))
    geodat_noe$"tree_cover_density:broadleaved" <- geodat_noe$tree_cover_density * geodat_noe$broadleaved
    geodat_noe$"tree_cover_density:coniferous" <- geodat_noe$tree_cover_density * geodat_noe$coniferous

    base_mod <- c("airports_buff",
                  "brd_areas*prt_areas",
                  "elevation",
                  "mil_areas",
                  "overnights",
                  "pastures",
                  "prx_grid",
                  "prx_greenbuildings",
                  "prx_roads",
                  "prx_wohnwidmung",
                  "prx_bauland",
                  "prx_turbines",
                  "preserve",
                  "slope",
                  "alpconv",
                  "broadleaved",
                  "coniferous",
                  "tree_cover_density:broadleaved",
                  "tree_cover_density:coniferous",
                  "tree_cover_density",
                  "lcoe_min_lcoe"
    )
    full_mod <- c("prx_greenzoning",
                  "waters"
    )

    # drop variables which are not used to save memory
    expl_variables_base <- unlist(str_split(base_mod, "\\*"))
    expl_variables_full <- unlist(str_split(full_mod, "\\*"))
    geodat_noe <- geodat_noe[,c(expl_variables_base, expl_variables_full, "zoning")]

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
          pars = pars
        )
    }

    # Parallel estimation
    # It is recommended to use no more than 50% of cores
    NMB_CORES = detectCores() / 2
    # configure number of logit model runs
    NMB_RUNS <- 1:2500

    cl <- makeCluster(NMB_CORES)
    clusterExport(cl, c("geodat_noe", "zone_out", "zone_in"))
    clusterEvalQ(cl, library("logitr", quietly = T))
    clusterEvalQ(cl, library("tidyverse"))
    clusterEvalQ(cl, library("broom"))

    specifications <- c("base", "full")
    for (specif in specifications){
        models <- list()
        if (specif == "full"){
            pars <- c(base_mod, full_mod)
        } else {
            pars <- base_mod
        }
        clusterExport(cl, c("pars"))
        models <- parLapply(cl, NMB_RUNS, fun = do_estimate)
        # save list of models
        saveRDS(models, file = sprintf("data/results/models_sig_%s_%s.rds", sc, specif))  # s))
        coefs <- lapply(models, tidy, simplify = F)
        coefs <- dplyr::bind_rows(coefs, .id = "mod")
        modstats <- lapply(models, glance) %>% dplyr::bind_rows(.id = "mod")
        write_csv(coefs, sprintf("data/results/lgtr_coefs_%s_%s.csv", sc, specif))  # s))
        write_csv(modstats, sprintf("data/results/lgtr_loglik_%s_%s.csv", sc, specif))  # s))
    }

    stopCluster(cl)

}


#if(!require(stargazer)){
#    install.packages('stargazer')
#    library(stargazer)  # for generating latex tables
#}
#
#stargazer(as.data.frame(geodat_noe))
