library(furrr)
# separate list for the 'large N truth'
largelist <- simlist
largelist$N <- Nlarge
largelist$fseed <- 20

# simulate large N truth-dataset
dflarge <- do.call(simfun, largelist)

# extract true treatment effects from large N-dataset
# (both counterfactuals are simulated)

truth <- 0:max(dflarge$time)%>%
  map_df( ~ realeff(dflarge%>%filter(time == .x))%>%mutate(time = .x))


saveRDS(truth, "truthseparate.rds")
saveRDS(dflarge, "dflargesep.rds")

# simulate the datasets
dfs <- 1:Nsim%>%map(function(isim){
  
  print(isim)
  tempsimlist <- simlist
  tempsimlist$fseed <- isim
  
  out <- do.call(simfun, tempsimlist)
  
})

plan(multisession)
Ns <- nbrOfWorkers()
Ns
plan(multisession, workers = Ns - 2)
nbrOfWorkers()

# analyse the datasets
results <- dfs%>%future_map_dfr(function(fdf){
  
    #fdf <- dfs[[x]]
    # bad idea to have dfs within this function: is then transfered at each session!!
    print(unique(fdf$simnr))
    x <- unique(fdf$simnr)
      # longitudinal - categorical time - three way interaction for iee
    tempres <- analysis_long(fdf, tmeas = 0:max(fdf$time),
                    regres_par = list(stand_form = ~ X, # separately per A, so essentially ...*A
                                      c_form = ~X, # here too, essentially ...*A
                                      y_form = ~ A*(ctime+X),
                                      surv_form = ~X), # and here as well!!
                    robust = TRUE, browse = FALSE,
                    nboot = Nboot)%>%mutate(Time = as.numeric(Time))%>%
        mutate(Setting = "Baseline Censoring Fitted only")%>%
        bind_rows(
      
      # longitudinal - categorical time - only two-way interactions with treatment
      analysis_long(fdf, tmeas = 0:max(fdf$time),
                    regres_par = list(stand_form = ~ X,
                                      c_form = ~  X+Y,
                                      y_form = ~ A*(ctime + X), 
                                      surv_form = ~X), robust = TRUE,
                    nboot = Nboot )%>%mutate(Time = as.numeric(Time))%>%
      mutate(Setting = "TV censoring fitted")
      
    )%>%bind_rows(
      
          
      analysis_long(fdf, tmeas = 0:max(fdf$time),
                    regres_par = list(stand_form = ~ X,
                                      c_form = ~ X + Y,
                                      y_form = ~ X + A*(ctime), 
                                      surv_form = ~X), robust = TRUE,
                    nboot = Nboot )%>%mutate(Time = as.numeric(Time))%>%
        mutate(Setting = "Wrong Y-model")
      
      
    )
    
    saveRDS(tempres, glue::glue("tempres{x}.rds"))

},  .options = furrr_options(seed = 102030))


rm(largelist, dflarge, truth, dfs, results)
