# separate list for the 'large N truth'
largelist <- simlist
largelist$N <- Nlarge
largelist$fseed <- 20

# simulate large N truth-dataset
dflarge <- do.call(simfun, largelist)

# extract true treatment effects from large N-dataset
# (both counterfactuals are simulated)

truth <- 1:max(dflarge$time)%>%
  map_df( ~ realeff(dflarge%>%filter(time == .x))%>%mutate(time = .x))

# simulate the datasets
dfs <- 1:Nsim%>%map(function(isim){
  
  print(isim)
  tempsimlist <- simlist
  tempsimlist$fseed <- isim
  
  out <- do.call(simfun, tempsimlist)
  
})


# analyse the datasets
results <- dfs%>%map_df(function(fdf){
  print(unique(fdf$simnum))
  unique(fdf$time)%>%
    map_df(
      # each timepoint separately
      ~ analysis_singletp(fdf = fdf,
                        tmeas = .x)
                        
                        
      )%>%
    bind_rows(
      
      # longitudinal - categorical time - three way interaction for iee
      analysis_long(fdf, tmeas = 1:max(fdf$time),
                   regres_par = list(stand_form = ~ X,
                                                c_form = ~X,
                                                y_form = ~ A*(X + ctime), # X and time interact with A - correct model
                                                surv_form = ~X) )%>%mutate(Time = as.numeric(Time))%>%
        mutate(Method = ifelse(Method != "IEE-IPTCW", "IEE-Reg-Int-Full", Method))
      
    )%>%bind_rows(
      
      # longitudinal - categorical time - only two-way interactions with treatment
      analysis_long(fdf, tmeas = 1:max(fdf$time),
                    regres_par = list(stand_form = ~ X,
                                      c_form = ~ X,
                                      y_form = ~ X + A*ctime, # No X*A interaction
                                      surv_form = ~X) )%>%mutate(Time = as.numeric(Time))%>%
        filter(Method != "IEE-IPTCW")%>% #IPTCW already in previous calc
        mutate(Method = "IEE-Reg-Int-Partial")
      
    )

})




out <- list(dflarge = dflarge,
            truth = truth,
            dfs = dfs,
            results = results,
            sim = sim_id)


saveRDS(out, sprintf("Sim_%s.rds", sim_id))

rm(largelist, dflarge, truth, dfs, results, out )


