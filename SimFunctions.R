library(tidyverse)
library(survival)
library(rsample)

# Simulate data -----

simfun <- function(N = 200, counterfactual = FALSE, 
                   tmeas = 3,
                   a_x_form = ~ X, a_x_coef = c(0,0),
                   y_form = ~ A*X + A*time, y_coef = c(0, 0,0, 0, 0),
                   c_form = ~ A*X+A*Y, c_coef = c(log(1), log(1), log(1), log(1), log(1)), c_basehaz = 0.1,
                   d_form = ~ A*X+A*Y, d_coef = c(log(1), log(1), log(1), log(1), log(1)), d_basehaz = 0.1,
                   fseed = "",
                   fdf = NULL){
  
  #browser()
  # check if arguments are properly filled in
  
  if( (length( labels(terms(a_x_form)) ) + 1) != length(a_x_coef) ){ stop("Dimensions of x-a formula and coefficient vector are not compatible (you need an intercept + one coefficient per covariate ) ") }
  
  if( length( labels(terms(y_form))) != length(y_coef)){ stop("For the outcome model, we assume the intercept to be 0. Length of y_coef should be equal to the number of predictors in the y_model")}
  
  if(c_form != ""){
    if(length(labels(terms(c_form))) != length(c_coef) ){"Dimension of c_form must match the number of coefficients in c_coef"}
  }
  if(d_form != ""){
    if(length(labels(terms(d_form))) != length(d_coef) ){"Dimension of d_form must match the number of coefficients in d_coef"}
  }
  
  if( fseed != ""){set.seed(fseed)}
  
  # trick to be able to reuse the same function for the counterfactuals
  if( is.null(fdf)){ # no df provided: new simulation
    
    Xsim <- runif(N, -1,1)
    RI <- rnorm(N)
    tempdf <- data.frame(X = Xsim, RI = RI) 
    
    x_p <-  myexpit(model.matrix(a_x_form, tempdf) %*% a_x_coef)
    tempdf$A <- rbinom(nrow(tempdf), 1, prob = x_p)
    
    tempdf$ID <- 1:nrow(tempdf)
    
    
    
  }else{  # df provided (counterfactual requested) => re-use baseline info
    
    tempdf <- fdf[c("A", "X", "ID", "RI")]
    tempdf$A <- 1-tempdf$A
    
  }
  
  tempdf <- expand_grid(tempdf,
                        time = 0:tmeas)
  
  # simulate outcome Y
  ps_y <- model.matrix(y_form, tempdf) %*% c(0,y_coef)
  tempdf$Y <- rnorm(nrow(tempdf), ps_y) + tempdf$RI
  
  
  # Censoring
  
  if(d_form == ""){d_basehaz <- 0}
  if(c_form == ""){c_basehaz <- 0}
  
  if(c_basehaz == 0){
    
    tempdf$Ctemp <- Inf
    
  }else{
    
    c_ps <- model.matrix(c_form, tempdf) %*% c(0, c_coef)
    
    tempdf$Ctemp <- rexp(nrow(tempdf), c_basehaz*exp(c_ps))
    
  }
  
  # Death 
  
  if(d_basehaz == 0){
    
    tempdf$Dtemp <- Inf
    
  }else{
    
    d_ps <- model.matrix(d_form, tempdf) %*% c(0, d_coef)
    tempdf$Dtemp <- rexp(nrow(tempdf), d_basehaz*exp(d_ps))
  }
  
  tempdf <- tempdf%>%
    mutate(Ctemp = ifelse(Ctemp < 1, Ctemp, Inf),
           Dtemp = ifelse(Dtemp < 1, Dtemp, Inf),
           C = time + Ctemp,
           D = time + Dtemp)%>%
    group_by(ID)%>%
    mutate(C = min(C),
           D = min(D))%>%
    mutate(Type = case_when(
      
      C > time & D > time ~ "Observed",
      C < time & C < D ~ "Censored",
      D < time ~ "Death",
      TRUE ~ "This option shouldn't exist"
      
    ),
    Yobs = case_when(
      
      Type == "Observed" ~ Y,
      TRUE ~ NaN
      
    ),
    ObsTime = pmin(C,D),
    Event = ifelse(ObsTime == D & !is.infinite(D), "Died", "Censored"),
    tstart = time,
    tstop = case_when(
      time+1 < ObsTime ~ time + 1,
      TRUE ~ ObsTime
    ),
    D_Ind = case_when(
      tstop == tstart +1 ~ FALSE,
      Event == "Died" ~ TRUE,
      TRUE ~ FALSE
    ),
    Death = (D < time))%>%
    mutate(simnr = fseed)
  
  
  if(counterfactual){ # if counterfactual, then redo
    
    cftemp <- simfun(fdf = tempdf%>%filter(time == 1), tmeas = tmeas ,counterfactual = FALSE,
                     a_x_form = a_x_form, a_x_coef = a_x_coef,
                     y_form = y_form, y_coef = y_coef,
                     c_form = c_form, c_coef = c_coef, c_basehaz = c_basehaz,
                     d_form = d_form, d_coef = d_coef, d_basehaz = d_basehaz,
                     fseed = fseed)
    
    
    tempdf$Y_cf <- cftemp$Y
    tempdf$Yobs_cf <- cftemp$Yobs
    tempdf$Type_cf <- cftemp$Type
    #tempdf$AdminCensTime_cf <- cftemp$AdminCensTime
    tempdf$CF_diff <- (-1)^(1-tempdf$A)*(tempdf$Y - tempdf$Y_cf)
    #tempdf$Type_cf <- cftemp$Type
    tempdf$D_Ind_cf <- cftemp$D_Ind
    tempdf$Death_cf <- cftemp$Death
    
    
  }
  
  
  tempdf%>%ungroup()%>%mutate(simnr = fseed)
  
}


## list with base-parameters for the simulation ----
## (update as needed for each new scenario)
baselist <- list(N= 200,
                  counterfactual = TRUE, 
                  tmeas = 1,
                  a_x_form = ~ X, a_x_coef = c(-1,2),
                  y_form = ~ A*X+A*time, y_coef = c(1, 1, 0.8),
                  c_form = ~ A*X+A*Y, c_coef = c(log(1),log(1),log(1),log(1),log(1)), c_basehaz = 0.1,
                  d_form = ~ A*X+A*Y, d_coef = c(log(1),log(1),log(1),log(1),log(1)), d_basehaz = 0.1,
                  fseed = "",
                  fdf = NULL)

### IPT ----

iptfun <- function(fdf, referencedf, form = ~X,
                   trunc = FALSE, truncp = 0.99){
  
  # fdf should be A = 0 or 1 in the sample. Must be weighted towards the referencedf
  
  # formula to be used
  tempform <- update(form, pop ~ .)
  
  # use WeightIt package to calculate weights to weight towards the reference
  tempW <- WeightIt::weightit(formula =  tempform,
                              data = bind_rows(fdf%>%mutate(pop = 0),
                                               referencedf%>%mutate(pop = 1)),
                              method = "glm",
                              estimand = "ATT")
  
  if(trunc){
    data.frame(ID = fdf$ID,
               W_ipt = WeightIt::trim(tempW$weights[1:nrow(fdf)]))
  }else{
    
    data.frame(ID = fdf$ID,
               W_ipt = tempW$weights[1:nrow(fdf)])  
    
  }
  
  
}



### IPC ----

ipcfun <- function(fdf, form = ~ X + Y, analysistimes){
  
  fdf <- fdf%>%filter(tstart <= max(analysistimes)+1)%>%
    filter(tstart < tstop)
  
  # change only last D_Ind for C_Ind (the actual C-event does not happen every intermediate tstop)
  
  fdf <- fdf%>%
    group_by(ID)%>%
    mutate(C_Ind = case_when(
      
      row_number() != max(row_number()) ~ D_Ind, # FALSE
      row_number() == max(row_number()) & tstop == C ~ ! D_Ind,
      TRUE ~ FALSE
      
    ))%>%ungroup()
  
  tempform <- update(form, Surv(tstart, tstop, C_Ind)~. )
  
  cfit <- coxph(tempform, data = fdf, id = ID, model = TRUE)
  
  res <- survfit(cfit, newdata = fdf, id = ID)
  
  # summary(res, times = analysistimes) would work as well
  
  pcdf <- data.frame(time = res$time,
                     prob = res$surv,
                     ID = addids(names(res$strata), res$strata))
  
  pcdf <- extractsurvps(pcdf, analysistimes)%>%
    mutate(W_ipc = 1/prob)
  
  pcdf
  
}


## Separate Longitudinal Analysis Functions ----

ieeiptcw <- function(dfy,   # longitudinal dataset
                     dfbase, # baseline dataset (to calculate ipt)
                     dfstand, # reference population
                     stand_form, # formula for iptw
                     analysistimes, # times at which Y needs to be analysed
                     type = c("iptcw", "ipconly", "iptonly", "crude"),
                     browse = FALSE){
  
  if(browse){browser()}
  #### Prep IPCW and IPTW  ----
  
  iptw <- c(0, 1)%>%map_df(function(x){
    iptfun(dfbase%>%filter(A == x), 
               referencedf = dfstand,
               form = stand_form, 
               trunc = FALSE,
               truncp = 0.99)%>%mutate(A = x)
  })

  
  fctemp <- dfy%>%left_join(iptw)%>%mutate(W_tot = W_ipc*W_ipt,
                                           W_crude = 1,
                                           W_ipconly = W_ipc,W_iptonly = W_ipt)%>%
    mutate(ctime = as.character(time))
  
  
  iptcres <- geepack::geeglm(Yobs ~ A*ctime,
                            weights = W_tot,
                            id = ID,
                            corstr = "independence",
                            data = fctemp)
  
  cruderes <- geepack::geeglm(Yobs ~ A*ctime,
                            weights = W_crude,
                            id = ID,
                            corstr = "independence",
                            data = fctemp)
  
  ipcres <- geepack::geeglm(Yobs ~ A*ctime,
                              weights = W_ipconly,
                              id = ID,
                              corstr = "independence",
                              data = fctemp)
  
  iptres <- geepack::geeglm(Yobs ~ A*ctime,
                            weights = W_iptonly,
                            id = ID,
                            corstr = "independence",
                            data = fctemp)
  
  out <- c("iptcres", "cruderes", "ipcres", "iptres")%>%
    map_df(function(res){
      
      fit <- get(res)
      
      tempests <- marginaleffects::comparisons(fit, variables = "A", by = "ctime")
      data.frame(Time = as.numeric(tempests$ctime),
                 Est = tempests$estimate,
                 SE = tempests$std.error)%>%mutate(Method = res)    
      
    })
  
  
out
  
}


ieeregstand <- function(dfy,
                        dfbase,
                        dfstand,
                        y_form,
                        surv_form,
                        tv_censoring = TRUE, 
                        analysistimes,
                        browse = FALSE){
  
  if(browse){browser()}
  
  dfy <- dfy%>%mutate(ctime = as.character(time))%>%
    left_join(dfbase%>%mutate(Ybase=Y)%>%select(ID, Ybase))%>%
    mutate(ID = as.character(ID))%>%
    filter(tstart < tstop)
  
  
  if(tv_censoring){
    
    # IPC weights calculated in separate step
    
  }else{dfy <- dfy%>%mutate(W_ipc = 1)} # Overwrite if not required
  
  df_stand_full <- expand_grid(ctime = as.character(analysistimes), 
                               dfstand%>%select(-A),
                               A = c(0,1))
  
  Yfit <-  geepack::geeglm(update(y_form, Yobs ~.),
                           id = ID,
                           corstr = "independence",
                           data = dfy%>%filter(time %in% analysistimes),
                           weights = W_ipc)
  
  ## For each observation in the target population sample: predicted Y
  ## (separately for each treatment)
  
  df_stand_full$ypred <- predict(Yfit, newdata = df_stand_full)
  
  
  survpart <- c(0,1)%>%map_df(  function(treat){
    
    survform <- update(surv_form, Surv(tstart, tstop, D_Ind) ~ .)
    
    # select data with A = treatarm
    tempdf <- dfy%>%filter(A == treat)
    
    sfit <- coxph(survform,
                  data = tempdf,
                  model = TRUE, weights = W_ipc)  
    # apply estimated survmodel to target population-sample
    tempsurv <- survfit(sfit, newdata = dfstand)
    
    tempsurv <- summary(tempsurv, times = analysistimes)
    tempsurv <- data.frame(cbind(tempsurv$surv, ctime = as.character(analysistimes)))  
    tempsurv%>%pivot_longer(-ctime, names_to = "delete", values_to = "psurv")%>%
      mutate(order = as.numeric(stringr::str_replace(delete, "X", "")))%>%
      arrange(order)%>%
      mutate(ID = rep(dfstand$ID, each = length(analysistimes)))%>%
      dplyr::select(ID, psurv, ctime)%>%mutate(A = treat)
    
  })
  
  
  ### Combine prepped target data with Y- and Survival predictions
  
  geeregstdata <- df_stand_full%>%
    left_join(survpart)%>%mutate(psurv = as.numeric(psurv))
  
  # lm on predicted outcomes, weighted by survival probabilities
  fittemp <- lm(ypred ~ A*ctime, data = geeregstdata, weights = psurv)
  esttemp <- marginaleffects::avg_comparisons(fittemp,
                                              variables = "A",
                                              by = "ctime",
                                               wts = "psurv")
  
  data.frame(Time = as.numeric(esttemp$ctime),
             Est = esttemp$estimate,
             SE = esttemp$std.error # better to bootstrap
  )%>%mutate(Method = "IEE-RegStand")
  
  
}

## Combined Longitudinal Analysis Function ----

analysis_long <- function(
    fdf,
    tmeas = 0:3,
    stand_type = "att",
    regres_par = list(stand_form = ~ (X+Ybase), # separately per A, so essentially ...*A
                      c_form = ~X+Y, # here too, essentially ...*A
                      y_form = ~ A*(ctime + X+Ybase),
                      surv_form = ~X+Ybase), # and here as well!!
    browse = FALSE,
    browseall = FALSE,
    robust = FALSE, robustN = 10,
    doboot = TRUE, nboot = 200
){
  
  if(browse){browser()}
  
  basefdf <- fdf%>%filter(time == 0)%>%select(- time)%>%mutate(Ybase = Y)%>%
    mutate(ID = as.character(ID))
  
  
  if(robust){
  # select only timepoints with enough observations  
  tsel <-  fdf%>%filter(!is.na(Yobs))%>%
      group_by(A, time)%>%
      summarise(N = n())%>%
      pivot_wider(names_from = "A", values_from = "N", names_prefix = "A")%>%
      filter(A1 > robustN & A0 > robustN & !is.na(A1) & !is.na(A0))%>%
      pull(time)
    
    
  fdf <- fdf%>%filter(time %in% tsel)
  
  tmeas <- tmeas[tmeas %in% tsel]
  
  }
  
  #### Target population ----
  
  if(stand_type == "ate"){
    
    df_stand <- basefdf
    
  }else if(stand_type == "att"){
    
    df_stand <- basefdf%>%filter(A == 1)
    
  }else if(stand_type == "atnt"){
    
    df_stand <- basefdf%>%filter(A == 0)
    
  }
  
  fdf <- fdf%>%filter(tstart < tstop & tstart <= max(tmeas))%>%mutate(ID = as.character(ID))
  
  fdf_orig <- fdf
  
  
  
  #### IPC (need it twice: for reg stand as well with tv censoring)----
  
  ipcw <- c(0, 1)%>%map_df(function(x){
    ipcfun(fdf%>%filter(A == x),
           analysistimes = tmeas,
           form = regres_par$c_form)%>%mutate(A = x)
  })
  
  
  fdf <- fdf%>%left_join(ipcw)
  
  
  #### IEE - IPCW -----
  ieeipct <- ieeiptcw(dfy = fdf, 
                     dfbase = basefdf,
                     dfstand = df_stand,
                     stand_form = regres_par$stand_form,
                     analysistimes = tmeas,
                     browse = browseall)
  
  
  ### Reg-Stand ----
  
  ieeregstand <- ieeregstand(dfy = fdf,
                             dfbase = basefdf,
                             dfstand = df_stand,
                             y_form = regres_par$y_form,
                             surv_form = regres_par$surv_form ,
                             analysistimes = tmeas,
                             tv_censoring = TRUE,
                             browse = browseall)
  
  
  ieeregstandnoipc <- ieeregstand(dfy = fdf,
                             dfbase = basefdf,
                             dfstand = df_stand,
                             y_form = regres_par$y_form,
                             surv_form = regres_par$surv_form ,
                             analysistimes = tmeas,
                             tv_censoring = FALSE,
                             browse = browseall)%>%
    mutate(Method = "IEE-Regstand-NoIPC")
    
  
  
  
  res <- ieeipct%>%bind_rows(ieeregstand)%>%bind_rows(ieeregstandnoipc)%>%
    mutate(simnr = unique(fdf$simnr))
  
  
  
  if(doboot){
    
    # %>%nest(-ID): one row per patient
    # to rsample::bootstraps
    # unnest, redefine dfbase (but not df_stand!)
    
    fdft <- fdf_orig%>%nest(data = -c(ID, A))
    
    bootsdf <- rsample::bootstraps(fdft,strata = A, times = nboot)
    
    
    
    resboot <- bootsdf$id%>%map_df(function(x){
      
       temp <- bootsdf%>%filter(id == x)%>%pull(splits)
       
       fdf <- temp[[1]]%>%as.data.frame()%>%mutate(ID = as.character(row_number()))%>%unnest(data)
       
       basefdf <- fdf%>%filter(time == 0)%>%select(- time)%>%mutate(Ybase = Y)%>%
         mutate(ID = as.character(ID))
       
       # now repeat all steps
       
       #### IPC (need it twice: for reg stand as well with tv censoring)----
       
       ipcw <- c(0, 1)%>%map_df(function(x){
         ipcfun(fdf%>%filter(A == x),
                analysistimes = tmeas,
                form = regres_par$c_form)%>%mutate(A = x)
       })
       
       
       fdf <- fdf%>%left_join(ipcw)
       
       
       #### IEE - IPCW -----
       ieeipct <- ieeiptcw(dfy = fdf, 
                           dfbase = basefdf,
                           dfstand = df_stand,
                           stand_form = regres_par$stand_form,
                           analysistimes = tmeas,
                           browse = browseall)
       
       
       ### Reg-Stand ----
       
       ieeregstand <- ieeregstand(dfy = fdf,
                                  dfbase = basefdf,
                                  dfstand = df_stand,
                                  y_form = regres_par$y_form,
                                  surv_form = regres_par$surv_form ,
                                  analysistimes = tmeas,
                                  tv_censoring = TRUE,
                                  browse = browseall)
       
       
       ieeregstandnoipc <- ieeregstand(dfy = fdf,
                                       dfbase = basefdf,
                                       dfstand = df_stand,
                                       y_form = regres_par$y_form,
                                       surv_form = regres_par$surv_form ,
                                       analysistimes = tmeas,
                                       tv_censoring = FALSE,
                                       browse = browseall)%>%
         mutate(Method = "IEE-Regstand-NoIPC")
       
       
       
       
       res <- ieeipct%>%bind_rows(ieeregstand)%>%bind_rows(ieeregstandnoipc)%>%
         mutate(simnr = unique(fdf$simnr),
                boot = x)
       
       res
      
    })
    
    
    
  seboot <- resboot%>%group_by(Time, Method)%>%
    left_join(res%>%rename(Or = Est)%>%select(-simnr, -SE))%>%
    summarise(Lqboot = quantile(Est-Or, prob = 0.025),
              Uqboot = quantile(Est-Or, prob = 0.975),
              Or = first(Or),
              LLperc = quantile(Est, prob = 0.025),
              ULperc = quantile(Est, prob = 0.975))%>%
    mutate(LLboot = Or - Uqboot,
           ULboot = Or - Lqboot
           )%>%select(-Lqboot, -Uqboot, -Or)
  
  res <- res%>%left_join(seboot)  
    
  }
  
}




# Miscellaneous ----


myexpit <- function(x){exp(x)/(1+exp(x))}
mylogit <- function(p){log(p/(1-p))}

addids <- function(labs, ns){
  # from survfit, combine strata with strata names to get to correct number of IDs repeated
  pmap(list(l = labs, n = ns), function(l,n){rep(l, n)})%>%unlist()
  
}

extractsurvps <- function(data, times){
  # extract the survival (censoring) probabilities from the survfit-based-result
  toadd <- expand.grid(ID = unique(data$ID),
                       time = times)%>%mutate(Extract = TRUE)
  
  out <- data%>%bind_rows(toadd)%>%
    group_by(ID)%>%
    arrange(time, .by_group = TRUE)%>%
    fill(prob, .direction = "down")%>%
    filter(Extract)%>%
    mutate(prob = ifelse(time == 0, 1, prob))%>%select(-Extract)
  
  out
  
  
}


realeff <- function(fdf){
  
  tempdf <- fdf%>%mutate(
    
    Y1 = case_when(
      
      A == 1 &  !Death ~ Y,
      A == 1 & Death ~ NaN,
      A == 0 & ! Death_cf ~ Y_cf,
      A == 0 & Death_cf ~ NaN
    ),
    
    Y0 = case_when(
      
      A == 0 &  !Death ~ Y,
      A == 0 & Death ~ NaN,
      A == 1 & ! Death_cf ~ Y_cf,
      A == 1 & Death_cf ~ NaN
      
    ))
  
  
  tempdf%>%group_by(A)%>%summarise(Eff = mean(Y1, na.rm=TRUE) - mean(Y0, na.rm=TRUE))%>%
    mutate(A = as.character(A))%>%
    bind_rows(
      
      
      tempdf%>%summarise(Eff = mean(Y1, na.rm=TRUE) - mean(Y0, na.rm=TRUE))%>%
        mutate(A = "Both")
      
      
    )
  
  
}


printsimlistfun <- function(fsimlist){
  
  
  fsimlist$y_present <- paste(
    paste( round(fsimlist$y_coef,2), attr(terms(fsimlist$y_form), "term.labels"), sep = " "),
    collapse = " + "
  )
  
  fsimlist$a_present <- paste(
    paste( round(fsimlist$a_x_coef,2), c("",attr(terms(fsimlist$a_x_form), "term.labels")), sep = " "),
    collapse = " + "
  )
  
  if(!is_formula(fsimlist$c_form)){
    fsimlist$c_coef <- "/"
    fsimlist$c_basehaz <- "/"
    fsimlist$c_form <- "/"
  }else{
    
    fsimlist$c_present <- paste(
      paste( round(fsimlist$c_coef,2), attr(terms(fsimlist$c_form), "term.labels"), sep = " "),
      collapse = " + "
    )
    
  }
  
  if(!is_formula(fsimlist$d_form)){
    fsimlist$d_coef <- "/"
    fsimlist$d_basehaz <- "/"
    fsimlist$d_form <- "/"
  }else{
    
    fsimlist$d_present <- paste(
      paste( round(fsimlist$d_coef,2), attr(terms(fsimlist$d_form), "term.labels"), sep = " "),
      collapse = " + "
    )
    
  }
  
  fsimlist
  
}




