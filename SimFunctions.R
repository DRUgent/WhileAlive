library(tidyverse)
library(survival)

# Simulate data -----

simfun <- function(N = 200, counterfactual = FALSE, 
                   tmeas = 3,
                   a_x_form = ~ X, a_x_coef = c(0,0),
                   y_form = ~ A*X + A*time, y_coef = c(0, 0,0, 0, 0),
                   c_form = ~ A*X, c_coef = c(log(1), log(1), log(1)), c_basehaz = 0.1,
                   d_form = ~ A*X, d_coef = c(log(1), log(1), log(1)), d_basehaz = 0.1,
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
    tempdf <- data.frame(X = Xsim) 
    
    x_p <-  myexpit(model.matrix(a_x_form, tempdf) %*% a_x_coef)
    tempdf$A <- rbinom(nrow(tempdf), 1, prob = x_p)
    
    tempdf$ID <- 1:nrow(tempdf)
    
    
    
  }else{  # df provided (counterfactual requested) => re-use baseline info
    
    tempdf <- fdf[c("A", "X", "ID")]
    tempdf$A <- 1-tempdf$A
    
  }
  
  tempdf <- expand_grid(tempdf,
                        time = 1:tmeas)
  
  # simulate outcome Y
  ps_y <- model.matrix(y_form, tempdf) %*% c(0,y_coef)
  tempdf$Y <- rnorm(nrow(tempdf), ps_y) 
  
  basedf <- tempdf%>%filter(time == 1)%>%select(A, X, ID)
  
  # Censoring
  
  if(d_form == ""){d_basehaz <- 0}
  if(c_form == ""){c_basehaz <- 0}
  
  if(c_basehaz == 0){
    
    tempdf$C0 <- Inf
    
  }else{
    
    c_ps <- model.matrix(c_form, basedf) %*% c(0, c_coef)
    
    basedf$C <- rexp(nrow(basedf), c_basehaz*exp(c_ps))
    
    tempdf <- tempdf%>%left_join(basedf)
    basedf <- basedf%>%select(-C)
  }
  
  # Death 
  
  if(d_basehaz == 0){
    
    tempdf$D <- Inf
    
  }else{
    
    d_ps <- model.matrix(d_form, basedf) %*% c(0, d_coef)
    basedf$D <- rexp(nrow(basedf), d_basehaz*exp(d_ps))
    
    tempdf <- tempdf%>%left_join(basedf)
    
    basedf <- basedf%>%select(-D)
  }
  
  
  
  tempdf <- tempdf%>%
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
    AdminCensTime = case_when(
      
      Type == "Observed" ~ time,
      TRUE ~ pmin(C, D)
      
    ),
    Ev_Time = AdminCensTime,
    C_Ind = (Type == "Censored"),
    D_Ind = (Type == "Death"),
    Death = (D < time))%>%
    mutate(simnum = fseed)
  
  
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
    tempdf$AdminCensTime_cf <- cftemp$AdminCensTime
    tempdf$CF_diff <- (-1)^(1-tempdf$A)*(tempdf$Y - tempdf$Y_cf)
    tempdf$Type_cf <- cftemp$Type
    tempdf$D_Ind_cf <- cftemp$D_Ind
    tempdf$Death_cf <- cftemp$Death
    
    
  }
  
  
  tempdf
  
}


## list with base-parameters for the simulation ----
## (update as needed for each new scenario)
baselist <- list(N= 200,
                  counterfactual = TRUE, 
                  tmeas = 1,
                  a_x_form = ~ X, a_x_coef = c(-1,2),
                  y_form = ~ A*X, y_coef = c(1, 1, 0.8),
                  c_form = ~ A*X, c_coef = c(log(1),log(1),log(1)), c_basehaz = 0.1,
                  d_form = ~ A*X, d_coef = c(log(1),log(1),log(1)), d_basehaz = "",
                  fseed = "",
                  fdf = NULL)


# Analysis - Single Timpoint -----

## Building blocks ----

### IPT -----

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

ipcfun <- function(fdf, form = ~ X, analysistimes){
  
  # otherwise we have 'number time points' rows per individual in the survival-dataset
  # (applying different administrative censoring dependent on time point)
  fdf <- fdf%>%filter(time == max(time))
  
  # update formula (pec::ipcw requires the use of Event=1 , no code as if censoring is the event)
  tempform <- update(form, as.formula(Surv(Ev_Time, !C_Ind) ~.))
  
  # use pec-package to calculate censoring-probabilities
  restemp <- pec::ipcw(formula = tempform, 
                       data = fdf,
                       times = analysistimes,
                       method = "cox",
                       what = "IPCW.times"
  )
  
  # tidy results
  censprobs <- as.data.frame(restemp$IPCW.times)
  names(censprobs) <- analysistimes
  censprobs <- censprobs%>%pivot_longer(cols = everything(), names_to = "time", values_to = "censprobs")
  
  # Put results in a dataframe
  data.frame(ID = rep(fdf$ID, each = length(analysistimes)),
             W_ipc = 1/censprobs$censprobs,
             time = as.numeric(censprobs$time))
  
}


## Separate single timepoint analyses ----

### IPTCW ----

iptcresfun <- function(fctemp){
  # would make sense to have the iptfun and ipcfun inside this function,
  # but in this simulation it is more efficient to have it outside,
  # as the weights are also used for iptw-only and ipcw-only
  
  fittemp <- glm(Yobs ~ A, data = fctemp, weights = W_tot)
  esttemp <- marginaleffects::avg_comparisons(fittemp,
                                              variables = "A",
                                              vcov = "HC3",
                                              wts = "W_tot")
  
  
  data.frame(Method = "IPTCW",
             Est = esttemp$estimate,
             SE = esttemp$std.error)
  
}

### Regression Standardization ---

regstandresfun <- function(ffctemp,
                            dfstand,
                            y_form = ~ X,
                           d_form = ~X){
  
  analysistime <- unique(ffctemp$time)
  # separately per treatment
  currtreat <- unique(ffctemp$A)
  # formula for Y
  tempqolform <- update(y_form, as.formula(Yobs ~ .))
  # formula for survival
  # AdminCensTime: 
  #    if death before time of interest: death time
  #    if censored before time of interest: censoring time
  #    no event nor censored before time of interest, censored just after time of interest (here t+1)
  #    So essentially observed time if study ran only up till 'time of interest'
  tempSform <- update(d_form, as.formula(Surv(AdminCensTime, D_Ind) ~ .))
  
  regstanddf <- currtreat%>%map_df(function(x){
    
    tempx <- ffctemp%>%filter(A == x)
  
  # fit Y
  qfit <- lm(tempqolform, data = tempx, model = TRUE)
  # fit survival
  sfit <- coxph(tempSform, data = tempx, model = TRUE)
  

    tempstand <- dfstand%>%mutate(A = x)  # set treatment of target pop to current treatment
    tempstand$Yest <- predict(qfit, newdata = tempstand, na.action = na.pass)  # estimate/predict qol while alive
    
    if(any(ffctemp$D_Ind == 1)){ # estimate survival if any events
      
      tempsurv <- survfit(sfit, newdata = tempstand)
      
      tempsurv <- summary(tempsurv, times = analysistime)
      
      
      tempstand$psurv <- as.vector(tempsurv$surv)
      
    }else{ # if no events, then no use (and potentially problems)
          # to fit a survival model      
      tempstand$psurv <- 1
      
    }
    
    tempstand%>%select(A, Yest, psurv)
    
  })
  
  
  fit <- lm(Yest ~ A, data = regstanddf, weights = psurv)
  
  data.frame(Method = "Regression-Standardization",
             Est = coef(summary(fit))[2,1],
             SE = coef(summary(fit))[2,2]) # better to bootstrap
  
 
  
}


### Crude differences ----


cruderesfun <- function(fctemp){
  
  fit_raw <- summary(lm(Yobs ~ A, data = fctemp))
  
  data.frame(Method = "Crude",
             Est = fit_raw$coefficients[2,1],
             SE = fit_raw$coefficients[2,2])
}


### IPC only ----

ipcresfun <- function(fctemp){
  
  fittemp <- summary(lm(Yobs ~ A, data = fctemp, weights = W_ipc))
  
  data.frame(Method = "IPC",
             Est = fittemp$coefficients[2,1],
             SE = fittemp$coefficients[2,2])
}

### IPT only ----

iptresfun <- function(fctemp){
  
  fittemp <- summary(lm(Yobs ~ A, data = fctemp, weights = W_ipt))
  
  data.frame(Method = "IPT",
             Est = fittemp$coefficients[2,1],
             SE = fittemp$coefficients[2,2])
}





## Combined single timepoint Analyses ----
# (crude, iptw only, ipcw only, 
#      iptcw, regression standardization )
analysis_singletp <- function(
    fdf,
    tmeas = 1,
    stand_type = "att",
    uncensored = TRUE,
    add_regres = TRUE,
    regres_par = list(stand_form = ~ X,
                      c_form = ~ X,
                      y_form = ~ X,
                      d_form = ~ X)
){
  
  #### ST-FUN Data Prep ----
  
  fdf <- fdf%>%filter(time == tmeas)
  
  #### ST-FUN Target population ----
  
  if(stand_type == "ate"){
    
    df_stand <- fdf
    
  }else if(stand_type == "att"){
    
    df_stand <- fdf%>%filter(A == 1)
    
  }else if(stand_type == "atnt"){
    
    df_stand <- fdf%>%filter(A == 0)
    
  }
  
  #### ST-FUN Prep IPCW and IPTW  ----
  # steps could be part of the individual iptcw-function,
  # but in this simulation, it is more efficient to have it outside,
  # as the weights are also used for iptw-only and ipcw-only
  
  iptw <- c(0, 1)%>%map_df(function(x){
    iptfun(fdf%>%filter(A == x), 
               referencedf = df_stand,
               form = regres_par$stand_form, 
               trunc = FALSE,
               truncp = 0.99)%>%mutate(A = x)
  })
  
  
  ipcw <- c(0, 1)%>%map_df(function(x, analysistimes = tmeas, form = regres_par$c_form){
    ipcfun(fdf%>%filter(A == x),
               form = form,
           analysistimes = tmeas)%>%mutate(A = x)
  })
  
  
  iptot <- iptw%>%left_join(ipcw)%>%mutate(W_tot = W_ipc*W_ipt)
  
  fctemp <- fdf%>%left_join(iptot)
  
  
  
  ## ST-FUN IPW Analyses ----
  
  res <- cruderesfun(fctemp)%>%
    bind_rows(
      
      ipcresfun(fctemp)
    )%>%
    bind_rows(
      iptresfun(fctemp)
    )%>%
    bind_rows(
      iptcresfun(fctemp)
      
    )
  
  ## ST-FUN Regression Standardization  ----
  
  regstand <- regstandresfun(fctemp, dfstand = df_stand, y_form = regres_par$y_form,
                             d_form = regres_par$d_form)

  
  res <- res%>%bind_rows(regstand)%>%mutate(Time = tmeas)
  rownames(res) <- NULL
  
  
  res
  
}



## Separate Longitudinal Analysis Functions ----

ieeiptcw <- function(dfy,   # longitudinal dataset
                     dfbase, # baseline dataset (for to calculate ipt)
                     dfstand, # reference population
                     stand_form, # formula for iptw
                     c_form, # formula for ipcw
                     analysistimes # times at which Y needs to be analyses
                     ){
  
  #browser()
  #### Prep IPCW and IPTW  ----
  
  iptw <- c(0, 1)%>%map_df(function(x){
    iptfun(dfbase%>%filter(A == x), 
               referencedf = dfstand,
               form = stand_form, 
               trunc = FALSE,
               truncp = 0.99)%>%mutate(A = x)
  })
  
  
  ipcw <- c(0, 1)%>%map_df(function(x){
    ipcfun(dfy%>%filter(A == x),
                   analysistimes = analysistimes,
                   form = c_form)%>%mutate(A = x)
  })
  
  
  iptot <- ipcw%>%left_join(iptw)%>%mutate(W_tot = W_ipc*W_ipt)
  
  fctemp <- dfy%>%left_join(iptot)%>%mutate(ctime = as.character(time))
  
  
  geeres <- geepack::geeglm(Yobs ~ A*ctime,
                            weights = W_tot,
                            id = ID,
                            corstr = "independence",
                            data = fctemp)
  
  tempests <- marginaleffects::comparisons(geeres, variables = "A", by = "ctime")
  data.frame(Time = as.numeric(tempests$ctime),
             Est = tempests$estimate,
             SE = tempests$std.error)%>%mutate(Method = "IEE-IPTCW")

  
}


ieeregstand <- function(dfy,
                        dfbase,
                        dfstand,
                        y_form,
                        surv_form,
                        analysistimes){
  
  #browser()
  dfy <- dfy%>%filter(time %in% analysistimes)%>%mutate(ctime = as.character(time))
  
  df_stand_full <- expand_grid(ctime = as.character(analysistimes), 
                               dfstand%>%select(-A),
                               A = c(0,1))
  
  Yfit <-  geepack::geeglm(update(y_form, Yobs ~.),
                           id = ID,
                           corstr = "independence",
                           data = dfy)
  
  ## For each observation in the target population sample: predicted Y
  ## (separately for each treatment)
  
  df_stand_full$ypred <- predict(Yfit, newdata = df_stand_full)
  
  
  survpart <- c(0,1)%>%map_df(  function(treat){
    
    survform <- update(surv_form, Surv(Ev_Time, D_Ind) ~ .)
    
    # select data with A = treatarm
    tempdf <- dfy%>%filter(time == max(time))%>%filter(A == treat)
    
    sfit <- coxph(survform,
                  data = tempdf,
                  model = TRUE)  
    # apply estimated survmodel to target population-sample
    tempsurv <- survfit(sfit, newdata = dfstand)
    # extract estimated survival probabilities (if tweeked - could have used pec::ipcw as well)
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
                                              vcov = "HC3",
                                              wts = "psurv")
  
  data.frame(Time = as.numeric(esttemp$ctime),
             Est = esttemp$estimate,
             SE = esttemp$std.error # better to bootstrap
  )%>%mutate(Method = "IEE-RegStand")
  
  
}

## Combined Longitudinal Analysis Function ----

analysis_long <- function(
    fdf,
    tmeas = 1:3,
    stand_type = "att",
    regres_par = list(stand_form = ~ X*A,
                      c_form = ~X,
                      y_form = ~ A*(ctime + X),
                      surv_form = ~X)
){
  
  #browser()
  
  basefdf <- fdf%>%filter(time == 1)%>%select(- time)
  
  
  #### Target population ----
  
  if(stand_type == "ate"){
    
    df_stand <- basefdf
    
  }else if(stand_type == "att"){
    
    df_stand <- basefdf%>%filter(A == 1)
    
  }else if(stand_type == "atnt"){
    
    df_stand <- basefdf%>%filter(A == 0)
    
  }
  
  fdf <- fdf%>%filter(time %in% tmeas)
  
  ieeipct <- ieeiptcw(dfy = fdf, 
                     dfbase = basefdf,
                     dfstand = df_stand,
                     stand_form = regres_par$stand_form,
                     c_form = regres_par$c_form,
                     analysistimes = tmeas)
  
  
  ### Reg-Stand ----
  
  ieeregstand <- ieeregstand(dfy = fdf,
                             dfbase = basefdf,
                             dfstand = df_stand,
                             y_form = regres_par$y_form,
                             surv_form = regres_par$surv_form ,
                             analysistimes = tmeas)
  
  
  
  
  ieeipct%>%bind_rows(ieeregstand)
  
}


## Update with Example-Analysis Code
# IEE-IPCW and IEE-Regression Standardization


analyse_long <- function(
    fdf,
    tmeas = 1:3,
    stand_type = "att",
    regres_par = list(stand_form = ~ X*A,
                      y_form = ~ X*A,
                      c_form = ~X, 
                      d_form = ~ X*A,
                      iee_form = ~ X + A*ctime,
                      surv_form = ~X)
){
  
  
  basefdf <- fdf%>%filter(time == 1)%>%select(- time)

  
  #### Target population ----
  
  if(stand_type == "ate"){
    
    df_stand <- basefdf
    
  }else if(stand_type == "att"){
    
    df_stand <- basefdf%>%filter(A == 1)
    
  }else if(stand_type == "atnt"){
    
    df_stand <- basefdf%>%filter(A == 0)
    
  }
  
  fdf <- fdf%>%filter(time %in% tmeas)

  geeRES <- ieeiptcw(dfy = fdf, 
                     dfbase = dfbase, dfstand = dfstand,
                     stand_form = regres_par$stand_form,
                     c_form = regres_par$c_form,
                     analysistimes = tmeas)
  
  
  ### Reg-Stand ----
  
  
  

  
  df_stand_full$ypred <- predict(Yfit, newdata = df_stand_full)
  
  survpart <- 0:1%>%map_df(function(x){
    
    survfunlong_sim(fdf%>%filter(A == x),
                    df_stand = df_stand,
                    tmeas = tmeas,
                    form = regres_par$surv_form)%>%mutate(A = x)
    
  })
  
  geeregst_tot <- df_stand_full%>%left_join(survpart%>%mutate(ctime = as.character(time)))
  geeres2 <- geeregst_tot%>%group_by(ctime, A)%>%
    summarise(M = weighted.mean(ypred, Psurv))%>%
    summarise(Est = diff(M))%>%
    dplyr::rename(Time = ctime)%>%
    mutate(Method = "GEE-Reg")
  
  
  geeRES%>%bind_rows(geeres2)%>%bind_rows(MMRES)
  
  
  
}  


# Miscellaneous ----


myexpit <- function(x){exp(x)/(1+exp(x))}
mylogit <- function(p){log(p/(1-p))}

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




