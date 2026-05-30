#note that version 2 of the function is modified for using the full prediction data and removes "indivPred_SAEM"
mAb_conc_protect_efficacy_plot <- function(monolix_pred_df, conc_type = "observed", d_before_infect = 10, startdate = "day_of_infusion", efficacy_threshold = 0.80, nboot = 1000, sigfigs = 2){
    if (!conc_type %in% c("observed", "indivPredMean", "indivPredMode", "popPred") == TRUE) {
      stop('conc_type must be "observed", "indivPredMean", "indivPredMode" or "popPred".')
  } else {
    message("installing dependencies")
    library(pacman)
    p_load(googledrive)
    p_load(tidyverse)
    p_load(survival)
    #change ID to lowercase
    names(monolix_pred_df) <- if_else(names(monolix_pred_df)=="ID", tolower(names(monolix_pred_df)), names(monolix_pred_df))
    if(conc_type == "observed"){
      conc_var <- "DV"
      monolix_pred_df <- monolix_pred_df %>%
        drop_na(DV)
    }
    else{
      conc_var <- conc_type
    }
    }
  if (any(!class(monolix_pred_df) %in% c("data.frame", "tbl", "spec_tbl_df", "tbl_df"))) {
    stop('monolix_pred_df must be a dataframe or tibble')
  } else {
    
    ## import time to event data
    temp<- tempfile(fileext = ".csv")
    local_drive_quiet()
    dl <- drive_download(
      as_id("1zpog1PWetOjKC46WeNyPJF0LNXUONiJq"), path = temp, overwrite = TRUE)
    dat_tf <- as.data.frame(read.csv(file = dl$local_path)[,-1])
    ## import CIS43LS study visit data
    temp<- tempfile(fileext = ".csv")
    dl <- drive_download(
      as_id("1zqZ4T1y6W4YImXldh9wfldFwBh3s9I8i"), path = temp, overwrite = TRUE)
    cis43ls_visits_dat <- as.data.frame(read.csv(file = dl$local_path)[,-1]) %>%
      mutate(date = as.Date(date, "%Y-%m-%d"))
    ## arrange monolix output and add interpolated, monolix time assumed to be in hours
    PK.matched <- monolix_pred_df %>%
      dplyr::select(id, time, DV, indivPredMode, indivPredMean, popPred) %>%
      ## key line to assign Conc from either observed, individual predicted, or pop predicted column
      dplyr::rename(Conc = conc_var) %>% 
      mutate(day = floor(time/24)) %>% #convert to days
      left_join(., cis43ls_visits_dat,
                by = c("id", "day")) %>%
      dplyr::select(id, date, date_type, day, time, everything()) %>%
      ## keep date as original date of infection to later merge with dat.tf using tmerge function)
      group_by(id) %>%
      mutate(date_at_days_before_for_conc = date - d_before_infect) %>%
      mutate(day_at_days_before_for_conc = day - d_before_infect) %>%
      mutate(time_at_days_before_for_conc = time - d_before_infect*24) %>%
      mutate(admin_date = first(date)) %>%
      ungroup() %>%
      arrange(id, day) %>%
      distinct() %>% #remove duplicates
      ## remove d140 time points for Subject 1401 and 1681 as these are aberrant outliers
      filter(!(day == 140 & id == 1404)) %>% 
      filter(!(day == 140 & id == 1681))
      if(conc_type == "observed"){
        PK.matched <- PK.matched %>%
          bind_rows(filter(., day > d_before_infect) %>% 
                      mutate(date_type = "PK_interprolated") %>%
                      mutate(day_at_days_before_for_conc = day - d_before_infect) %>% ## line to subtract d_before_infect days from study day
                      mutate(date_at_days_before_for_conc = date - d_before_infect) %>% ## line to subtract d_before_infect days from study day
                      mutate(time_at_days_before_for_conc = time - d_before_infect*24) %>% ## line to subtract 240 from time
                      mutate(Conc = NA))  %>%
          mutate(date_type = ifelse(!is.na(Conc), "PK_obs_fit", date_type)) %>%
          arrange(id, time) %>% 
          mutate(Conc_interp = NA) %>%
          group_by(id) %>%
          mutate(time = ifelse(is.na(time), as.numeric(Date-first(Date)), time)) %>%
          mutate(Conc_interp = ifelse(is.na(Conc),
                                      (((lead(Conc)-lag(Conc))/(lead(day)-lag(day)))*(day-lag(day)))+lag(Conc),
                                      Conc)) %>% #estimate the concentration by multiplying slope by time in question then adding the y intercept
          mutate(Conc_interp = ifelse(is.na(Conc_interp),
                                      (((lag(Conc)-lag(Conc,n=2L))/(lag(day)-lag(day,n=2L)))*(day-lag(day,n=2L)))+lag(Conc,n=2L),
                                      Conc_interp)) %>% #for rare instances where there are the malaria event is the last visit and there is no actual conc.
          ungroup() %>%
          dplyr::select(id, date_type, admin_date, date, day, time, contains("days_before_for_conc"), everything()) %>%
          arrange(id, date, Conc) %>%
          drop_na(Conc)
      } else {
        PK.matched <- PK.matched %>%
          arrange(id, time) %>% 
          dplyr::select(id, date_type, admin_date, date, day, time, everything()) %>%
          arrange(id, date, Conc) %>%
          drop_na(Conc)
      }
    
    # plot to double check
    # PK.matched %>%
    #   ggplot(., aes(x = day, y = Conc_interp, group = id)) +
    #   geom_line(aes(color = id)) +
    #   scale_color_viridis_c() +
    #   theme_bw()
    
  }
  if (any(!startdate %in% c("study_start", "day_of_infusion")) == TRUE){
    stop('startdate must be either "study_start" or "day_of_infusion"')
    } else {
      #### time should use PK.matched$date (original date of events [infections])
      #calendar day of evaluation (1st date of administration was 2021-05-05)
      if(startdate == "study_start"){
        time <- as.numeric(difftime(as.character(PK.matched$date), "2021-05-05", units="days"))
        }
      #date of administration for each individual
      if(startdate == "day_of_infusion"){
        time <- as.numeric(difftime(as.character(PK.matched$date), as.character(PK.matched$admin_date), units="days")) 
        }
      #Remove subjects who had an infection within the first 7 days and subtract "d_before_infect" from postdat or postdat.0
      dat_tf <- dat_tf %>%
        filter(!(infection == 1 & posdat <= 7)) %>%
        mutate(posdat = ifelse(infection == 1 & d_before_infect > 0,
                                posdat - d_before_infect, posdat)) %>%
        mutate(posdat.0 = ifelse(infection == 1 & d_before_infect > 0,
                                 posdat.0 - d_before_infect, posdat.0))
      #only use interpolated concentrations if using observed concentrations
      if(conc_type == "observed"){
        df.td = data.frame(id=PK.matched$id, time=time, conc = PK.matched$Conc_interp, StudyDay = PK.matched$day)
      }
      else{
        df.td = data.frame(id=PK.matched$id, time=time, conc = PK.matched$Conc, StudyDay = PK.matched$day)
      }
      
      
      #####################################
      ### keep posdat or posdat.0 unchanged
      #####################################
      #calendar day of evaluation (1st date of administration was 2021-05-05)
      if(startdate == "study_start"){
        df <- tmerge(dat_tf, dat_tf, id=id, endpt=event(posdat.0,infection)) 
        }
      #date of administration for each individual
      if(startdate == "day_of_infusion"){
        df <- tmerge(dat_tf, dat_tf, id=id, endpt=event(posdat,infection))  #admin 
        }
      df.m <- tmerge(df, df.td, id=id, conc=tdc(time,conc), StudyDay=tdc(time, StudyDay)) #adding time dependent covariates by merging time with true study day
      df.m$conc.log = log(1 + df.m$conc)  # change conc to log-scale
      df.m$conc[df.m$Arm=="C"] = df.m$conc.log[df.m$Arm=="C"] = 0  # assign conc of placebo recipients 0 
      
      #determine AUC up until interval
      df.m.auc <- df.m %>%
        arrange(id, tstart) %>%
        group_by(id) %>%
        mutate(AUC_interval = ifelse(is.na(lag(conc)) & lag(tstart) == 0,
                                     ((conc+0)/2)*(tstop-tstart),
                                     ((conc+lag(conc))/2)*(tstop-tstart))) %>%
        mutate(AUC_cumsum = cumsum(ifelse(is.na(AUC_interval), 0, AUC_interval)) + AUC_interval*0) %>%
        mutate(AUC_cumsum_event = last(AUC_cumsum)) %>% #in Zonghui's dadta, the data ends at the last event or censor time, so you can use the last cumsum AUC as your AUC variable
        ungroup()
      # association analysis based on time to first infection 
      # Cox regression with time-varying covariate is performed with randomization arm (mAb or placebo, time-fixed) and mAb concentration (10 days before each pf evaluation, time-varying) as regressors.  
      
      # Model 1: h(t) = h0(t) exp{ Z [ B0 + B1 mAb(t)]}, VE(mAb) =  1 – exp(B0+B1 mAb) 
      # To be consistent with the VRC 612 Part C study (Lyke et al Lancet Infect Dis 2023), we will use 90% protection here.
      # Cox regression 1:
      df.m$arm = 1
      df.m$arm[df.m$Arm=="C"] = 0
      
      # conc in log scale fits better 
      # df.m <- df.m %>%
      #  mutate(arm = factor(ifelse(arm == 1, "CIS43LS", "placebo"), levels = c("placebo", "CIS43LS")))
      fit <- coxph(Surv(tstart,tstop,endpt)~ arm + arm:conc.log + cluster(id),data = df.m)
      summary.fit <- summary(fit)
      cox.fit.log = summary.fit$coefficients[,1]
      conc.plot = range(df.m$conc[df.m$conc>0], na.rm=TRUE)
      conc.plot = seq(conc.plot[1], conc.plot[2], length=3000)
      # conc corresponding to XX% efficacy
      x <- exp((cox.fit.log[1] - log(1-efficacy_threshold))/(-cox.fit.log[2])) -1
      PE = 1- exp(cox.fit.log[1] + cox.fit.log[2]*log(1+conc.plot))
      # plot(conc.plot, PE, xlab="CIS43LS concentration", ylab="Protective efficacy", type="l", ylim=c(0, 1))
      # abline(h=efficacy_threshold, lty=3, col="red")
      # abline(v=x, lty=3, col="blue")
      
      # find confidence band
      ID = unique(df.m$id)
      VEC <- matrix(rep(0,nboot*length(conc.plot)),ncol=length(conc.plot))
      for(ijk in 1:nboot){
        samp = sample(ID, replace=TRUE)
        dat = df.m[df.m$id %in%samp,]
        coef <- summary(coxph(Surv(tstart,tstop,endpt)~ arm + conc.log+cluster(id),data = dat))$coef[,1]
        VEC[ijk,] = 1- exp(coef[1] + coef[2]*log(1+conc.plot))
        }
      LOWER = apply(VEC, 2, quantile, p=0.025)
      UPPER = apply(VEC, 2, quantile, p=0.975)
      # lines(conc.plot,LOWER,lty=2,lwd=2, col="green")
      # lines(conc.plot,UPPER,lty=2,lwd=2, col="green")
      
      ggplot_df <- data.frame(conc = conc.plot, pe = PE, lower_ci = LOWER, upper_ci = UPPER)
      mybreaks <- seq(28,168, by=28)
      max_limit <- 1000
      PE_conc_forward <- ggplot_df %>%
        filter(conc <= max_limit) %>%
        mutate(lower_ci = ifelse(lower_ci<0,0, lower_ci)) %>%
        mutate(lower_ci = ifelse(pe<0,0, lower_ci)) %>%
        drop_na(conc, pe) %>%
        ggplot(., aes(x = conc, y = pe)) +
        geom_ribbon(aes(ymin=lower_ci, ymax=upper_ci), linetype=2, alpha=0.3) +
        geom_line() +
        ylab("protective efficacy") +
        xlab("CIS43LS concentration (µg/mL)") +
        scale_x_continuous(expand=c(0, 0), breaks = c(0,as.numeric(x),seq(200,max_limit, by = 200)), labels = c(0,"",seq(200,max_limit, by = 200)), limits = c(-500,max_limit)) +
        scale_y_continuous(expand=c(0, 0), breaks = c(0,0.25,0.50,0.75,efficacy_threshold,1), limits = c(-0.05,1)) +
        geom_vline(xintercept = x, linetype = "dotted", color = "red") +
        geom_hline(yintercept = efficacy_threshold, linetype = "dotted", color = "blue") +
        annotate(geom="text", x=x+60, y=0.03, label=signif(as.numeric(x),sigfigs), size = 4,color="red") +
        theme_bw() +
        theme(legend.title=element_blank(),
              axis.text = element_text(size = 9),
              axis.title = element_text(size=11),
              strip.background = element_blank(),
              strip.text = element_text(size = 12),
              legend.text = element_text(size=12),
              legend.position = "none",
              plot.margin = margin(0.5,0.5,0.5,0.5, 'cm')) +
        coord_cartesian(xlim = c(-50,max_limit+10))
      results <- list("long_data" = PK.matched,
                      "event_data" = df.m,
                      "auc2event_data" = df.m.auc,
                      "cox_res" = summary.fit,
                      "conc_at_PE" = x,
                      "plot" = PE_conc_forward)
      }
}