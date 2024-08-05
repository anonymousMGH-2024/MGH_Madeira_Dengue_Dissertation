#### Before running this code please run 00_run.R lines 2-11, 16-25, 44-72
#### This code produces the model fitting trace plots and histograms, as well as plots of parameters


library(ggpubr)
library(cowplot)
require(tidyverse)
require(ggplot2)

##these have to be the same as in the chains named below in the filenames
date_start<- as.Date("2012-01-01","%Y-%m-%d") -1
time_start <- 0
time_stop <- 485
delta_T<- 0.0625 
tps <- seq(time_start , time_stop , by = delta_T)

source("parameter_functions.R")
source("MCMC_functions.R")
source("set_up_data_for_modelling.R")
source('MODEL_C_odeintr.R')  ##this sets up the model in C

# ############################
# ############################
# ############################

############################

NAME_RUN<- "2millRESULT"
BURNIN_PROPORTION<- 0.65
BURNED_SAMPLES_N<- 200
BURNED_SAMPLES_SHUFFLE<- FALSE
chain_files<- c(
        "Outputs/2millRESULT.csv_fullchains.Rdata")

############################
############################
############################
###########################
## join all MCMC chains, export chains and posteriors

        entire_chain<- c()
        for(chain in chain_files){
            load(chain) ##saved_chains  
            if("like" %in% colnames(saved_chains)) saved_chains<- saved_chains %>% select(-"like")
            entire_chain<- rbind(entire_chain, saved_chains)
        }
        head(entire_chain)
        dim(entire_chain)

        chain_data<- entire_chain
        chain_data$step<- 1:nrow(chain_data)

        ######################### chains

        gg_chains_k<- ggplot(chain_data) + theme_classic() +
                geom_line(aes(x=step, y=K)) +
                geom_vline(xintercept=max(chain_data$step)*BURNIN_PROPORTION, color="red") 

        gg_chains_psi<- ggplot(chain_data) + theme_classic() +
                geom_line(aes(x=step, y=psi)) +
                geom_vline(xintercept=max(chain_data$step)*BURNIN_PROPORTION, color="red") 

        gg_chains_alpha<- ggplot(chain_data) + theme_classic() +
                geom_line(aes(x=step, y=alpha)) +
                geom_vline(xintercept=max(chain_data$step)*BURNIN_PROPORTION, color="red") 

        gg_chains_eta<- ggplot(chain_data) + theme_classic() +
                geom_line(aes(x=step, y=eta)) +
                geom_vline(xintercept=max(chain_data$step)*BURNIN_PROPORTION, color="red") 

        # gg_chains_delta<- ggplot(chain_data) + theme_classic() +
        #         geom_line(aes(x=step, y=delta)) +
        #         geom_vline(xintercept=max(chain_data$step)*BURNIN_PROPORTION, color="red") 

        gg_chains_t<- ggplot(chain_data) + theme_classic() +
                geom_line(aes(x=step, y=T)) +
                geom_vline(xintercept=max(chain_data$step)*BURNIN_PROPORTION, color="red") 

        # png(paste0("Outputs/",NAME_RUN,"_all_chain.png"),w=1000,h=600, res=100)
        # print(cowplot::plot_grid(gg_chains_k,
        #                             gg_chains_psi,
        #                             gg_chains_alpha,
        #                             gg_chains_eta,
        #                             #gg_chains_delta,
        #                             gg_chains_t, nrow=2))
        # dev.off()

        # data_sub <- subset(chain_data, chain_data$step > max(chain_data$step*BURNIN_PROPORTION))
        
######################### posteriors

        posterior_chains<- saved_chains[(nrow(saved_chains)*BURNIN_PROPORTION):nrow(saved_chains),]
        posterior_chains$step<- 1:nrow(posterior_chains)

        xi<- 0.8; xf<- 1.2
        bins<- 50

        xran<- c(min(posterior_chains$K)*xi,max(posterior_chains$K)*xf)
        gg_post_k<- ggplot(posterior_chains) + theme_classic() +
                    geom_histogram(aes(x=K),bins=bins) + scale_x_continuous(lim=xran) 

        xran<- c(min(posterior_chains$psi)*xi,max(posterior_chains$psi)*xf)
        gg_post_psi<- ggplot(posterior_chains) + theme_classic() +
                    geom_histogram(aes(x=psi),bins=bins) + scale_x_continuous(lim=xran) 

        xran<- c(min(posterior_chains$alpha)*xi,max(posterior_chains$alpha)*xf)
        gg_post_alpha<- ggplot(posterior_chains) + theme_classic() +
                    geom_histogram(aes(x=alpha),bins=bins) + scale_x_continuous(lim=xran) 

        xran<- c(min(posterior_chains$eta)*xi,max(posterior_chains$eta)*xf)
        gg_post_eta<- ggplot(posterior_chains) + theme_classic() +
                    geom_histogram(aes(x=eta),bins=bins) + scale_x_continuous(lim=xran) 

        # xran<- c(min(posterior_chains$delta)*xf,max(posterior_chains$delta)*xi)
        # gg_post_delta<- ggplot(posterior_chains) + theme_classic() +
        #             geom_histogram(aes(x=delta),bins=bins) + scale_x_continuous(lim=xran) 

        xran<- c(min(posterior_chains$T)*xi,max(posterior_chains$T)*xf)
        gg_post_t<- ggplot(posterior_chains) + theme_classic() +
                    geom_histogram(aes(x=T),bins=bins) + scale_x_continuous(lim=xran) 

        # png(paste0("Outputs/",NAME_RUN,"_post.png"),w=1000,h=600, res=100)
        # print(cowplot::plot_grid(gg_post_k,
        #                             gg_post_psi,
        #                             gg_post_alpha,
        #                             gg_post_eta,
        #                             #gg_post_delta,
        #                             gg_post_t, nrow=2))
        # dev.off()

############################
###########################
## use burned posteriors to make sims and get final dynamical results

    samples<- sample(1:nrow(posterior_chains),BURNED_SAMPLES_N)
    samples_eta<- posterior_chains$eta[samples]
    samples_alpha<- posterior_chains$alpha[samples]
    samples_a<- posterior_chains$a[samples]
    samples_sigH<- posterior_chains$sigH[samples]
    samples_K<- posterior_chains$K[samples]
    samples_gamH<- posterior_chains$gamH[samples]
    samples_T<- posterior_chains$T[samples]
    samples_i<- posterior_chains$i[samples]
    samples_psi<- posterior_chains$psi[samples]
    samples_delta<- posterior_chains$delta[samples]

    if(BURNED_SAMPLES_SHUFFLE){
        samples_eta<- sample(samples_eta)
        samples_alpha<- sample(samples_alpha)
        samples_a<- sample(samples_a)
        samples_sigH<- sample(samples_sigH)
        samples_K<- sample(samples_K)
        samples_gamH<- sample(samples_gamH)
        samples_T<- sample(samples_T)
        samples_i<- sample(samples_i)
        samples_psi<- sample(samples_psi)
        samples_delta<- sample(samples_delta)
    }
    ##########################################################################
    ##########################################################################

    ncol_run_model<- length(tps)
    ncol_run_pars<- length(tps)
    nrow_run<-  BURNED_SAMPLES_N

    sample_phiVH<- data.frame(matrix(rep(0,ncol_run_pars*nrow_run), ncol=ncol_run_pars,nrow=nrow_run))
    sample_phiHV<- data.frame(matrix(rep(0,ncol_run_pars*nrow_run), ncol=ncol_run_pars,nrow=nrow_run))
    sample_cV<- data.frame(matrix(rep(0,ncol_run_pars*nrow_run), ncol=ncol_run_pars,nrow=nrow_run))
    sample_aV<- data.frame(matrix(rep(0,ncol_run_pars*nrow_run), ncol=ncol_run_pars,nrow=nrow_run))
    sample_thetaV<- data.frame(matrix(rep(0,ncol_run_pars*nrow_run), ncol=ncol_run_pars,nrow=nrow_run))
    sample_epsA<- data.frame(matrix(rep(0,ncol_run_pars*nrow_run), ncol=ncol_run_pars,nrow=nrow_run))
    sample_muA<- data.frame(matrix(rep(0,ncol_run_pars*nrow_run), ncol=ncol_run_pars,nrow=nrow_run))
    sample_muV<- data.frame(matrix(rep(0,ncol_run_pars*nrow_run), ncol=ncol_run_pars,nrow=nrow_run))
    sample_gamV<- data.frame(matrix(rep(0,ncol_run_pars*nrow_run), ncol=ncol_run_pars,nrow=nrow_run))
    sample_gamV<- data.frame(matrix(rep(0,ncol_run_pars*nrow_run), ncol=ncol_run_pars,nrow=nrow_run))

    sample_R0<- data.frame(matrix(rep(0,ncol_run_model*nrow_run), ncol=ncol_run_model,nrow=nrow_run))
    sample_A<- data.frame(matrix(rep(0,ncol_run_model*nrow_run), ncol=ncol_run_model,nrow=nrow_run))
    sample_V<- data.frame(matrix(rep(0,ncol_run_model*nrow_run), ncol=ncol_run_model,nrow=nrow_run))
    sample_INC_obs<- data.frame(matrix(rep(0,ncol_run_model*nrow_run), ncol=ncol_run_model,nrow=nrow_run))
    sample_CSINC_obs<- data.frame(matrix(rep(0,ncol_run_model*nrow_run), ncol=ncol_run_model,nrow=nrow_run))
    sample_model_time<- data.frame(matrix(rep(0,ncol_run_model*nrow_run), ncol=ncol_run_model,nrow=nrow_run))
    
    sample_CSINC<- data.frame(matrix(rep(0,ncol_run_model*nrow_run), ncol=ncol_run_model,nrow=nrow_run))

    for(rr in 1:BURNED_SAMPLES_N){
        print(paste("running",rr,"of",BURNED_SAMPLES_N,"..."))

        p_eta<- samples_eta[rr]
        p_alpha<- samples_alpha[rr]
        p_a<- samples_a[rr]
        p_sigH<- samples_sigH[rr]
        p_K<- samples_K[rr]
        p_gamH<- samples_gamH[rr]
        p_T<- samples_T[rr]
        p_i<- samples_i[rr]
        p_psi<- samples_psi[rr]
        p_delta<- samples_delta[rr]

        ##this is how the parameters are calculated in the C model

        phiVH<- model_input$phi_VH
        phiHV<- model_input$phi_HV
        cV<- model_input$c_temp
        thetaV<- model_input$theta
        R<- model_input$prec
        epsA<- model_input$epsilon
        muA<- (p_eta * model_input$mu_A) * (1 + model_input$mu_A_rain)^p_eta
        muV<- (p_eta * model_input$mu_V_temp) * (1 + model_input$mu_V_hum)^p_eta
        gamV<- (p_alpha * model_input$gamma_V) #* (1 + model_input$gam_V_temp)^p_alpha
        aV<- (p_a) #* (1 + model_input$a_hum)^p_eta
        muH<- 0

        ###same checks as in the C model code

        NEAR_ZERO= 0.00000000001
        NEAR_ONE= 0.999999999
        MIN_MUV= 1/60 ##assume max life exp is 60 days
        MIN_MUA= 1/120 ##assume max life exp is 120 days
        phiVH[phiVH<0.0] <- NEAR_ZERO
        phiVH[phiVH>1.0] <- NEAR_ONE
        phiHV[phiHV<0.0] <- NEAR_ZERO
        phiHV[phiHV>1.0] <- NEAR_ONE
        cV[cV<0.0] <- NEAR_ZERO
        cV[cV>1.0] <- NEAR_ONE
        epsA[epsA<0.0] <- NEAR_ZERO
        muA[muA<MIN_MUA] <- MIN_MUA
        # muA[muA>1.0] <- NEAR_ONE
        muV[muV<MIN_MUV] <- MIN_MUV
        # muV[muV>1.0] <- NEAR_ONE
        gamV[gamV<0.0] <- NEAR_ZERO
        aV[aV<0.0] <- NEAR_ZERO   
        aV[aV>1.0] <- NEAR_ONE  

        ###save

        set(sample_phiVH, i=rr, j=1:ncol(sample_phiVH), value=as.list(as.numeric(phiVH)))
        set(sample_phiHV, i=rr, j=1:ncol(sample_phiHV), value=as.list(as.numeric(phiHV)))
        set(sample_cV, i=rr, j=1:ncol(sample_cV), value=as.list(as.numeric(cV)))
        set(sample_aV, i=rr, j=1:ncol(sample_aV), value=as.list(as.numeric(aV)))
        set(sample_thetaV, i=rr, j=1:ncol(sample_thetaV), value=as.list(as.numeric(thetaV)))
        set(sample_epsA, i=rr, j=1:ncol(sample_epsA), value=as.list(as.numeric(epsA)))
        set(sample_muA, i=rr, j=1:ncol(sample_muA), value=as.list(as.numeric(muA)))
        set(sample_muV, i=rr, j=1:ncol(sample_muV), value=as.list(as.numeric(muV)))
        set(sample_gamV, i=rr, j=1:ncol(sample_gamV), value=as.list(as.numeric(gamV)))
        

        ##simulate

            model_init_bf_intro<- c(NH,0,0,0, p_K,0,0,0, 0,0, 0,0,0) ##initial conditions before intro
            ODEmodel_set_params(ts=delta_T, alpha=p_alpha, eta=p_eta, gamH=p_gamH, sigH=p_sigH, a=p_a, f=prior_f, K=p_K, delta=p_delta)
            model_res_bf_intro<- ODEmodel(init=model_init_bf_intro, start=time_start, duration=model_input$timestep[p_T], step_size=delta_T)
            colnames(model_res_bf_intro)<- c('time', "SH","EH","IH","RH",  "AV","SV","EV","IV",  "Inc","IncCS", "R0","Vls","EIC")
            last_state<- model_res_bf_intro[nrow(model_res_bf_intro),]
            
            model_init_af_intro<- c(NH-p_i,0,p_i,0, last_state$AV,last_state$SV-p_i,last_state$EV+p_i,last_state$IV, p_i,p_i, last_state$R0, last_state$Vls, last_state$EIC) ##initial conditions a_T intro
            ODEmodel_set_params(ts=delta_T, alpha=p_alpha, eta=p_eta, gamH=p_gamH, sigH=p_sigH, a=p_a, f=prior_f, K=p_K, delta=p_delta)
            model_res_af_intro<- ODEmodel(init=model_init_af_intro, start=last_state$time+delta_T, duration=time_stop-(last_state$time+delta_T), step_size=delta_T)
            colnames(model_res_af_intro)<- c('time', "SH","EH","IH","RH",  "AV","SV","EV","IV",  "Inc","IncCS", "R0","Vls","EIC")
            
            model_res<- rbind(model_res_bf_intro, model_res_af_intro)

            R0<- model_res$R0
            V<- model_res$SV+model_res$EV+model_res$IV
            A<- model_res$AV
            obs_inc<- model_res$Inc * p_psi
            obs_cs_inc<- model_res$IncCS * p_psi
            cs_inc<- model_res$IncCS

            set(sample_R0, i=rr, j=1:ncol(sample_R0), value=as.list(R0))
            set(sample_V, i=rr, j=1:ncol(sample_V), value=as.list(V))
            set(sample_A, i=rr, j=1:ncol(sample_A), value=as.list(A))
            set(sample_INC_obs, i=rr, j=1:ncol(sample_INC_obs), value=as.list((obs_inc)))
            set(sample_CSINC_obs, i=rr, j=1:ncol(sample_CSINC_obs), value=as.list((obs_cs_inc)))
            set(sample_CSINC, i=rr, j=1:ncol(sample_CSINC), value=as.list((cs_inc)))
    }

    model_time<- model_res$time ##time should be same across models

    #MATCH ESTIMATED CASES WITH REPORTED CASES IN THE RIGHT TIME SCALE
    timeframe<- data.frame(time=model_time, timestep=round(model_time,1))
    timeframe_matches<- match(timeframe$timestep,reportedcases$timestep)
    FITTING_MATCHES<- which(!is.na(timeframe_matches))
    MODEL_REPORTED_TIMEMATCHES<- timeframe$timestep[FITTING_MATCHES]
    MODEL_REPORTED_TIMEMATCHES<- unique(MODEL_REPORTED_TIMEMATCHES)

    sample_match_CSINC_obs<- sample_CSINC_obs[,c(FITTING_MATCHES[1]-diff(FITTING_MATCHES)[1],FITTING_MATCHES)]
    sample_match_INC_obs<- t(apply(sample_match_CSINC_obs, MARG=1, FUN=function(x){ x<-as.numeric(x); diff(x) }))
    sample_match_CSINC_obs<- sample_match_CSINC_obs[,-1]
    
    sample_match_CSINC<- sample_CSINC[,c(FITTING_MATCHES[1]-diff(FITTING_MATCHES)[1],FITTING_MATCHES)]
    sample_match_INC<- t(apply(sample_match_CSINC, MARG=1, FUN=function(x){ x<-as.numeric(x); diff(x) }))
    sample_match_CSINC<- sample_match_CSINC[,-1]
    
    sample_match_INC_obs_mean<- as.numeric(colMeans(sample_match_INC_obs))
    sample_match_INC_obs_sd<- as.numeric(apply(sample_match_INC_obs, MARG=2, FUN=sd))
    sample_match_INC_obs_upper<- as.numeric(apply(sample_match_INC_obs, MARG=2, FUN=quantile, probs=0.975))
    sample_match_INC_obs_lower<- as.numeric(apply(sample_match_INC_obs, MARG=2, FUN=quantile, probs=0.025))

    sample_match_INC_mean<- as.numeric(colMeans(sample_match_INC))
    sample_match_INC_sd<- as.numeric(apply(sample_match_INC, MARG=2, FUN=sd))
    sample_match_INC_upper<- as.numeric(apply(sample_match_INC, MARG=2, FUN=quantile, probs=0.975))
    sample_match_INC_lower<- as.numeric(apply(sample_match_INC, MARG=2, FUN=quantile, probs=0.025))

    sample_match_CSINC_obs_mean<- as.numeric(colMeans(sample_match_CSINC_obs))
    sample_match_CSINC_obs_sd<- as.numeric(apply(sample_match_CSINC_obs, MARG=2, FUN=sd))
    sample_match_CSINC_obs_upper<- as.numeric(apply(sample_match_CSINC_obs, MARG=2, FUN=quantile, probs=0.975))
    sample_match_CSINC_obs_lower<- as.numeric(apply(sample_match_CSINC_obs, MARG=2, FUN=quantile, probs=0.025))

    sample_match_CSINC_mean<- as.numeric(colMeans(sample_match_CSINC))
    sample_match_CSINC_sd<- as.numeric(apply(sample_match_CSINC, MARG=2, FUN=sd))
    sample_match_CSINC_upper<- as.numeric(apply(sample_match_CSINC, MARG=2, FUN=quantile, probs=0.975))
    sample_match_CSINC_lower<- as.numeric(apply(sample_match_CSINC, MARG=2, FUN=quantile, probs=0.025))

    
    sample_match<- data.frame(timestep=MODEL_REPORTED_TIMEMATCHES,
                              model_INC_obs_mean=sample_match_INC_obs_mean,
                              model_INC_obs_sd=sample_match_INC_obs_sd,
                              model_INC_obs_upper=sample_match_INC_obs_upper,
                              model_INC_obs_lower=sample_match_INC_obs_lower,
                              model_CSINC_obs_mean=sample_match_CSINC_obs_mean,
                              model_CSINC_obs_sd=sample_match_CSINC_obs_sd,
                              model_CSINC_obs_upper=sample_match_CSINC_obs_upper,
                              model_CSINC_obs_lower=sample_match_CSINC_obs_lower,
                              model_INC_mean=sample_match_INC_mean,
                              model_INC_sd=sample_match_INC_sd,
                              model_INC_upper=sample_match_INC_upper,
                              model_INC_lower=sample_match_INC_lower,
                              model_CSINC_mean=sample_match_CSINC_mean,
                              model_CSINC_sd=sample_match_CSINC_sd,
                              model_CSINC_upper=sample_match_CSINC_upper,
                              model_CSINC_lower=sample_match_CSINC_lower)

    fitting<- reportedcases %>% dplyr::left_join(sample_match, by="timestep")

    seropositivity<- (sample_CSINC_obs[,ncol(sample_CSINC_obs)]*1/samples_psi)/NH

    inc_mean<- as.numeric(apply(sample_INC_obs, MARG=2, FUN=mean))
    inc_sd<- as.numeric(apply(sample_INC_obs, MARG=2, FUN=sd))

    phiVH_mean<- as.numeric(apply(sample_phiVH, MARG=2, FUN=mean))
    phiVH_sd<- as.numeric(apply(sample_phiVH, MARG=2, FUN=sd))

    phiHV_mean<- as.numeric(apply(sample_phiHV, MARG=2, FUN=mean))
    phiHV_sd<- as.numeric(apply(sample_phiHV, MARG=2, FUN=sd))

    V_mean<- as.numeric(apply(sample_V, MARG=2, FUN=mean))
    V_sd<- as.numeric(apply(sample_V, MARG=2, FUN=sd))

    V_all<- as.numeric(unlist(sample_V))
    VNH_all<- V_all/NH
    print(paste("mean V",mean(V_all), "SD V",sd(V_all)))
    print(paste("mean V/NH",mean(VNH_all), "SD V",sd(VNH_all)))

    A_mean<- as.numeric(apply(sample_A, MARG=2, FUN=mean))
    A_sd<- as.numeric(apply(sample_A, MARG=2, FUN=sd))

    R0_mean<- as.numeric(apply(sample_R0, MARG=2, FUN=mean))
    R0_sd<- as.numeric(apply(sample_R0, MARG=2, FUN=sd))

    cV_mean<- as.numeric(apply(sample_cV, MARG=2, FUN=mean))
    cV_sd<- as.numeric(apply(sample_cV, MARG=2, FUN=sd))

    aV_mean<- as.numeric(apply(sample_aV, MARG=2, FUN=mean))
    aV_sd<- as.numeric(apply(sample_aV, MARG=2, FUN=sd))

    thetaV_mean<- as.numeric(apply(sample_thetaV, MARG=2, FUN=mean))
    thetaV_sd<- as.numeric(apply(sample_thetaV, MARG=2, FUN=sd))

    epsA_mean<- as.numeric(apply(sample_epsA, MARG=2, FUN=mean))
    epsA_sd<- as.numeric(apply(sample_epsA, MARG=2, FUN=sd))

    muA_mean<- as.numeric(apply(sample_muA, MARG=2, FUN=mean))
    muA_sd<- as.numeric(apply(sample_muA, MARG=2, FUN=sd))

    muV_mean<- as.numeric(apply(sample_muV, MARG=2, FUN=mean))
    muV_sd<- as.numeric(apply(sample_muV, MARG=2, FUN=sd))

    gamV_mean<- as.numeric(apply(sample_gamV, MARG=2, FUN=mean))
    gamV_sd<- as.numeric(apply(sample_gamV, MARG=2, FUN=sd))

    Q_mean<- (epsA_mean/(epsA_mean+muA_mean)) * cV_mean*prior_f*thetaV_mean/muV_mean

    




    #pdf(paste0("Outputs/",NAME_RUN,"_mechanistic_summary.pdf"), w=7, h=5)

    
        ##### plot of number of cases - with confidence interval first and model cases first
    
        mdates<- date_start + fitting$timestep
        plot(mdates, fitting$model_INC_upper*p_psi, t='l', col="transparent", xlab = "Month (2012-2013)", ylab = "Number of Cases")
        lines(mdates, fitting$model_INC_mean*p_psi, t='b', col="tomato")
        points(mdates, fitting$Cases, t='b', col="black")
        #lines(mdates, fitting$model_INC_mean*p_psi-fitting$model_INC_obs_sd*p_psi, t='l', col="tomato")
        #lines(mdates, fitting$model_INC_mean*p_psi+fitting$model_INC_obs_sd*p_psi, t='l', col="tomato")
        lines(mdates, fitting$model_INC_lower*p_psi, t='l', col="transparent")
        abline(v=as.numeric(date_start +mean(samples_T)/(1/delta_T)))
        polygon(c(mdates, rev(mdates)), c(fitting$model_INC_upper*p_psi, rev(fitting$model_INC_lower*p_psi)),
                col=rgb(1, 0, 0,0.1), border = NA)

      
        #### cases without reporting rate subtracted
        
        
        mdates<- date_start + fitting$timestep
        plot(mdates, fitting$model_INC_upper, t='l', col="transparent", xlab = "Month (2012-2013)", ylab = "Number of Cases")
        lines(mdates, fitting$model_INC_mean, t='b', col="tomato")
        points(mdates, fitting$Cases, t='b', col="black")
        #lines(mdates, fitting$model_INC_mean-fitting$model_INC_obs_sd, t='l', col="tomato")
        #lines(mdates, fitting$model_INC_mean+fitting$model_INC_obs_sd, t='l', col="tomato")
        lines(mdates, fitting$model_INC_lower, t='l', col="transparent")
        abline(v=date_start +mean(samples_T)/(1/delta_T))
        polygon(c(mdates, rev(mdates)), c(fitting$model_INC_upper, rev(fitting$model_INC_lower)),
                col=rgb(1, 0, 0,0.1), border = NA)


        ##### plot of cumulative incidence cases - with model cases first and confidence interval first
        
        mdates<- date_start + fitting$timestep
        plot(mdates, fitting$model_CSINC_upper*p_psi, t='b', col="transparent", xlab = "Month (2012-2013)", ylab = "Cumulative Incidence")
        points(mdates, fitting$CSCases, t='b', col="black")
        #lines(mdates, fitting$model_CSINC_mean*p_psi-fitting$model_CSINC_obs_sd*p_psi, t='l', col="tomato")
        #lines(mdates, fitting$model_CSINC_mean*p_psi+fitting$model_CSINC_obs_sd*p_psi, t='l', col="tomato")
        lines(mdates, fitting$model_CSINC_mean*p_psi, t='b', col="tomato")
        lines(mdates, fitting$model_CSINC_lower*p_psi, t='l', col="transparent")
        abline(v=date_start +mean(samples_T)/(1/delta_T))
        polygon(c(mdates, rev(mdates)), c(fitting$model_CSINC_upper*p_psi, rev(fitting$model_CSINC_lower*p_psi)),
                col=rgb(1, 0, 0,0.1), border = NA)
        
        
        #### R0 plot
         
        mdates<- date_start + model_time
        plot(mdates, R0_mean, t='l', col="black", xlab = "Date", ylab = "Reproductive number (R0)")
        lines(mdates, R0_mean-R0_sd, t='l', col="darkgrey")
        lines(mdates, R0_mean+R0_sd, t='l', col="darkgrey")
        abline(h=1,col="#4d4deb", lty = 2)
        lines(mdates, R0_mean, t='l', col="black", xlab = "Date", ylab = "Reproductive number (R0)")
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")
        
        # R0<- unlist(sample_R0)
        # print(paste("R0 mean",mean(R0),"median",median(R0),"sd",sd(R0)))
        
        
        #### R0 and temperature plot
        
        par(mar=c(5, 4, 4, 5) + 0.1)
        plot(mdates, R0_mean, t='l', col="black", xlab = "Date", ylab = "Reproductive number (R0)")
        abline(h=1,col="#4d4deb", lty = 2)
        par(new=TRUE)
        plot(mdates, new_temperature, col="#971fcc", t='l', ylab = "", xlab = "", axes = FALSE)
        axis(side=4)
        mtext("Temperature (Celcius)", side=4, line=3)
        legend("topleft", inset=.04,legend=c("R0", "Temperature", "R0 of 1", "Introduction date"),
               col=c("black", "#971fcc", "#4d4deb", "tomato"), lty=c(1,1,2), cex=0.7)
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")
        
        #### Humidity and precipitation plot
        
        par(mar=c(5, 4, 4, 5) + 0.1)
        plot(mdates, new_humidity*100, t='l', col="#f7ae60", xlab = "Date", ylab = "Humidity (%)")
        abline(h=1,col="#4d4deb")
        par(new=TRUE)
        plot(mdates,new_precipitation, col="#f760c5", t='l', ylab = "", xlab = "", axes = FALSE)
        axis(side=4)
        mtext("Precipitation (mm)", side=4, line=3)
        legend("topleft", inset=.04,legend=c("Humidity", "Precipitation", "Introduction date"),
               col=c("#f7ae60", "#f760c5", "tomato"), lty=1, cex=0.7)
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")
        
        
        #### Histogram of seroprevelance (count and density)
        
        {phist <- gghistogram(
          as.data.frame(seropositivity), x = "seropositivity", fill="#fa8e8e", ylab = "Count") +
          theme_half_open(11, rel_small = 1) +
          rremove("x.axis")+
          rremove("xlab") +
          rremove("x.text") +
          rremove("x.ticks") +
          rremove("legend")
        
        pdensity <- ggdensity(
          as.data.frame(seropositivity), x = "seropositivity", fill="#fcf0f0", color = "#6e6d6d", title = "Histogram of Seroprevelance", xlab = "Seropositivity", ylab = "Density") +
          scale_y_continuous(position = "right") +
          theme(plot.title = element_text(hjust = 0.5))
        
        aligned_plots <- align_plots(pdensity, phist, align="hv", axis="tblr")
        ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])}
        
        
        ####  phi plots
        mdates<- date_start + tps
        plot(mdates, phiVH_mean, t='l', col="blue", xlab = "Date", ylab = "Probability of transmission from vector to human")
        # lines(mdates, phiVH_mean-phiVH_sd, t='l', col="grey33")
        # lines(mdates, phiVH_mean+phiVH_sd, t='l', col="grey33")
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")

        mdates<- date_start + tps
        plot(mdates, phiHV_mean, t='l', col="blue", xlab = "Date", ylab = "Probability of transmission from human to vector")
        # lines(mdates, phiHV_mean-phiHV_sd, t='l', col="grey33")
        # lines(mdates, phiHV_mean+phiHV_sd, t='l', col="grey33")
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")

        #### mosquito population plots
        mdates<- date_start + model_time
        plot(mdates, V_mean/NH, t='l', col="blue", xlab = "Date", ylab = "Adult mosquitos per human")
        # lines(mdates, (V_mean-V_sd)/NH, t='l', col="transparent")
        # lines(mdates, (V_mean+V_sd)/NH, t='l', col="transparent")
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")
        polygon(c(mdates, rev(mdates)), c((V_mean-V_sd)/NH, rev((V_mean+V_sd)/NH)),
                col=rgb(0, 0, 1,0.08), border = NA)
        

        mdates<- date_start + model_time
        plot(mdates, A_mean/NH, t='l', col="blue", xlab = "Date", ylab = "Larval mosquitos per human")
        lines(mdates, (A_mean-A_sd)/NH, t='l', col="transparent")
        lines(mdates, (A_mean+A_sd)/NH, t='l', col="transparent")
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")
        polygon(c(mdates, rev(mdates)), c((A_mean-A_sd)/NH, rev((A_mean+A_sd)/NH)),
                col=rgb(0, 0, 1,0.08), border = NA)

        mdates<- date_start + model_time
        plot(mdates, A_mean, t='l', col="blue", xlab = "Date", ylab = "Number of larval mosquitos")
        lines(mdates, (A_mean-A_sd), t='l', col="transparent")
        lines(mdates, (A_mean+A_sd), t='l', col="transparent")
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")
        polygon(c(mdates, rev(mdates)), c((A_mean-A_sd), rev((A_mean+A_sd))),
                col=rgb(0, 0, 1,0.08), border = NA)
        
        
        #### c plot
        mdates<- date_start + tps
        plot(mdates, cV_mean, t='l', col="blue", xlab = "Date", ylab = "Fraction of eggs hatching")
        # lines(mdates, cV_mean-cV_sd, t='l', col="grey33")
        # lines(mdates, cV_mean+cV_sd, t='l', col="grey33")
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")

        
        #### Q plot
        mdates<- date_start + tps
        plot(mdates, Q_mean, t='l', col="blue", xlab = "Date", ylab = "Basic offspring number (Q)")
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")

        #### Theta plot
        mdates<- date_start + tps
        plot(mdates, thetaV_mean, t='l', col="blue", xlab = "Date", ylab = "Oviposition rate")
        # lines(mdates, thetaV_mean-thetaV_sd, t='l', col="grey33")
        # lines(mdates, thetaV_mean+thetaV_sd, t='l', col="grey33")
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")


        
        #### muA and eps together plot
        
        mdates<- date_start + tps
        plot(mdates, 1/(epsA_mean), t='l', col="blue", xlab = "Date", ylab = "Days")
        lines(mdates, 1/(muA_mean), t='l', col="orange")
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")
        legend("topright", inset=.04,legend=c("epsilon", "muA", "Introduction date"),
               col=c("blue", "orange", "tomato"), lty=1:1, cex=0.7)
        
        
        
        #### muV and gamV together plot
        
        mdates<- date_start + tps
        plot(mdates, 1/(muV_mean), t='l', col="blue", ylim=c(0,25), xlab = "Date", ylab = "Days")
        lines(mdates, 1/(gamV_mean), t='l', col="orange")
        abline(v=date_start +mean(samples_T)/(1/delta_T), col = "tomato")
        legend("topright", inset=.04,legend=c("muV", "gammaV", "Introduction date"),
               col=c("blue", "orange", "tomato"), lty=1:1, cex=0.7)
        
        

    #dev.off()
    
    
    
