#### Please do not run this file directly as it will not work.
#### Please refer to file 00_run.R

if(RUN_MCMC){


      ##liu et al 2020 for temperate regions overall
        prior_R0_estimates<- c(1, 1.8, 7.78, 0.82, 0.36, 0.4, 2.07, 0.68, 0.46, 0.48, 3.9, 2.035, 1.96, 2.36, 2.94)
        prior_R0_me<- mean(prior_R0_estimates)
        prior_R0_sd<- sd(prior_R0_estimates)

      ###Nakase 2023
      prior_a_me<- 0.25

      ##taken from the paper Auerswald et al 2019
      seropositive_i<- 0.078 ##proportion of popula_Tion
      seropositive_ii<- 0.089 ##proportion of popula_Tion
      prior_seropositive_me<- mean(c(seropositive_i,seropositive_ii))
      prior_seropositive_sd<- sd(c(seropositive_i,seropositive_ii))
      

      ## taken from taishi paper
      prior_gammaH_me_LOGNORMAL<- 5.94
      prior_gammaH_sd_LOGNORMAL <- 1.80
      ## note tha_T numbers reported in papers are not in log scale, here we convert to
      ## gaussian, since gam is lognormal when log(gam) = gaussian
      simdist<- log(rlnorm(1000000, prior_gammaH_me_LOGNORMAL, prior_gammaH_sd_LOGNORMAL))
      # hist(simdist,breaks=1000) ##can check visually
      prior_gammaH_me<- mean(simdist)
      rm(simdist)

      
      ## taken from taishi paper
      prior_sigmaH_me<- 4
      prior_sigmaH_sd <- 0.51
    


      jj_p<- 1.1
      JUMP_N_p_aRS<- 3 ##how many parameters to jump per MCMC step
      jump_step_K<- 20000*jj_p
      jump_step_eta<- 0.006*jj_p
      jump_step_alpha<- 0.006*jj_p
      jump_step_T<- 3*(1/delta_T)*jj_p #equivalent in days = days*1/delta_T
      jump_step_psi<- 0.006*jj_p
      jump_step_a<- 0.01*jj_p
      jump_step_i<- 0.5*jj_p
      jump_step_delta<- 1*jj_p

      if(STARTLOAD_MCMC){
          load(MCMC_STARTLOAD_FILE_NAME) ##saved_chains
          state<- saved_chains[nrow(saved_chains),]
          rm(saved_chains)
          ##initial proposals
          p_K<- state$K
          p_eta<- state$eta
          p_alpha<- state$alpha
          p_gamH<- state$gamH
          p_sigH<- state$sigH
          p_T<- state$T
          p_psi<- state$psi 
          p_a<- state$a
          p_i<- state$i
          p_delta<- state$delta
          p_data_for_ll<- data.frame()
      }else{
        if(TEMP=="min"){
          ##initial proposals
          p_K<- 3000000
          p_eta<- 2
          p_alpha<- 1
          p_gamH<- 1/prior_gammaH_me
          p_sigH<- 1/prior_sigmaH_me
          p_T<- 3100
          p_psi<- 0.05
          p_a<- prior_a_me
          p_i<- 1
          p_delta<- 0
          p_data_for_ll<- data.frame()
        }else{
          ##initial proposals
          p_K<- 500000
          p_eta<- 1
          p_alpha<- 2
          p_gamH<- 1/prior_gammaH_me
          p_sigH<- 1/prior_sigmaH_me
          p_T<- 3100
          p_psi<- 0.5
          p_a<- prior_a_me
          p_i<- 1
          p_delta<- 0
          p_data_for_ll<- data.frame()
        }
      }

      ##cycle of the MCMC
      MCMCst<- as.integer(1) ##count current MCMC state, integer keeps performance of 'set' call below
      MCMCaccept<- 1 ##we will accept first try anyway
      alikelihood<- -10000 ##previous initial MCMwhich(model_input$timestep>200)[1]C proposal likelihood (keep super low)
      ##saves the accepted results
      saved_chains<- data.frame(delta=rep(0,MCMCstepsN), i=0, K=0, alpha=0, eta=0, gamH=0, sigH=0, T=0, psi=0, a=0, like=0)

      p_Tm <- proc.time() ##allows to measure run time
      repeat{
        
        ##simulate the pre introduction period
        model_init_bf_intro<- c(NH,0,0,0, p_K,0,0,0, 0,0, 0,0,0) ##initial conditions before intro
        ODEmodel_set_params(ts=delta_T, alpha=p_alpha, eta=p_eta, gamH=p_gamH, sigH=p_sigH, a=p_a, f=prior_f, K=p_K, delta=p_delta)
        model_res_bf_intro<- ODEmodel(init=model_init_bf_intro, start=time_start, duration=model_input$timestep[p_T], step_size=delta_T)
        colnames(model_res_bf_intro)<- c('time', "SH","EH","IH","RH",  "AV","SV","EV","IV",  "Inc","IncCS", "R0", "Vls","EIC")
        last_state<- model_res_bf_intro[nrow(model_res_bf_intro),]
        
        model_init_af_intro<- c(NH-p_i,0,p_i,0, last_state$AV,last_state$SV-p_i,last_state$EV+p_i,last_state$IV, p_i,p_i, last_state$R0, last_state$Vls, last_state$EIC) ##initial conditions a_T intro
        ODEmodel_set_params(ts=delta_T, alpha=p_alpha, eta=p_eta, gamH=p_gamH, sigH=p_sigH, a=p_a, f=prior_f, K=p_K, delta=p_delta)
        model_res_af_intro<- ODEmodel(init=model_init_af_intro, start=last_state$time+delta_T, duration=time_stop-(last_state$time+delta_T), step_size=delta_T)
        colnames(model_res_af_intro)<- c('time', "SH","EH","IH","RH",  "AV","SV","EV","IV",  "Inc","IncCS", "R0", "Vls","EIC")

        ##join the pre and post introduction
        model_res<- rbind(model_res_bf_intro, model_res_af_intro)

        ##extract R0
        ##if looking at entire period
        p_R0<- mean(model_res$R0)
        
        ##extract Vls
        Vls<- model_res$Vls
        p_Vls<- mean(Vls)
      

        ###this checks if the total N is OK, if not, then the model has an error
        N <- model_res[nrow(model_res)-1,colnames(model_res) %in% c('SH','EH','IH','RH')] ##used to check solver
        if((sum(N)-NH)>1) stop('ODE did not solve OK. Check N.')

          timeframe<- data.frame(time=model_res$time, timestep=round(model_res$time,1))
          timeframe_matches<- match(timeframe$timestep,reportedcases$timestep)
          FITTING_MATCHES<- which(!is.na(timeframe_matches))
          MODEL_REPORTED_TIMEMATCHES<- timeframe$timestep[FITTING_MATCHES]
          MODEL_REPORTED_TIMEMATCHES<- unique(MODEL_REPORTED_TIMEMATCHES)
          ##get the model values for those time points
          model_CSINC_obs<- model_res$IncCS[c(FITTING_MATCHES[1]-diff(FITTING_MATCHES)[1],FITTING_MATCHES)]* p_psi
          model_INC_obs<- diff(model_CSINC_obs)
          model_CSINC_obs<- model_CSINC_obs[-1]
          EIC<- model_res$EIC[c(FITTING_MATCHES[1]-diff(FITTING_MATCHES)[1],FITTING_MATCHES)][-1]
          VLS<- model_res$Vls[c(FITTING_MATCHES[1]-diff(FITTING_MATCHES)[1],FITTING_MATCHES)][-1]
          predicted<- data.frame(timestep=MODEL_REPORTED_TIMEMATCHES, 
                                 model_INC_obs=model_INC_obs, 
                                 model_CSINC_obs=model_CSINC_obs, 
                                 model_CSINC=model_CSINC_obs*1/p_psi,
                                 model_EIC=EIC,
                                 model_VLS=VLS)
          ##join
          p_data_for_ll<- reportedcases %>% dplyr::left_join(predicted, by="timestep")

        alpha<- 1 ##force first state to be accepted
        plikelihood<- -1000000

        if(MCMCst>1){

          #poisson
          data_probs1<- dpois( round(p_data_for_ll$model_INC_obs), lambda=p_data_for_ll$Cases)          
          data_probs2<- dpois( round(p_data_for_ll$model_CSINC_obs), lambda=p_data_for_ll$CSCases)


          ##gaussian
          prior_probs_sero<- dnorm( tail(p_data_for_ll$model_CSINC,1)/NH, prior_seropositive_me, prior_seropositive_sd)
          prior_probs_R0<- dnorm( p_R0, prior_R0_me, prior_R0_sd)

          
          prior_probs_K<- dweibull( p_K, shape=0.787625929070352, scale=3238597.02247263)

          
          prior_expphase_ti<- 275
          prior_expphase_tf<- 303
          t_period<- p_data_for_ll$timestep>=prior_expphase_ti & p_data_for_ll$timestep<=prior_expphase_tf
          prior_probs_VLSEIC<- p_data_for_ll$model_VLS[t_period]>=p_data_for_ll$model_EIC[t_period]
          prior_probs_VLSEIC<- as.numeric(prior_probs_VLSEIC) #+0.0001

          
          probs<- c(data_probs1, data_probs2, prior_probs_sero, prior_probs_VLSEIC)

          
          plikelihood<- sum(log(probs),na.rm=T)
          ##R has a lot of numerical problems here when likelihoods are too large
          if(is.infinite(plikelihood) & plikelihood<0) plikelihood<- -10000
          alpha<- exp(plikelihood - alikelihood)

          if(is.nan(alpha)) stop('LH was nan...')
        }

        ##before we update the states, show the current and proposed differences

        if(MCMCst %% 200 == 0) {

            mdates<- date_start + p_data_for_ll$timestep
            par(mfrow=c(1,3))
            plot(mdates, p_data_for_ll$Cases, ylim=c(0,500),t ='b', col="black")
            points(mdates, p_data_for_ll$model_INC_obs, col="tomato", pch=20)
            points(mdates, a_data_for_ll$model_INC_obs, col="skyblue", pch=20)

            plot(mdates, p_data_for_ll$CSCases, ylim=c(0,2500),t ='b', col="black")
            points(mdates, p_data_for_ll$model_CSINC_obs, col="tomato", pch=20)
            points(mdates, a_data_for_ll$model_CSINC_obs, col="skyblue", pch=20)

            plot(date_start + model_res$time, model_res$EIC, t ='l', col="black", ylim=c(0,20))
            points(date_start + model_res$time, model_res$Vls, col="tomato", pch=20, t='l')

        }

        ##accept vs not accept
        if(runif(1)<alpha){
          
          ##accept this state
          
          ##saving into data.frames is very slow and becomes slower as MCMC grows
          ##save by reference using data.table; performace increases EXP with MCMCsteps
          set(saved_chains, i=MCMCst, j=1:ncol(saved_chains), value=list(delta=p_delta, i=p_i, K=p_K, alpha=p_alpha, eta=p_eta, gamH=p_gamH, sigH=p_sigH, T=p_T, psi=p_psi, a=p_a, like=plikelihood))
          
          ##make proposals the accepted state
          a_K<- p_K 
          a_alpha<- p_alpha 
          a_eta<- p_eta 
          a_gamH<- p_gamH 
          a_sigH<- p_sigH 
          a_T<- p_T 
          a_psi<- p_psi
          a_a<- p_a
          a_i<- p_i
          a_delta<- p_delta
          ##also likelihood
          alikelihood<- plikelihood
          a_data_for_ll<- p_data_for_ll
          
          ##make jumps: we do not jump all parameters per step
          sample_par<- sample(1:5,JUMP_N_p_aRS)
          if(1 %in% sample_par) p_K<- round(sampleDistGaussian(mean=p_K, sd=jump_step_K, min=min(prior_range_K), max=max(prior_range_K)))
          if(2 %in% sample_par) p_eta<- (sampleDistGaussian(mean=p_eta, sd=jump_step_eta, min=min(prior_range_eta), max=max(prior_range_eta)))
          if(3 %in% sample_par) p_T<- round(sampleDistGaussian(mean=p_T, sd=jump_step_T, min=min(prior_range_T), max=max(prior_range_T)))
          if(4 %in% sample_par) p_psi<- (sampleDistGaussian(mean=p_psi, sd=jump_step_psi, min=min(prior_range_psi), max=max(prior_range_psi)))
          if(5 %in% sample_par) p_alpha<- (sampleDistGaussian(mean=p_alpha, sd=jump_step_alpha, min=min(prior_range_alpha), max=max(prior_range_alpha)))

          MCMCaccept<- MCMCaccept+1
        }else{
          
          ##reject this state >> accept previous
          
          ##saving into data.frames is very slow and becomes slower as MCMC grows
          ##save by reference using data.table; performace increases EXP with MCMCsteps
          set(saved_chains, i=MCMCst, j=1:ncol(saved_chains), value=list(delta=a_delta, i=a_i, K=a_K, alpha=a_alpha, eta=a_eta, gamH=a_gamH, sigH=a_sigH, T=a_T, ps=a_psi, a=a_a, like=alikelihood))
          
          ##force all next proposals to be current values by default (in case not jumped below)
          p_K<- a_K 
          p_alpha<- a_alpha
          p_eta<- a_eta 
          p_gamH<- a_gamH 
          p_sigH<- a_sigH 
          p_T<- a_T 
          p_psi<- a_psi
          p_a<- a_a
          p_i<- a_i
          p_delta<- a_delta
          p_data_for_ll<- a_data_for_ll
          
          ##make jumps: we do not jump all parameters per step
          sample_par<- sample(1:5,JUMP_N_p_aRS)
          if(1 %in% sample_par) p_K<- round(sampleDistGaussian(mean=p_K, sd=jump_step_K, min=min(prior_range_K), max=max(prior_range_K)))
          if(2 %in% sample_par) p_eta<- (sampleDistGaussian(mean=p_eta, sd=jump_step_eta, min=min(prior_range_eta), max=max(prior_range_eta)))
          if(3 %in% sample_par) p_T<- round(sampleDistGaussian(mean=p_T, sd=jump_step_T, min=min(prior_range_T), max=max(prior_range_T)))
          if(4 %in% sample_par) p_psi<- (sampleDistGaussian(mean=p_psi, sd=jump_step_psi, min=min(prior_range_psi), max=max(prior_range_psi)))
          if(5 %in% sample_par) p_alpha<- (sampleDistGaussian(mean=p_alpha, sd=jump_step_alpha, min=min(prior_range_alpha), max=max(prior_range_alpha)))
          
        }
        
        acceptrate<- round(MCMCaccept/MCMCst,2) ##calculate acceptance rate in real time
        if(MCMCst>= MCMCstepsN) break;

        if(MCMCst %% 200 == 0) {
          print(paste(">>> %",round(100*MCMCst/MCMCstepsN,5), "AR", acceptrate))
          print(paste("R0",round(p_R0,3),"Vls",round(p_Vls,3),"delta",p_delta,"i",round(p_i,3),"K",round(p_K,3), "aV",round(p_a,3), "alp",round(p_alpha,3), "eta",round(p_eta,3), "T",round(p_T,3),"[",model_res$time[p_T]+date_start,"]" , "psi",round(p_psi,3), "[ sero", round(tail(predicted$model_CSINC,1)/NH,4), "]" ))
        }
        
        MCMCst<- as.integer(MCMCst+1) ##integer keeps performance of using 'set' call above
        
      }
      gc() ##garbage collection
      print(proc.time() - p_Tm)

      save(saved_chains, file=paste0(MCMC_OUTPUT_FILE_NAME,"_fullchains.Rdata"))

}

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

