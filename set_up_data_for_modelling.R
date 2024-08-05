#### Please do not run this file directly as it will not work.
#### Please refer to file 00_run.R

##load data

  df_inter <- read.csv(DATA_FILE_NAME)
  df_inter$days<- 1:nrow(df_inter)
  head(df_inter)

  TOTAL_CASES_REPORTED<- sum(df_inter$Cases, na.rm=T)

### interpolate the climate to match the model's step size

  interpolator_temperature<- approxfun(df_inter$days, df_inter$Temp, rule = 2)
  new_temperature<- interpolator_temperature(tps)
  mean_temperature<- mean(new_temperature, na.rm=TRUE)

  ##normalize hum + prec to 0-1
  df_inter$Hum<- df_inter$Hum/max(df_inter$Hum, na.rm=TRUE)
  interpolator_humidity<- approxfun(df_inter$days, df_inter$Hum, rule = 2)
  new_humidity<- interpolator_humidity(tps)
  mean_humidity<- mean(new_humidity, na.rm=TRUE)

  ##normalize hum + prec to 0-1
  df_inter$Prec<- df_inter$Prec/max(df_inter$Prec, na.rm=TRUE)
  interpolator_precipitation<- approxfun(df_inter$days, df_inter$Prec, rule = 2)
  new_precipitation<- interpolator_precipitation(tps)
  mean_precipitation<- mean(new_precipitation, na.rm=TRUE)

  if(SMOOTH_CLIMATE){
    
    x<- SMOOTH_DAYS #7 #number of days to smooth around
    n<- ((1/delta_T))*x

    sides<- 2 ## 1= into past, 2=centered
		temp<- stats::filter(new_temperature, rep(1/n, n), sides=sides)
		hum<- stats::filter(new_humidity, rep(1/n, n), sides=sides)
		prec<- stats::filter(new_precipitation, rep(1/n, n), sides=sides)

    ##because we loose data on the start or/and end, just replace it with original values
    temp[which(is.na(temp))]<- new_temperature[which(is.na(temp))]
    hum[which(is.na(hum))]<- new_humidity[which(is.na(hum))]
    prec[which(is.na(prec))]<- new_precipitation[which(is.na(prec))]

    # pdf(paste0(MCMC_OUTPUT_FILE_NAME,"_smoothed_climate.pdf"), w=5, h=3)

    plot(12*tps/365, new_temperature, col="red", t='l', main="SMOOTHED")
    lines(12*tps/365,temp, col="black", lw=2)

    plot(12*tps/365,new_humidity, col="red", t='l', main="SMOOTHED")
    lines(12*tps/365,hum, col="black", lw=2)

    plot(12*tps/365,new_precipitation, col="red", t='l', main="SMOOTHED")
    lines(12*tps/365, prec, col="black", lw=2)

    dev.off()

    new_temperature<- temp
    new_humidity<- hum
    new_precipitation<- prec

  }else{

    # pdf(paste0(MCMC_OUTPUT_FILE_NAME,"_smoothed_climate.pdf"), w=5, h=3)

    plot(12*tps/365, new_temperature, col="red", t='l', main="NOT SMOOTHED")

    plot(12*tps/365,new_humidity, col="red", t='l', main="NOT SMOOTHED")

    plot(12*tps/365,new_precipitation, col="red", t='l', main="NOT SMOOTHED")

    dev.off()

  }

## generate and interpolate all of the parameters

  C_epsilon <- epsilon_function(new_temperature)
  C_mu_A <- mu_A_function(new_temperature)
  C_theta <- theta_function(new_temperature)
  C_gamma_V <- gamma_V_function(new_temperature) ##
  C_phi_HV <- phi_HV_function(new_temperature)
  C_phi_VH <- phi_VH_function(new_temperature)
  C_mu_V_temp <- mu_V_temp_function(new_temperature)
  C_c_temp <- c_temp_function(new_temperature)
  C_mu_V_hum <-mu_V_hum_function(new_humidity, mean_humidity) 
  C_a_hum <- a_hum_function(new_humidity, mean_humidity)
  C_c_rain <- c_rain_function(new_precipitation, mean_precipitation)
  mu_A_rain <- mu_A_prec_function(new_precipitation, mean_precipitation)
  gam_V_temp <- gamV_temp_function(new_temperature, mean_temperature) ##


## create model input data

  model_input<- data.frame(timestep=tps, 
                          temp=new_temperature, 
                          hum=new_humidity, 
                          prec=new_precipitation,
                          epsilon= C_epsilon,
                          mu_A= C_mu_A,
                          theta= C_theta,
                          gamma_V= C_gamma_V,
                          phi_HV= C_phi_HV,
                          phi_VH= C_phi_VH,
                          mu_V_temp= C_mu_V_temp,
                          mu_V_hum= C_mu_V_hum,
                          c_temp= C_c_temp,
                          c_rain= C_c_rain,
                          mu_A_rain= mu_A_rain,
                          a_rain= C_c_rain,
                          gam_V_temp=gam_V_temp,
                          a_hum=C_a_hum)

  model_input<- round(model_input,4)
  head(model_input)
  
  write.csv(model_input, file="data/model_input.csv", row.names=FALSE)

  ##curate the data to be fitted

  reportedcases<- df_inter %>% dplyr::select(days, Cases)
  reportedcases<- reportedcases %>% dplyr::filter(!is.na(Cases))
  reportedcases<- reportedcases %>% dplyr::rename(timestep="days")
  reportedcases$CSCases<- cumsum(reportedcases$Cases)
