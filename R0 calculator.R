#### This code runs the R0 calculations from 2002-2013.
#### Lines 5 to 140 run the model. Select the year you wish to do this for to run this section of code.
#### The input and output files have already been created and saved using this code if you wish to skip this section.
#### Lines 143 to 285 load these previously saved inputs and outputs and creates the graphs seen in the thesis and supplementary information.
#### For the individual year graphs, change the year written in lines 255, 256 and 257 to the desired year (2002-2013)


source("parameter_functions.R")


##load data
DATA_FILE <- '2002'
DATA_FILE <- '2003'
DATA_FILE <- '2004'
DATA_FILE <- '2005'
DATA_FILE <- '2006'
DATA_FILE <- '2007'
DATA_FILE <- '2008'
DATA_FILE <- '2009'
DATA_FILE <- '2010'
DATA_FILE <- '2011'
DATA_FILE <- '2012'
DATA_FILE <- '2013'

df_inter <- read.csv(paste0('Data/Climate R0/',DATA_FILE,'.csv'))
df_inter$days<- 1:nrow(df_inter)
head(df_inter)


time_start <- 0
time_stop <- max(df_inter$days) - 1
delta_T<- 0.0625 
tps <- seq(time_start , time_stop , by = delta_T)
length(tps)

SMOOTH_CLIMATE<- TRUE
SMOOTH_DAYS<- 7

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
if (max(df_inter$Prec, na.rm = TRUE) > 0){
    df_inter$Prec<- df_inter$Prec/max(df_inter$Prec, na.rm=TRUE)
    print("done")
}
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
  
  plot(12*tps/365, new_temperature, col="red", t='l', main="SMOOTHED")
  lines(12*tps/365,temp, col="black", lw=2)
  
  plot(12*tps/365,new_humidity, col="red", t='l', main="SMOOTHED")
  lines(12*tps/365,hum, col="black", lw=2)
  
  plot(12*tps/365,new_precipitation, col="red", t='l', main="SMOOTHED")
  lines(12*tps/365, prec, col="black", lw=2)
  
  new_temperature<- temp
  new_humidity<- hum
  new_precipitation<- prec
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
mu_A_rain <- mu_A_prec_function(new_precipitation, mean_precipitation)



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
                         mu_A_rain= mu_A_rain)

model_input<- round(model_input,4)
# write.csv(model_input, file=paste0('Data/R0 model input/model_input_', DATA_FILE, '.csv'), row.names=FALSE)

NH<- 270000
p_alpha <- 2.480908
p_eta <- 1.657568
p_gamH <- 1/5.93686
p_sigH <- 1/4
p_a <- 0.25
prior_f <- 0.5
p_K <- 1753862
p_delta <- 0

ODEmodel_set_params(ts=delta_T, alpha=p_alpha, eta=p_eta, gamH=p_gamH, sigH=p_sigH, a=p_a, f=prior_f, K=p_K, delta=p_delta)
model_starting <- c(NH,0,0,0, p_K,0,0,0, 0,0, 0,0,0)
model_output<- ODEmodel(init=model_starting, start=time_start, duration = max(model_input$timestep), step_size=delta_T)
colnames(model_output)<- c('time', "SH","EH","IH","RH",  "AV","SV","EV","IV",  "Inc","IncCS", "R0", "Vls","EIC")

# write.csv(model_output, file=paste0('Outputs/model_output_', DATA_FILE, '.csv'), row.names=FALSE)

##################################################################################################################################

model_input_file_2002 <- read.csv(paste0('Data/R0 model input/model_input_2002.csv'))
model_input_file_2003 <- read.csv(paste0('Data/R0 model input/model_input_2003.csv'))
model_input_file_2004 <- read.csv(paste0('Data/R0 model input/model_input_2004.csv'))
model_input_file_2005 <- read.csv(paste0('Data/R0 model input/model_input_2005.csv'))
model_input_file_2006 <- read.csv(paste0('Data/R0 model input/model_input_2006.csv'))
model_input_file_2007 <- read.csv(paste0('Data/R0 model input/model_input_2007.csv'))
model_input_file_2008 <- read.csv(paste0('Data/R0 model input/model_input_2008.csv'))
model_input_file_2009 <- read.csv(paste0('Data/R0 model input/model_input_2009.csv'))
model_input_file_2010 <- read.csv(paste0('Data/R0 model input/model_input_2010.csv'))
model_input_file_2011 <- read.csv(paste0('Data/R0 model input/model_input_2011.csv'))
model_input_file_2012 <- read.csv(paste0('Data/R0 model input/model_input_2012.csv'))
model_input_file_2013 <- read.csv(paste0('Data/R0 model input/model_input_2013.csv'))

mean_input_2002 <- colMeans(model_input_file_2002)
mean_input_2003 <- colMeans(model_input_file_2003)
mean_input_2004 <- colMeans(model_input_file_2004)
mean_input_2005 <- colMeans(model_input_file_2005)
mean_input_2006 <- colMeans(model_input_file_2006)
mean_input_2007 <- colMeans(model_input_file_2007)
mean_input_2008 <- colMeans(model_input_file_2008)
mean_input_2009 <- colMeans(model_input_file_2009)
mean_input_2010 <- colMeans(model_input_file_2010)
mean_input_2011 <- colMeans(model_input_file_2011)
mean_input_2012 <- colMeans(model_input_file_2012)
mean_input_2013 <- colMeans(model_input_file_2013)

mean_input_parameters <- data.frame(mean_input_2002, mean_input_2003, mean_input_2004, mean_input_2005,mean_input_2006, mean_input_2007, mean_input_2008, mean_input_2009, mean_input_2010, mean_input_2011, mean_input_2012, mean_input_2013)


model_output_file_2002 <- read.csv(paste0('Outputs/model_output_2002.csv'))
model_output_file_2003 <- read.csv(paste0('Outputs/model_output_2003.csv'))
model_output_file_2004 <- read.csv(paste0('Outputs/model_output_2004.csv'))
model_output_file_2005 <- read.csv(paste0('Outputs/model_output_2005.csv'))
model_output_file_2006 <- read.csv(paste0('Outputs/model_output_2006.csv'))
model_output_file_2007 <- read.csv(paste0('Outputs/model_output_2007.csv'))
model_output_file_2008 <- read.csv(paste0('Outputs/model_output_2008.csv'))
model_output_file_2009 <- read.csv(paste0('Outputs/model_output_2009.csv'))
model_output_file_2010 <- read.csv(paste0('Outputs/model_output_2010.csv'))
model_output_file_2011 <- read.csv(paste0('Outputs/model_output_2011.csv'))
model_output_file_2012 <- read.csv(paste0('Outputs/model_output_2012.csv'))
model_output_file_2013 <- read.csv(paste0('Outputs/model_output_2013.csv'))

#########################################################################################################################################################################################
#### Mean R0 results

mean_R0 <- c(mean(model_output_file_2002$R0),
                  mean(model_output_file_2003$R0),
                  mean(model_output_file_2004$R0),
                  mean(model_output_file_2005$R0),
                  mean(model_output_file_2006$R0),
                  mean(model_output_file_2007$R0),
                  mean(model_output_file_2008$R0),
                  mean(model_output_file_2009$R0),
                  mean(model_output_file_2010$R0),
                  mean(model_output_file_2011$R0),
                  mean(model_output_file_2012$R0),
                  mean(model_output_file_2013$R0))
year <- c(2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013)

mean_data_R0 <- data.frame(year, mean_R0)

plot(mean_data_R0)
  
mean_input_parameters <- rbind(mean_input_parameters, mean_R0)
row.names(mean_input_parameters)[15] <- "R0"

t_params <- t(mean_input_parameters)
R0_all <- as.data.frame(cbind(t_params, year))

year <- R0_all$year
mean_R0 <- R0_all$R0
temp <- R0_all$temp
hum <- R0_all$hum
prec <- R0_all$prec

{climate = list(temp,hum,prec)
colors = c("#f7ae60", "#f760c5", "#6e6ef0")
par(oma = c(0, 2, 2, 3), cex.axis = 0.85)

plot(year, mean_R0, t='b', col="black", ylim=c(2.5,2.65), ylab = "")

sides <- list(2, 4, 4) 
lines <- list(2, NA, 2)

par(new = TRUE)
plot(year, climate[[1]], axes = FALSE, col = "#f7ae60", xlab = "", ylab = "", type = "b")
axis(at = pretty(climate[[1]]), side = 2, line = 2, col = "#f7ae60")


par(new = TRUE)
plot(year, climate[[2]], axes = FALSE, col = "#f760c5", xlab = "", ylab = "", type = "b")
axis(at = pretty(climate[[2]]), side = 4, line = NA, 
     col = "#f760c5")


par(new = TRUE)
plot(year, climate[[3]], axes = FALSE, col = "#6e6ef0", xlab = "", ylab = "", type = "b")
axis(at = pretty(climate[[3]]), side = 4, line = 2, 
     col = "#6e6ef0")


par(new = TRUE)
plot(year, mean_R0, t='b', col="black", ylim=c(2.5,2.65), ylab = "", xlab = "")


legend("topleft", inset=.02, legend=c("R0", "Temperature (Celcius)", "Humidity (%)", "Precipitation (mm)"), col=c("black", "#f7ae60", "#f760c5", "#6e6ef0"), lty=1, cex=0.65)}

#########################################################################################################################################################################################
#### R0 yearly trend
year_to_plot <- 2013
input_year <- model_input_file_2013
output_year <- model_output_file_2013

#### R0 and temperature plot
date_start<- as.Date(paste0(year_to_plot, "-01-01"),"%Y-%m-%d") -1
mdates<- date_start + input_year$timestep
par(mar=c(5, 4, 4, 5) + 0.1)
plot(mdates, output_year$R0, t='l', col="black", xlab = "Date", ylab = "Reproductive number (R0)", main = year_to_plot)
abline(h=1,col="#4d4deb", lty = 2)
par(new=TRUE)
plot(mdates, input_year$temp, col="#971fcc", t='l', ylab = "", xlab = "", axes = FALSE)
axis(side=4)
mtext("Temperature (Celcius)", side=4, line=3)
legend("topleft", inset=.04,legend=c("R0", "Temperature", "R0 of 1", "Introduction date"),
       col=c("black", "#971fcc", "#4d4deb", "tomato"), lty=c(1,1,2), cex=0.7)

#### Humidity and precipitation plot

par(mar=c(5, 4, 4, 5) + 0.1)
plot(mdates, input_year$hum*100, t='l', col="#f7ae60", xlab = "Date", ylab = "Humidity (%)", main = year_to_plot)
abline(h=1,col="#4d4deb")
par(new=TRUE)
plot(mdates,input_year$prec, col="#f760c5", t='l', ylab = "", xlab = "", axes = FALSE)
axis(side=4)
mtext("Precipitation (mm)", side=4, line=3)
legend("topleft", inset=.04,legend=c("Humidity", "Precipitation", "Introduction date"),
       col=c("#f7ae60", "#f760c5", "tomato"), lty=1, cex=0.7)

#################################################################################################################################
