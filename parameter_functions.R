#### Please do not run this file directly as it will not work.
#### Please refer to file 00_run.R

#########################################################################

epsilon_function <- function(x) { ##checked = Madeira article
  eps <- 0.131 - (0.05723 * x) + (0.01164 * x^2) - (0.001341 * x^3) + (0.00008723 * x^4) - (0.000003017 * x^5) + ((5.153*10^-8) * x^6) - ((3.42*10^-10) * x^7)
  return (eps)
} 

#########################################################################

mu_A_function <- function(x){ ##checked = Madeira article
  muA <- 2.13 - (0.3797 * x) + (0.02457 * x^2) - (0.0006778 * x^3) + (0.000006794 * x^4)
  return(muA)
}

# x1<- seq(1,40) #model_input$temp
# x2<- seq(1,100,length.out=length(x1))/100 #model_input$prec
# xx<- expand_grid(x1,x2)
# x1<- xx$x1
# x2<- xx$x2
# mu_A_prec_function1 <- function(w, v){  ###FEIRA SANTANA
#   # muV_h <- v - (w - v) / sqrt(1 + (w - v)^2) 
#   muV_h <- (w - v) / sqrt(1 + (w - v)^2)
#   return(muV_h)
# }
# mu_A_prec_function2 <- function(w, v){  ###FEIRA SANTANA
#   muV_h <- v - (w - v) / sqrt(1 + (w - v)^2) 
#   # muV_h <- (w - v) / sqrt(1 + (w - v)^2)
#   return(muV_h)
# }
# df<- data.frame()
# for(rho in c(0.5,1,3)){
#   y1<- mu_A_function(x1) * (1 + mu_A_prec_function1(x2, mean(x2)))^rho
#   y2<- mu_A_function(x1) * (1 + mu_A_prec_function2(x2, mean(x2)))^rho
#   df<- rbind(data.frame(x1=x1, x2=x2, y1=y1, y2=y2, rho=rho), df)
# }
# ggplot(df %>% filter(rho==0.5)) + geom_point(aes(x=x1, y=x2, fill=1/y1, color=1/y1), pch=21) +theme_classic()
# ggplot(df %>% filter(rho==0.5)) + geom_point(aes(x=x1, y=x2, fill=1/y2, color=1/y2), pch=21) +theme_classic()
# ggplot(df %>% filter(rho==3)) + geom_point(aes(x=x1, y=x2, fill=1/y1, color=1/y1), pch=21) +theme_classic()
# ggplot(df %>% filter(rho==3)) + geom_point(aes(x=x1, y=x2, fill=1/y2, color=1/y2), pch=21) +theme_classic()


#########################################################################

mu_V_temp_function <- function(x){ ##checked = Madeira article
  muV_t <- 0.8692 - (0.1599 *x) + (0.01116 * x^2) - (0.0003408 * x^3) + (0.000003809 * x^4)
  return(muV_t)
}


# x1<- seq(1,40) #model_input$temp
# x2<- seq(1,100,length.out=length(x1))/100 #model_input$prec
# xx<- expand_grid(x1,x2)
# x1<- xx$x1
# x2<- xx$x2
# mu_V_hum_function1 <- function(w, v){  ###FEIRA SANTANA
#   # muV_h <- v - (w - v) / sqrt(1 + (w - v)^2) 
#   muV_h <- (w - v) / sqrt(1 + (w - v)^2)
#   return(muV_h)
# }
# mu_V_hum_function2 <- function(w, v){  ###FEIRA SANTANA
#   muV_h <- v - (w - v) / sqrt(1 + (w - v)^2) 
#   # muV_h <- (w - v) / sqrt(1 + (w - v)^2)
#   return(muV_h)
# }
# df<- data.frame()
# for(rho in c(0.5,1,3)){
#   y1<- mu_V_temp_function(x1) * (1 + mu_V_hum_function1(x2, mean(x2)))^rho
#   y2<- mu_V_temp_function(x1) * (1 + mu_V_hum_function2(x2, mean(x2)))^rho
#   y1[y1<0]<- 0
#   y2[y2<0]<- 0
#   df<- rbind(data.frame(x1=x1, x2=x2, y1=y1, y2=y2, rho=rho), df)
# }
# ggplot(df %>% filter(rho==0.5)) + geom_point(aes(x=x1, y=x2, fill=1/y1, color=1/y1), pch=21) +theme_classic()
# ggplot(df %>% filter(rho==0.5)) + geom_point(aes(x=x1, y=x2, fill=1/y2, color=1/y2), pch=21) +theme_classic()
# ggplot(df %>% filter(rho==3)) + geom_point(aes(x=x1, y=x2, fill=1/y1, color=1/y1), pch=21) +theme_classic()
# ggplot(df %>% filter(rho==3)) + geom_point(aes(x=x1, y=x2, fill=1/y2, color=1/y2), pch=21) +theme_classic()


#########################################################################
#########################################################################

theta_function <- function(x){ ##checked = Madeira article
  the <- -5.4 + (1.8 * x) - (0.2124 * x^2) + (0.01015 * x^3) - (0.0001515 *x^4)
  return(the)
}

#########################################################################

gamma_V_function <- function(x){

    gamV<- rep(NA, length(x))
    gamV[which(x<15)]<- 1/1228
    gamV[which(x >= 15 & x < 17.5)]<- 1/232
    gamV[which(x >= 17.5 & x < 20)]<- 1/72.4
    gamV[which(x >= 20 & x < 22.5)]<- 1/28.9
    gamV[which(x >= 22.5 & x < 25)]<- 1/14.7
    gamV[which(x >= 25 & x < 27.5)]<- 1/8.68
    gamV[which(x >= 27.5 & x < 30)]<- 1/5.76
    gamV[which(x >= 30 & x < 32.5)]<- 1/4.14
    gamV[which(x >= 32.5)]<- 1/3.19

    if(sum(is.na(gamV))) stop("gamV error.")

    return(gamV)
}

##old gamma_V function
MADEIRA_temp_effect_gammaV<- function(T){ ##checked = Madeira article
  Tk<- T+273.15 
  R<- 1.987
  ef<- (24.0*( 0.003359* (Tk/298.) * exp((15000./R)*(1/298.-1./Tk)) / (1.+ exp((6.203*(10^21)/R)*(1./(-2.176*(10^30))-1./Tk))) ))
  ef[which(ef<0)]<- 0 #fix negative numbers as bio -> trim to zero
  return(ef)
}

##set the one we want = the original one from Madeira and Feira de Santana
gamma_V_function<- MADEIRA_temp_effect_gammaV

# ##this is just an approximation for Madeira
# MADEIRA_temp_effect_gammaV<- function(x){
#  return((4+exp(5.15-0.123*x)))
# }

#########################################################################


phi_HV_function<- function(x){ ##checked = Madeira + FSA article
    pHV<- (0.001044*x)*(x-12.286)*(32.461-x)^(1/2)
    return(pHV)
}

#########################################################################

phi_VH_function<- function(x){ ##Madeira article
  pVH<- 0.0729*x - 0.97
  return(pVH)
}

#########################################################################
## these are true additions to the model, since they were not modelled before for MADEIRA

c_temp_function <- function(x){ ###FEIRA SANTANA
  c_t_calc <- (-184.8 + (27.94 * x) - (0.9254 * x^2) + (0.009226 * x^3)) / 100
  return(c_t_calc)
}


###

c_rain_function <- function(y, z){ ###FEIRA SANTANA
  # c_r_calc <- z - (y - z) / sqrt(1 + (y - z)^2) ##reduces parameter when climate var is higher
  c_calc <- (y - z) / sqrt(1 + (y - z)^2) ##increases parameter when climate var is higher
  return(c_calc)
}

a_hum_function <- function(w, v){ ###FEIRA SANTANA
  # a_h_calc <- v - (w - v) / sqrt(1 + (w - v)^2) ##reduces parameter when climate var is higher
  a_calc <- (w - v) / sqrt(1 + (w - v)^2) ##increases parameter when climate var is higher
  return(a_calc)
}

# a_rain_function <- function(w, v){ ###FEIRA SANTANA
#   # a_h_calc <- v - (w - v) / sqrt(1 + (w - v)^2) ##reduces parameter when climate var is higher
#   a_calc <- (w - v) / sqrt(1 + (w - v)^2) ##increases parameter when climate var is higher
#   return(a_calc)
# }

mu_V_hum_function <- function(w, v){  ###FEIRA SANTANA
  muV_h <- v - (w - v) / sqrt(1 + (w - v)^2) ##reduces parameter when climate var is higher
  # muV_h <- (w - v) / sqrt(1 + (w - v)^2) ##increases parameter when climate var is higher
  return(muV_h)
}

mu_A_prec_function <- function(w, v){  ###FEIRA SANTANA
  muV_h <- v - (w - v) / sqrt(1 + (w - v)^2) ##reduces parameter when climate var is higher
  # muV_h <- (w - v) / sqrt(1 + (w - v)^2) ##increases parameter when climate var is higher
  return(muV_h)
}


gamV_temp_function <- function(y, z){ ###NEW
    # c_calc <- z - (y - z) / sqrt(1 + (y - z)^2) ##reduces parameter when climate var is higher
    c_calc <- (y - z) / sqrt(1 + (y - z)^2) ##increases parameter when climate var is higher
    return(c_calc)
}
