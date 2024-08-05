
#### This code runs the spatial analysis including the MFA, spatial autocorrelation and spatial cross-correlation
#### The only required packages are listed at the top here:

library(FactoMineR)
library(factoextra)
library(tidyverse)
library(sf)
library(sp)
library(spdep)
library(spData)
library(spatialEco)
library(tmap)
library(tmaptools)
library(spgwr)


###################################################################################################################################################################################################################################################################################################################################################
## MFA http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/116-mfa-multiple-factor-analysis-in-r-essentials/


spatial_data <- read.csv("Data/Spatial data.csv")
spatial_data <- na.omit(spatial_data)

spatial_data[,1] <- as.character(spatial_data[,1])
spatial_data[,2] <- as.character(spatial_data[,2])
spatial_data[,3] <- as.character(spatial_data[,3])
spatial_data[,4] <- as.character(spatial_data[,4])

mfa_spatial <- MFA(spatial_data, group = c(4,1,1,8,91), type = c("n", "s", "s", "s", "s"), name.group = c("Location", "CIR", "Population", "Land type", "Climate"), num.group.sup = c(1), graph = FALSE)

eig_val_mfa <- get_eigenvalue(mfa_spatial)

mfa_spatial_more <- MFA(spatial_data, group = c(4,1,1,8,13,13,13,13,13,13,13), type = c("n", "s", "s", "s", "s", "s", "s", "s","s", "s", "s"), name.group = c("Location", "CIR", "Population", "Land type", "Moisture", "Humidity", "Precipitation", "Wind", "Max temp", "Min temp", "Mean temp"), num.group.sup = c(1), graph = FALSE)

fviz_mfa_var(mfa_spatial_more, "quanti.var", geom = "arrow")
fviz_mfa_var(mfa_spatial_more, "quanti.var", geom = c("point"))
fviz_mfa_var(mfa_spatial_more, "quanti.var", geom = "text")
fviz_mfa_var(mfa_spatial_more, "quanti.var")

###################################################################################################################################################################################################################################################################################################################################################
##Spatial correlation https://rpubs.com/quarcs-lab/spatial-autocorrelation

shape_file <- read_sf("Data/Madeira/MADEIRA_adm3.shp")
non_shape_file <- na.omit(read.csv("Data/Spatial data.csv"))
merged_file <- merge(shape_file, non_shape_file, by.x="ID_3", by.y = "Parish..ID.3.")
merged_file_sf <- st_as_sf(merged_file)

bbox_new <- st_bbox(merged_file_sf) # current bounding box

xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values


bbox_new[3] <- bbox_new[3] + (0.25 * xrange) # xmax - right

bbox_new <- bbox_new %>%  # take the bounding box ...
  st_as_sfc() # ... and make it a sf polygon

tm_shape(merged_file_sf, bbox = bbox_new) + 
  tm_fill("Cases..mean.",
          palette = "Reds", 
          style = "cat",
          n = 6,
          title = "") +
  tm_borders(alpha=.4)  +
  tm_layout(legend.position = c(0.8,0.1),
            legend.height = 0.7,
            title = "CIR",
            title.position = c(0.81, 0.83),
            title.size = 0.85)

tm_shape(merged_file_sf, bbox = bbox_new) + 
  tm_fill("Urban",
          palette = "Blues", 
          style = "equal", 
          title = "",
          n = 6) +
  tm_borders(alpha=.4)  +
  tm_layout(legend.position = c(0.8,0.1),
            legend.height = 0.7,
            title = "Urban Land",
            title.position = c(0.81, 0.83),
            title.size = 0.85)

tm_shape(merged_file_sf, bbox = bbox_new) + 
  tm_fill("maxtemp_ago_mean",
          palette = "Greens", 
          style = "equal", 
          title = "") +
  tm_borders(alpha=.4)  +
  tm_layout(legend.position = c(0.8,0.1),
            legend.height = 0.7,
            title = "Max Temp Aug",
            title.position = c(0.81, 0.83),
            title.size = 0.85)

tm_shape(merged_file_sf, bbox = bbox_new) + 
  tm_fill("popden_mean",
          palette = "Purples", 
          style = "equal",
          n = 6,
          title = "") +
  tm_borders(alpha=.4)  +
  tm_layout(legend.position = c(0.8,0.1),
            legend.height = 0.7,
            title = "Population Density",
            title.position = c(0.81, 0.83),
            title.size = 0.7)

tm_shape(merged_file_sf, bbox = bbox_new) + 
  tm_fill("maxtemp_jul_mean",
          palette = "Greens", 
          style = "equal",
          n = 4,
          title = "") +
  tm_borders(alpha=.4)  +
  tm_layout(legend.position = c(0.8,0.1),
            legend.height = 0.7,
            title = "Max Temp July",
            title.position = c(0.81, 0.83),
            title.size = 0.85)

tm_shape(merged_file_sf, bbox = bbox_new) + 
  tm_fill("hurs_may_mean",
          palette = "Greens", 
          style = "equal", 
          n = 4,
          title = "") +
  tm_borders(alpha=.4)  +
  tm_layout(legend.position = c(0.8,0.1),
            legend.height = 0.7,
            title = "Mean Humidity May",
            title.position = c(0.81, 0.83),
            title.size = 0.7)

tm_shape(merged_file_sf, bbox = bbox_new) + 
  tm_fill("fcover_mean",
          palette = "Blues", 
          style = "quantile", 
          n = 5,
          title = "") +
  tm_borders(alpha=.4)  +
  tm_layout(legend.position = c(0.8,0.1),
            legend.height = 0.7,
            title = "Forest Cover",
            title.position = c(0.81, 0.83),
            title.size = 0.85)

tm_shape(merged_file_sf, bbox = bbox_new) + 
  tm_fill("cmi_ago_mean",
          palette = "Greens", 
          style = "quantile", 
          n = 5,
          title = "") +
  tm_borders(alpha=.4)  +
  tm_layout(legend.position = c(0.8,0.1),
            legend.height = 0.7,
            title = "Mean Moisture Aug",
            title.position = c(0.81, 0.83),
            title.size = 0.7)

tm_shape(merged_file_sf, bbox = bbox_new) + 
  tm_fill("Wall",
          palette = "Greens", 
          style = "quantile",
          n = 5,
          title = "") +
  tm_borders(alpha=.4)  +
  tm_layout(legend.position = c(0.8,0.1),
            legend.height = 0.7,
            title = "Mean wind",
            title.position = c(0.81, 0.83),
            title.size = 0.85)

###### make queen neighbourhoods
neighbours_queen <- poly2nb(merged_file)

{plot(st_geometry(merged_file), border = 'lightgrey')
plot.nb(neighbours_queen, st_geometry(merged_file), add = TRUE, col = "red")
title(main="Queen Neighbourhoods",line = -2)}

##### calculate morans I and p value
weights_def <- nb2listw(neighbours_queen)
globalMoran_cases <- moran.test(merged_file$Cases..mean., weights_def)
globalMoran_urban <- moran.test(merged_file$Urban, weights_def)
globalMoran_max_aug <- moran.test(merged_file$maxtemp_ago_mean, weights_def)
globalMoran_population <- moran.test(merged_file$popden_mean, weights_def)
globalMoran_max_july <- moran.test(merged_file$maxtemp_jul_mean, weights_def)
globalMoran_hum_may <- moran.test(merged_file$hurs_may_mean, weights_def)
globalMoran_forest <- moran.test(merged_file$fcover_mean, weights_def)
globalMoran_moisture_aug <- moran.test(merged_file$cmi_ago_mean, weights_def)
globalMoran_wind <- moran.test(merged_file$Wall, weights_def)


###################################################################################################################################################################################################################################################################################################################################################
#### spatial cross-correlation

MORAN_res<- c()
MAPS_res<- list()

coord <- st_as_sf(sf::st_centroid(merged_file$geometry))
coord<- as(coord,"Spatial")
coord_correlation <- coord
coord_correlation$Cases..mean.<- merged_file$Cases..mean.
coord_correlation$Urban<- merged_file$Urban
coord_correlation$maxtemp_ago_mean<- merged_file$maxtemp_ago_mean
coord_correlation$popden_mean<- merged_file$popden_mean
coord_correlation$maxtemp_jul_mean<- merged_file$maxtemp_jul_mean
coord_correlation$hurs_may_mean<- merged_file$hurs_may_mean
coord_correlation$fcover_mean<- merged_file$fcover_mean
coord_correlation$cmi_ago_mean<- merged_file$cmi_ago_mean
coord_correlation$Wall<- merged_file$Wall

##### calculation for cross correlation
I <- crossCorrelation(coord_correlation$Cases..mean., coord_correlation$Urban, coords=sp::coordinates(coord_correlation), clust = TRUE, k=99)
MORAN_res<- rbind(MORAN_res,data.frame(comb="cases-urban",moran=I$I))
coord_correlation$cases_urban_LISA <- I$SCI[,"lsci.xy"] ##local
coord_correlation$cases_urban_LISA_K <- as.factor(I$cluster)

I <- crossCorrelation(coord_correlation$Cases..mean., coord_correlation$maxtemp_ago_mean, coords=sp::coordinates(coord_correlation), clust = TRUE, k=99)
MORAN_res<- rbind(MORAN_res,data.frame(comb="cases-maxaug",moran=I$I))
coord_correlation$cases_maxaug_LISA <- I$SCI[,"lsci.xy"] ##local
coord_correlation$cases_maxaug_LISA_K <- as.factor(I$cluster)

I <- crossCorrelation(coord_correlation$Cases..mean., coord_correlation$popden_mean, coords=sp::coordinates(coord_correlation), clust = TRUE, k=99)
MORAN_res<- rbind(MORAN_res,data.frame(comb="cases-pop",moran=I$I))
coord_correlation$cases_pop_LISA <- I$SCI[,"lsci.xy"] ##local
coord_correlation$cases_pop_LISA_K <- as.factor(I$cluster)

I <- crossCorrelation(coord_correlation$Cases..mean., coord_correlation$maxtemp_jul_mean, coords=sp::coordinates(coord_correlation), clust = TRUE, k=99)
MORAN_res<- rbind(MORAN_res,data.frame(comb="cases-maxjul",moran=I$I))
coord_correlation$cases_maxjul_LISA <- I$SCI[,"lsci.xy"] ##local
coord_correlation$cases_maxjul_LISA_K <- as.factor(I$cluster)

I <- crossCorrelation(coord_correlation$Cases..mean., coord_correlation$hurs_may_mean, coords=sp::coordinates(coord_correlation), clust = TRUE, k=99)
MORAN_res<- rbind(MORAN_res,data.frame(comb="cases-hursmay",moran=I$I))
coord_correlation$cases_hursmay_LISA <- I$SCI[,"lsci.xy"] ##local
coord_correlation$cases_hursmay_LISA_K <- as.factor(I$cluster)

I <- crossCorrelation(coord_correlation$Cases..mean., coord_correlation$fcover_mean, coords=sp::coordinates(coord_correlation), clust = TRUE, k=99)
MORAN_res<- rbind(MORAN_res,data.frame(comb="cases-fcover",moran=I$I))
coord_correlation$cases_fcover_LISA <- I$SCI[,"lsci.xy"] ##local
coord_correlation$cases_fcover_LISA_K <- as.factor(I$cluster)

I <- crossCorrelation(coord_correlation$Cases..mean., coord_correlation$cmi_ago_mean, coords=sp::coordinates(coord_correlation), clust = TRUE, k=99)
MORAN_res<- rbind(MORAN_res,data.frame(comb="cases-cmiago",moran=I$I))
coord_correlation$cases_cmiago_LISA <- I$SCI[,"lsci.xy"] ##local
coord_correlation$cases_cmiago_LISA_K <- as.factor(I$cluster)

I <- crossCorrelation(coord_correlation$Cases..mean., coord_correlation$Wall, coords=sp::coordinates(coord_correlation), clust = TRUE, k=99)
MORAN_res<- rbind(MORAN_res,data.frame(comb="cases-Wall",moran=I$I))
coord_correlation$cases_wall_LISA <- I$SCI[,"lsci.xy"] ##local
coord_correlation$cases_wall_LISA_K <- as.factor(I$cluster)


shape_file_LISA<- shape_file
spatialcorr<- data.frame(ID_3=shape_file$ID_3,
                         cases_urban_LISA=coord_correlation$cases_urban_LISA,
                         cases_maxaug_LISA=coord_correlation$cases_maxaug_LISA,
                         cases_pop_LISA=coord_correlation$cases_pop_LISA,
                         cases_maxjul_LISA=coord_correlation$cases_maxjul_LISA,
                         cases_hursmay_LISA=coord_correlation$cases_hursmay_LISA,
                         cases_fcover_LISA=coord_correlation$cases_fcover_LISA,
                         cases_cmiago_LISA=coord_correlation$cases_cmiago_LISA,
                         cases_wall_LISA=coord_correlation$cases_wall_LISA,
                         
                         cases_urban_LISA_K=coord_correlation$cases_urban_LISA_K,
                         cases_maxaug_LISA_K=coord_correlation$cases_maxaug_LISA_K,
                         cases_pop_LISA_K=coord_correlation$cases_pop_LISA_K,
                         cases_maxjul_LISA_K=coord_correlation$cases_maxjul_LISA_K,
                         cases_hursmay_LISA_K=coord_correlation$cases_hursmay_LISA_K,
                         cases_fcover_LISA_K=coord_correlation$cases_fcover_LISA_K,
                         cases_cmiago_LISA_K=coord_correlation$cases_cmiago_LISA_K,
                         cases_wall_LISA_K=coord_correlation$cases_wall_LISA_K)
shape_file_LISA<- left_join(shape_file_LISA, spatialcorr, by=c("ID_3"="ID_3"))

gs.palb <- colorRampPalette(c("#6495ED","#e8e8e8"),bias=.1,space="rgb")
gs.palp <- colorRampPalette(c("#ed376d","#e8e8e8"),bias=.1,space="rgb")

##### plot cross correlation result
ggplot_URBAN <- ggplot() + theme_void() + theme(legend.position = c(0.1,0.25), plot.title = element_text(hjust = 0.5)) +
                geom_sf(data=shape_file_LISA, aes(fill=cases_urban_LISA_K), color= "#E8E8E8", size=.01)+
                geom_sf(data=shape_file, fill = "transparent", colour = "black", size = 2) +
                labs(fill = "") +
                ggtitle(paste("CIR - Urban")) +
                scale_fill_manual(values=gs.palb(3))

ggplot_MAXAUG <- ggplot() + theme_void() + theme(legend.position = c(0.1,0.25), plot.title = element_text(hjust = 0.5)) +
                geom_sf(data=shape_file_LISA, aes(fill=cases_maxaug_LISA_K), color= "#E8E8E8", size=.01)+
                geom_sf(data=shape_file, fill = "transparent", colour = "black", size = 2) +
                labs(fill = "") +
                ggtitle(paste("CIR - Max Temp August")) +
                scale_fill_manual(values=gs.palb(3))

ggplot_POP <- ggplot() + theme_void() + theme(legend.position = c(0.1,0.25), plot.title = element_text(hjust = 0.5)) +
              geom_sf(data=shape_file_LISA, aes(fill=cases_pop_LISA_K), color= "#E8E8E8", size=.01)+
              geom_sf(data=shape_file, fill = "transparent", colour = "black", size = 2) +
              labs(fill = "") +
              ggtitle(paste("CIR - Population Density")) +
              scale_fill_manual(values=gs.palb(3))

ggplot_MAXJUL <- ggplot() + theme_void() + theme(legend.position = c(0.1,0.25), plot.title = element_text(hjust = 0.5)) +
                geom_sf(data=shape_file_LISA, aes(fill=cases_maxjul_LISA_K), color= "#E8E8E8", size=.01)+
                geom_sf(data=shape_file, fill = "transparent", colour = "black", size = 2) +
                labs(fill = "") +
                ggtitle(paste("CIR - Max Temp July")) +
                scale_fill_manual(values=gs.palb(3))

ggplot_HUMMAY <- ggplot() + theme_void() + theme(legend.position = c(0.1,0.25), plot.title = element_text(hjust = 0.5)) +
                geom_sf(data=shape_file_LISA, aes(fill=cases_hursmay_LISA_K), color= "#E8E8E8", size=.01)+
                geom_sf(data=shape_file, fill = "transparent", colour = "black", size = 2) +
                labs(fill = "") +
                ggtitle(paste("CIR - Mean Humidity May")) +
                scale_fill_manual(values=gs.palp(4))

ggplot_FCOVER <- ggplot() + theme_void() + theme(legend.position = c(0.1,0.25), plot.title = element_text(hjust = 0.5)) +
                geom_sf(data=shape_file_LISA, aes(fill=cases_fcover_LISA_K), color= "#E8E8E8", size=.01)+
                geom_sf(data=shape_file, fill = "transparent", colour = "black", size = 2) +
                labs(fill = "") +
                ggtitle(paste("CIR - Forest Cover")) +
                scale_fill_manual(values=gs.palb(3))

ggplot_MOISTAUG <- ggplot() + theme_void() + theme(legend.position = c(0.1,0.25), plot.title = element_text(hjust = 0.5)) +
                  geom_sf(data=shape_file_LISA, aes(fill=cases_cmiago_LISA_K), color= "#E8E8E8", size=.01)+
                  geom_sf(data=shape_file, fill = "transparent", colour = "black", size = 2) +
                  labs(fill = "") +
                  ggtitle(paste("CIR - Mean Moisture August")) +
                  scale_fill_manual(values=gs.palp(4))

ggplot_WIND <- ggplot() + theme_void() + theme(legend.position = c(0.1,0.25), plot.title = element_text(hjust = 0.5)) +
              geom_sf(data=shape_file_LISA, aes(fill=cases_wall_LISA_K), color= "#E8E8E8", size=.01)+
              geom_sf(data=shape_file, fill = "transparent", colour = "black", size = 2) +
              labs(fill = "") +
              ggtitle(paste("CIR - Mean wind")) +
              scale_fill_manual(values=gs.palb(3))

