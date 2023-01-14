
setwd('C:/workspace/R/dp/dp_program')


library(sp) # for data meuse
library(gstat) # for function variogram
library(spatstat) # for runifpoint, ppp, swedishpines
library(gatepoints) # for free hand selection (fhs)
library(geoR) # for simulation function grf (Gaussian random field)
library(pracma) # for function |dot|
library(Matrix) # for function |rankMatrix|
library(assist) # for the dataset |climate|
library(tidyverse)
# Sys.setenv(PATH="%PATH%;C:/Rtools/gcc-4.6.3/bin;c:/Rtools/bin")
library(maps)
library(rgl)
# library(rworldmap)
library(raster) # for function brick
library(rgdal)
library(rasterVis)  
library(usdm)
library(quadmesh)
library(spacetime)
library(CompRandFld) # for the simulation of spatio-temporal data
library(fields) # for |image.plot| function
require(scatterplot3d)
# implementation part - Czech data
library(RCzechia) 
library(readxl)
library(tibble)
library(ggplot2)
library(sf)
library(Rce)
library(RCzechia)
library(reshape)
library(data.table)
library(mgcv) # for GAM
library(RColorBrewer)
library(gridExtra) # for arranging plots made by ggplot
library(plotly)
library(deldir)
library(dplyr)



##################################
# Implementation part
#################################

path_excel <- "DATA/Archive/import_stations_parameters_V1.20.xlsx"
table_stations <- read_excel(path_excel)

head(table_stations)
summary(table_stations)

colnames(table_stations) <- c("id", "id_domain", "station_name", "measure_begin",
                              "measure_end", "longitude", "latitude", "elevation")
table_stations <- tail(table_stations, -1)

options(digits=5)
table_stations$longitude <- as.double(table_stations$longitude)
table_stations$latitude <- as.double(table_stations$latitude)

# locations to add to the map
len_df_station <- nrow(table_stations)

points <- vector(mode="list", length=len_df_station)

for (i in 1:len_df_station) {
  x <- as.double(table_stations[i, "longitude"])
  y <- as.double(table_stations[i, "latitude"])
  points[[i]] <- st_point(c(x, y))
}


# reference: https://cran.r-project.org/web/packages/RCzechia/vignettes/vignette.html


# ggplot does not play nice with {raster} package; a data frame is required
relief <- vyskopis("rayshaded") %>%
  as("SpatialPixelsDataFrame") %>% # old style format - {sp}
  as_tibble()

# report results
ggplot() +
  geom_raster(data = relief, aes(x = x, y  = y, alpha = -raytraced), # relief
              fill = "gray30",  show.legend = F) + # no legend is necessary
  geom_sf(data = subset(RCzechia::reky(), Major == T), # major rivers
          color = "steelblue", alpha = .7) +
  geom_sf(data = df_geo,  color = "black", pch = 4) +  # X marks the spot!
  theme_bw() +
  theme(axis.title = element_blank())


res_vyskopis <-  vyskopis("rayshaded")
sp.df <- vyskopis("rayshaded") %>% as("SpatialPixelsDataFrame")

coords.cz <- coordinates(sp.df)
coords.cz <- as.data.frame(coords.cz)
coords.cz$z <- sp.df@data[["raytraced"]] # add ealtitude

# read the data
# average temperature
df_Temp <- readRDS("DATA/CHMU_Output/df_T_Hodnota_V1.32.rds")
df_Temp_n <- na.omit(df_Temp)
df_Temp_n[, 2]


# percitipation
df_T <- readRDS("DATA/CHMU_Output/df_SRA_Hodnota_V1.32.rds")


num_cols_T <- ncol(df_T)
mat_T <- as.matrix(df_T[, c(6: (num_cols_T-5))])
min_T <- min(mat_T, na.rm=TRUE)
max_T <- max(mat_T, na.rm=TRUE)

# check the last column
df_T[, ncol(df_T)]
nrow(df_T)

df_T_n <- na.omit(df_T)

indices_selected <- which(df_T_n$`Jméno stanice` %in% df_Temp_n$`Jméno stanice`)
df_T_n <- df_T_n[indices_selected,]


df_info_indices <- c('station_id', 'station_name', 'longitude', 'latitude', 'altitude')

coords.cz.re <- coords.cz[seq(1, nrow(coords.cz), 100), ]


#' Function makes a Raster object for ploting
#'
get_r_data <- function(df_pred){
  
  long_max <- max(coords.cz$x)
  long_min <- min(coords.cz$x)
  la_max <- max(coords.cz$y)
  la_min <- min(coords.cz$y)
  
  r_obj <- raster(xmn=long_min, xmx=long_max, ymn=la_min, ymx=la_max, ncol=155, nrow=112) # 300, 250
  r_data <- rasterize(x=df_pred[, 1:2], # lon-lat data
                      y=r_obj, # raster object
                      field=df_pred$value, # vals to fill raster with
                      fun=mean) # aggregate function
  return(r_data)
}

stations_T <- na.omit(df_T)[, 1:5]
colnames(stations_T) <- df_info_indices 


# https://www.datanovia.com/en/blog/the-a-z-of-rcolorbrewer-palette/
# color template
# myPalette <- colorRampPalette(rev(brewer.pal(n=9, "Spectral")))
# sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(-3, 22))

# for rain drop 
myPalette <- colorRampPalette(brewer.pal(n=9, "Blues"))
sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(0, 16))

#######################
# 
# implementation 
# Spline
#
#######################


#' Function calculates a GAM model with splines, by default thin-plate
#'
#' year: if seasonal is set to FALSE, choose the data from a whole year and make a one-day-head prediction
gam_model <- function(df, bs='tp', datetime_cz='6-20', seasonal=TRUE, year=2020) {
  
  # prepare the points for interpolation
  # select a subset to low down the resolution
  
  # reindex
  rownames(coords.cz.re) <- NULL
  
  
  # preprocessing
  
  # get the subset for a concrete season
  if (seasonal) {
    df_sub<- df[, grep(datetime_cz, names(df))]
  } else { # take data of one year
    df_sub<- df[, grep(year, names(df))]
  }
  
  # columns of information about the stations
  df_info <- df[, seq(1:5)]
  
  colnames(df_info) <- df_info_indices
  
  df_sub <- cbind(df_info, df_sub)
  
  df_sub_long <- melt(setDT(df_sub), id.vars = df_info_indices, variable.name = "year")
  
  
  # remove the na values
  df_sub_long <- na.omit(df_sub_long)
  
  
  df_sub_long$year <- as.character(df_sub_long$year)
  df_sub_long$year <- as.POSIXct(df_sub_long$year, format="%Y-%m-%d")
  df_sub_long$year <- as.numeric(df_sub_long$year)
  
  # summary(df_T_sub_long)
  # str(df_T_sub_long)
  
  
  # GAM model
  
  # document: s function: https://www.rdocumentation.org/packages/gam/versions/1.20.1/topics/s
  # document: ti function: https://www.rdocumentation.org/packages/tis/versions/1.39/topics/ti (compare with te - add to text)
  #            more details: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/te.html
  
  # model I
  # tensor product interactive term
  sub_model = gam(value ~ s(longitude, latitude, altitude, bs=bs)+ s(year) + ti(longitude, latitude, altitude, year, d=c(3, 1)),
                  data=df_sub_long, family=gaussian(link="identity"), method="REML")
  
  # 
  # thin - plate without interaction # TODO: rerun this for the four seasons
  # separate the altitude from the longitude and latitude
  # to make model I an extension of this model
  # sub_model = gam(value ~ s(longitude, latitude, bs=bs) + s(altitude, bs=bs) + s(year),
  #                   data=df_sub_long, family=gaussian(link="identity"), method="REML")
  
  return(sub_model)
}



# write.csv(stations_T, "stations_T.csv", row.names = FALSE)

#' Funtion makes prediction using input GAM model
#'
gam_prediction  <- function(datetime_cz, gam_model) {
  # interpolation
  # prepare data for prediction
  df_pred_T_sub <- coords.cz.re
  colnames(df_pred_T_sub) <- c('longitude', 'latitude', 'altitude')
  
  # predict the weather in 2021 for thr given day
  date_pred_sub <- as.POSIXct(paste('2021-', datetime_cz, sep=""), format="%Y-%m-%d")
  date_pred_sub <- as.numeric(date_pred_sub, format="%Y-%m-%d")
  
  df_pred_T_sub$year <- date_pred_sub
  
  # prediction
  res_pred_T_sub = predict(gam_model,
                           df_pred_T_sub, type = "response")
  df_pred_T_sub$value <- res_pred_T_sub
  
  return(df_pred_T_sub)
}




# Spring
gam_T_spring <- gam_model(df_T, datetime_cz='3-20')
df_pred_gam_spring <- gam_prediction('3-20', gam_T_spring)

summary(gam_T_spring)
# get the degree of freedom 
gam_T_spring_f <- family(gam_T_spring)
gam_T_spring_f$getTheta(trans=TRUE)


r_data_spring <- get_r_data(df_pred_gam_spring)

r_data_spring_pts <- rasterToPoints(r_data_spring, spatial = TRUE)
r_data_spring_df  <- data.frame(r_data_spring_pts)


# Summer
gam_T_summer <- gam_model(df_T, datetime_cz='6-20')
df_pred_gam_summer <- gam_prediction('6-20', gam_T_summer)

r_data_summer <- get_r_data(df_pred_gam_summer)


r_data_summer_pts <- rasterToPoints(r_data_summer, spatial = TRUE)
r_data_summer_df  <- data.frame(r_data_summer_pts)


# Autumn
gam_T_autumn <- gam_model(df_T, datetime_cz='9-22')
df_pred_gam_autumn <- gam_prediction('9-22', gam_T_autumn)

r_data_autumn <- get_r_data(df_pred_gam_autumn)

r_data_autumn_pts <- rasterToPoints(r_data_autumn, spatial = TRUE)
r_data_autumn_df  <- data.frame(r_data_autumn_pts)


# Winter
gam_T_winter <- gam_model(df_T, datetime_cz='12-21')
df_pred_gam_winter <- gam_prediction('12-21', gam_T_winter)

r_data_winter <- get_r_data(df_pred_gam_winter)

r_data_winter_pts <- rasterToPoints(r_data_winter, spatial = TRUE)
r_data_winter_df  <- data.frame(r_data_winter_pts)

# save the models
saveRDS(gam_T_spring, file = "models/gam_T_spring.rds")
saveRDS(gam_T_summer, file = "models/gam_T_summer.rds")
saveRDS(gam_T_autumn, file = "models/gam_T_sutumn.rds")
saveRDS(gam_T_winter, file = "models/gam_T_winter.rds")



sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(-0.4, 13))

g_spring <- ggplot() +
  geom_raster(data=r_data_spring_df_SR , aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("3-20-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")


sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(-0.4, 13))

g_summer <- ggplot() +
  geom_raster(data=r_data_summer_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) + 
  ggtitle("6-20-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")


sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(-0.4, 3))

g_autumn <- ggplot() +
  geom_raster(data=r_data_autumn_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) + 
  ggtitle("9-22-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")

sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(-0.4, 13))

g_winter <- ggplot() +
  geom_raster(data=r_data_winter_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("12-22-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")


plot(g_spring)
plot(g_summer)
plot(g_autumn)
plot(g_winter)

grid.arrange(g_spring, g_summer, g_autumn, g_winter, nrow = 2, ncol=2)


# plot the model info of gam_T_spring
par(mfrow=c(2, 2))
par(mar = c(2, 1.5, 2, 1))
plot(gam_T_spring, all.terms=TRUE)


# one year ahead prediction
gam_year_T <- gam_model(df_T_n, year=2020, seasonal=FALSE)
df_pred_year_T <- gam_prediction('01-01', gam_year_T)

r_data_year_T <- get_r_data(df_pred_year_T)

r_data_year_pts <- rasterToPoints(r_data_year_T, spatial = TRUE)
r_data_year_df  <- data.frame(r_data_year_pts)


sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(-0.4, 4))
ggplot() +
  geom_raster(data=r_data_year_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("1-1-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")




# seasonal again
# for comparing with the one-day head version

date_seasonal <- '01-01' # set the date, when the values should be predicted
# the year is hard coded, set to be 2021

gam_T_seasonal <- gam_model(df_T, datetime_cz=date_seasonal)
df_pred_gam_seasonal <- gam_prediction(date_seasonal, gam_T_seasonal)

r_data_seasonal_T <- get_r_data(df_pred_gam_seasonal)

r_data_seasonal_pts <- rasterToPoints(r_data_seasonal_T, spatial = TRUE)
r_data_seasonal_df  <- data.frame(r_data_seasonal_pts)

ggplot() +
  geom_raster(data=r_data_seasonal_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("1-1-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")

# plot the model info of gam_T_spring
par(mfrow=c(2, 2))
par(mar = c(3, 2, 2, 1))
plot(gam_T_seasonal, all.terms=TRUE)



############################
#
# Application :
# Kriging
#
##############################

# ---------------------------
# e.g. gstat
data("DE_RB_2005")


sp = cbind(x = c(0,0,1), y = c(0,1,1))
row.names(sp) = paste("point", 1:nrow(sp), sep="")

sp = SpatialPoints(sp)
time = as.POSIXct("2010-08-05")+3600*(10:13)
m = c(10,20,30) # means for each of the 3 point locations
mydata = rnorm(length(sp)*length(time),mean=rep(m, 4))

IDs = paste("ID",1:length(mydata))
mydata = data.frame(values = signif(mydata,3), ID=IDs)

mydata = data.frame(values = signif(mydata,3), ID=IDs)
stfdf = STFDF(sp, time, mydata)

# ----------------------
df <- df_T_n

# total number of stations
nrow(df_T)

# number of stations with observations
nrow(df_T_n)

###############################################################
# seasonal prediction

datetime_cz <- '-01-01'

df_sub<- df[, grep(datetime_cz, names(df))]
###############################################################

#################################################################
# one-day-ahead prediction
year <- '2020' # year, the observations in which are used for predictions

df_sub<- df[, grep(year, names(df))]

# pick every 5th row, starting from the first
df_sub <- df_sub %>% 
  select(seq(1, last_col(), by = 6))

#################################################################


# columns of locations about the stations
df_info <- df[, seq(1:5)]

df_info_indices

colnames(df_info) <- df_info_indices
nrow(df_info)

df_info_sp <- SpatialPoints(df_info[3:5])
df_info_sp

df_sub <- cbind(df_info, df_sub)
lst_times <- colnames(df_sub[, 6:ncol(df_sub)])
times <- as.POSIXct(lst_times) # the times points
length(times)

id_vars <- c("longitude", "latitude", "altitude")
df_sub_long <- melt(setDT(df_sub), id.vars = df_info_indices, variable.name = "year") # convert to long form
nrow(df_sub_long)


values <- df_sub_long$value
values_df <- data.frame(values)

stfdf_T <- STFDF(df_info_sp, times, values_df)

# reference: gstat.pdf
# calculate the variogram
sample_vgm <- variogramST(formula=values~1, data=stfdf_T)

# plot the empirical ST variogram

z <- as.numeric(sample_vgm$gamma)
z[1] <- 0
x <- as.numeric(sample_vgm$timelag)
y <- as.numeric(sample_vgm$dist)
y[1] <- 0

col <- colorRampPalette(c("blue", "white"))(20)[1 + round(19*(z - min(z))/diff(range(z)))]
dxyz <- deldir::deldir(x, y, z = z, suppressMsge = TRUE)
persp3d(dxyz, col = col, front = "lines", back = "lines", xlab="time", ylab="dist", zlab="gamma")

model_T <- vgmST("productSum",
                 space=vgm(39, "Sph", 343, 0),
                 time=vgm(36, "Exp", 3, 0),
                 k=15)


# seems to be better
# how to determine the initial values?

# ---------------------------------------------------------------------------

arr_kappa_test <- seq(0.1, 10, by=0.5)
# arr_kappa_test <- seq(0.005, 0.01, by=0.001)

min_MSE <- Inf
mse <- Inf

arr_mse <- list()

arr_kappa <- list()
arr_kappa2 <- list()
arr_kappa3 <- list()

for (kappa in arr_kappa_test) {
  for (kappa2 in arr_kappa_test) {
    for (kappa3 in arr_kappa_test) {
      
      model_T <- vgmST("sumMetric",
                       space = vgm(4.4, "Lin", 300, kappa),
                       time = vgm(2.2, "Lin", 365.25, kappa2),
                       joint = vgm(26, "Exp", 9000, kappa3),
                       stAni = 51.7)
      
      fitted_vgm <- fit.StVariogram(object=sample_vgm, model=model_T, fit.method=0)
      # attributes(fitted_vgm)
      # print(attr(x=fitted_vgm, which="MSE"))
      
      mse <- round(attr(x=fitted_vgm, which="MSE"), digits=6)
      
      min_MSE <- min(mse, min_MSE)
      
      
      arr_mse <- append(arr_mse, list(mse))
      
      arr_kappa <- append(arr_kappa, list(kappa))
      arr_kappa2 <- append(arr_kappa, list(kappa2))
      arr_kappa3 <- append(arr_kappa, list(kappa3))
    }
  }
}

min_MSE

index_min <- which(arr_mse==min_MSE)
index_min

arr_kappa[1]

# ---------------------------------------------------------------------------


model_T <- vgmST("sumMetric",
                 space = vgm(4.5, "Lin", 300, 0.05),
                 time = vgm(2.5, "Lin", 365.25, 0.05),
                 joint = vgm(25, "Exp", 9000, 0.05),
                 stAni = 51.7)

# prodSumModel <- vgmST("productSum",
#                       space=vgm(39, "Sph", 343, 0),
#                       time= vgm(36, "Exp", 3, 0),
#                       k=15)
# 
fitted_vgm <- fit.StVariogram(object=sample_vgm, model=model_T, fit.method=0, fit.kappa=TRUE)
attributes(fitted_vgm)



coords.cz.re <- coords.cz[seq(1, nrow(coords.cz), 100), ]

# the data on which interpolation is conducted
df_info_new_sp <- SpatialPoints(coords.cz.re)
df_info_new_sp

times_new <- rep(as.POSIXct(paste("2021", datetime_cz, sep="")), times=nrow(coords.cz.re))

stfdf_T_new = STI(df_info_new_sp, times_new, times_new)

gc()

krige_res <- krigeST(values~1, data=stfdf_T, computeVar=TRUE, newdata=stfdf_T_new, modelList=fitted_vgm)


kri_preds <- krige_res@data[["var1.pred"]]

pred_df_T <- coords.cz.re
pred_df_T$value <- kri_preds

r_data_T <- get_r_data(pred_df_T)

r_data_pts <- rasterToPoints(r_data_T, spatial = TRUE)
r_data_df  <- data.frame(r_data_pts)


sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(0.5, 2.7))

kri_preds <- krige_res_SR_spring@data[["var1.pred"]]

summary(kri_preds) # (0.9676, 1.4736)

pred_df_T <- coords.cz.re
pred_df_T$value <- kri_preds

r_data_T <- get_r_data(pred_df_T)

r_data_pts <- rasterToPoints(r_data_T, spatial = TRUE)
r_data_df  <- data.frame(r_data_pts)

g_kri_spring <- ggplot() +
  geom_raster(data=r_data_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("03-20-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")

plot(g_kri_spring_SR)


sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(5, 9))

kri_preds <- krige_res_SR_summer@data[["var1.pred"]]

summary(kri_preds) # (6.795, 7.974 )

pred_df_T <- coords.cz.re
pred_df_T$value <- kri_preds

r_data_T <- get_r_data(pred_df_T)

r_data_pts <- rasterToPoints(r_data_T, spatial = TRUE)
r_data_df  <- data.frame(r_data_pts)

g_kri_summer <- ggplot() +
  geom_raster(data=r_data_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("06-20-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")


plot(g_kri_summer)


sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(4, 7))

kri_preds <- krige_res_SR_autumn@data[["var1.pred"]]

summary(kri_preds) # (4.761   , 5.942  )

pred_df_T <- coords.cz.re
pred_df_T$value <- kri_preds

r_data_T <- get_r_data(pred_df_T)

r_data_pts <- rasterToPoints(r_data_T, spatial = TRUE)
r_data_df  <- data.frame(r_data_pts)

g_kri_autumn <- ggplot() +
  geom_raster(data=r_data_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("09-22-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")


plot(g_kri_autumn)


sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(0.5, 3))

kri_preds <- krige_res_SR_winter@data[["var1.pred"]]

summary(kri_preds) # (0.9401     , 1.7671   )

pred_df_T <- coords.cz.re
pred_df_T$value <- kri_preds

r_data_T <- get_r_data(pred_df_T)

r_data_pts <- rasterToPoints(r_data_T, spatial = TRUE)
r_data_df  <- data.frame(r_data_pts)

g_kri_winte <- ggplot() +
  geom_raster(data=r_data_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("12-22-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")


plot(g_kri_winter)


grid.arrange(g_kri_spring, g_kri_summer , g_kri_autumn , g_kri_winter ,  nrow=2, ncol=2)


# grid.arrange(g_kri_spring_SR  , g_kri_summer_SR , g_kri_autumn_SR , g_kri_winter ,  nrow=2, ncol=2)




sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(0.2, 1.7))

kri_preds <- krige_res@data[["var1.pred"]]

summary(kri_preds) # (0.372     , 1.149    )

pred_df_T <- coords.cz.re
pred_df_T$value <- kri_preds

r_data_T <- get_r_data(pred_df_T)

r_data_pts <- rasterToPoints(r_data_T, spatial = TRUE)
r_data_df  <- data.frame(r_data_pts)

g_kri_seasonal <- ggplot() +
  geom_raster(data=r_data_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("12-22-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")


plot(g_kri_seasonal)

# variances



myPalette2 <- colorRampPalette(brewer.pal(n=9, "BuPu"))

sc <- scale_fill_gradientn(colours = myPalette2(50), limits=c(4.9, 6.1))

kri_preds_var <- sqrt(krige_res@data[["var1.var"]])

summary(kri_preds_var) # (5.27         , 5.55     )

pred_df_T <- coords.cz.re
pred_df_T$value <- kri_preds_var

r_data_T <- get_r_data(pred_df_T)

r_data_pts <- rasterToPoints(r_data_T, spatial = TRUE)
r_data_df  <- data.frame(r_data_pts)

g_kri_seasonal_var <- ggplot() +
  geom_raster(data=r_data_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("standard error") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")


plot(g_kri_seasonal_var)



grid.arrange(g_kri_seasonal, g_kri_seasonal_var ,  nrow=1, ncol=2)


# one day ahead


sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(-0.7, 1))

kri_preds <- krige_res@data[["var1.pred"]]

summary(kri_preds) # (-0.60076     , 0.69236     )

pred_df_T <- coords.cz.re
pred_df_T$value <- kri_preds

r_data_T <- get_r_data(pred_df_T)

r_data_pts <- rasterToPoints(r_data_T, spatial = TRUE)
r_data_df  <- data.frame(r_data_pts)

g_kri_one_day_ahead <- ggplot() +
  geom_raster(data=r_data_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("12-22-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")


plot(g_kri_one_day_ahead)

# variances



myPalette2 <- colorRampPalette(brewer.pal(n=9, "BuPu"))

sc <- scale_fill_gradientn(colours = myPalette2(50), limits=c(1, 2.4))

kri_preds_var <- sqrt(krige_res@data[["var1.var"]])

summary(kri_preds_var) # (1.31      ,2.28  )

pred_df_T <- coords.cz.re
pred_df_T$value <- kri_preds_var

r_data_T <- get_r_data(pred_df_T)

r_data_pts <- rasterToPoints(r_data_T, spatial = TRUE)
r_data_df  <- data.frame(r_data_pts)

g_kri_one_day_ahead_var <- ggplot() +
  geom_raster(data=r_data_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("standard error") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")


plot(g_kri_one_day_ahead_var)



grid.arrange(g_kri_one_day_ahead, g_kri_one_day_ahead_var ,  nrow=1, ncol=2)








# plot variances
# krige_res_var <- krige_res_spring_sumMetric_new@data[["var1.var"]]
# krige_res_var <- krige_res_one_day_sumMetric@data[["var1.var"]]
krige_res_var <- krige_res_one_day_sumMetric_season_Jan@data[["var1.var"]]

# try square variance
n <- length(krige_res_var)

krige_res_var_sqr <- sqrt(krige_res_var)

min(krige_res_var_sqr)
max(krige_res_var_sqr)

# myPalette2 <- colorRampPalette(rev(brewer.pal(n=9, "Purples")))
myPalette2 <- colorRampPalette(brewer.pal(n=9, "BuPu"))
sc2 <- scale_fill_gradientn(colours = myPalette2(50), limits=c(5, 5.7))

pred_df_T_var <- coords.cz.re
pred_df_T_var$value <- krige_res_var_sqr

r_data_T <- get_r_data(pred_df_T_var)

r_data_pts <- rasterToPoints(r_data_T, spatial = TRUE)
r_data_df  <- data.frame(r_data_pts)


g_var <- ggplot() +
  geom_raster(data=r_data_df, aes(x=x, y=y, fill=layer)) +
  sc2 +
  geom_point(aes(stations_T$longitude, stations_T$latitude), size=1.4, color="#474747", pch=16) +
  ggtitle("standard error") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")

# plot(g_spring_var)

grid.arrange(g_one_day , g_var, nrow=1, ncol=2)

# compare two densities

krige_res1 <- krige_res_one_day_sumMetric@data[["var1.pred"]] 
krige_res2 <- krige_res_one_day_sumMetric_season_Jan@data[["var1.pred"]]


kri_res1_df <- data.frame(krige_res1)
colnames(kri_res1_df) <- "temperature"

kri_res2_df <- data.frame(krige_res2)
colnames(kri_res2_df) <- "temperature"


kri_df <- kri_res1_df
kri_df['2'] <- kri_res2_df

colnames(kri_df) <- c("one-day-ahead", "seasonal")

kri_df_melt <- melt(setDT(kri_df), id.var = c()) 

colnames(kri_df_melt) <- c("type", "temperature")

ggplot(kri_df_melt, aes(temperature, colour=type, fill=type)) +
  geom_density(alpha=0.1)