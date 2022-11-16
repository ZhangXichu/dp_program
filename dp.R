# check all the installed packages
ip = as.data.frame(installed.packages()[,c(1,3:4)])
ip = ip[is.na(ip$Priority),1:2,drop=FALSE]
ip

setwd('C:/workspace/R/dp')


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

##########################################
# Spatial variogram
##########################################

# resource: http://gsp.humboldt.edu/OLM/R/04_01_Variograms.html
data(meuse)
meuse
# create SpatialPointsDataFrame
# resource: https://cengel.github.io/R-spatial/intro.html
coordinates(meuse) = ~x + y
plot(meuse$x, meuse$y, xlab="x", ylab="y",pch=21, bg=20, cex=0.8)

vario.meuse <- variogram(zinc ~ x + y, data=meuse)
vario.meuse
plot(vario.meuse, cex=0.8, pch=21)

# fit the empirical model with theoretical model
model.vario.meuse <- vgm(model="Sph", nugget=0.01)
fit.vario.meuse <- fit.variogram(vario.meuse, model=model.vario.meuse)
fit.vario.meuse
plot(vario.meuse, model=fit.vario.meuse)


###########################################
# Point patterns generation
###########################################

# uniform distribution - random
pp <- runifpoint(n=500, win=owin(c(0, 100), c(0, 100)))
par(mar=c(0.5, 0.5, 1, 0.5))
       
# gaussian distribution - clustered
x <- rnorm(200, mean=40, sd=10)
y <- rnorm(200, mean=40, sd=10)
X <- ppp(x, y, window=owin(c(0,100),c(0,100)))
X

# possion dis. along x, uniform along y
s.p <- seq(5, 95, 10)
x.poi <- rpois(200, lambda=70)
y.poi <- runif(200, min=0, max=100)
X.poi <- ppp(x.poi, y.poi, window=owin(c(0, 100), c(0, 100)))
X.poi

# regular
s <- seq(from=5, to=95, by = 10)
x.reg <- rep(s, each=10)
y.reg <- rep(s, times=10)
X.reg <- ppp(x.reg, y.reg, window=owin(c(0,100),c(0,100)))
X.reg

par(mfrow=c(2, 2))
plot(pp, main="(a) Random", pch=20, cex=1, lwd=0.1)  
plot(X, main="(b) Gaussian", pch=20, cex=1, lwd=0.1)
plot(X.poi, main="(c) Poisson along x-axis", pch=20, cex=1, lwd=0.1)
plot(X.reg, main="(d) Regular", pch=20, cex=1, lwd=0.1)

# Quadrat count
par(mfrow=c(1, 1))

data(swedishpines)
plot(swedishpines)

Q <- quadratcount(swedishpines, nx = 3, ny = 3)
plot(Q, add=TRUE, size=15)

summary(swedishpines)
ei <- 71 / 9
ei
res <- sum((Q - ei)^2 / ei)
# --> acording to the table of chi-square distribution, the p-value is between 0.1 and 0.9 
# or directly use the function from theh package
res.q <- quadrat.test(swedishpines, 3, 3)
res.q
# --> p-value is 0.4169

# variance to mean, t-test
n.g <- 10
Q2 <- quadratcount(swedishpines, nx=n.g, ny=n.g)
plot(swedishpines, pch=19)
plot(Q2, add=TRUE, size=15)

n.T <- n.g * n.g # number of cells
n.bar <- 71 / (n.g*n.g)
s2 <- 1/(n.T-1) * sum((Q2 - n.bar)^2)
ratio <- s2/n.bar
t <- (ratio-1) / sqrt(2/(n.T-1))
t
pt(t, df=n.T-1, lower.tail=TRUE)

###########################
# Spatial Kriging
###########################

# set.seed(123)
set.seed(125)
sim.p <-  grf(7, cov.pars = c(1, .25), nug=0.5,
              xlim=c(0, 10), ylim=c(0, 10), mean=10)
par(mar = c(5, 1, 2, 1))
points(sim.p, lwd=5, xlab="x", ylab="y")  
sim.p

coord.df <- as.data.frame(sim.p$coords)
distances <- dist(coord.df, method="euclidean")
ma <- ceiling(max(distances))
mi <- floor(min(distances))

# automatic semivariogram calculation using |variog|
plot(sim.p)
# but this is the variogram of all the points
v.aut <- variog(sim.p, max.dist=11, estimator.type="classical", option="bin", uvec=seq(1, ma, by=1))
plot(v.aut)
values <- sim.p$data

# find the coordinates of the selected data
plot(coord.df)
# target <- fhs(coord.df)

vario.df <- as.data.frame(seq(mi, ma))
names(vario.df)[1] <- "distance" 

# manual calculation of the variogram (necessary to exclude the target point)
n <- nrow(coord.df)-1
# Exclude point s7, which is the point we are going to predict
indices <- t(combn(1:n, m=2))

# helper function
find.dist <- function(p.i1, p.i2) {
  # the points are rows in coord.df
  p1 <- coord.df[p.i1,]
  p2 <- coord.df[p.i2,]
  
  p1.x <- p1[1]
  p1.y <- p1[2]
  p2.x <- p2[1]
  p2.y <- p2[2]
  dp <- sqrt((p1.x - p2.x)^2 + (p1.y - p2.y)^2)
  
  return(dp)
}

i <- 1
for (d in vario.df$distance) {
  n.indices <- nrow(indices)
  n.pairs <- 0
  s.sum <- 0
  for (j in 1:n.indices) {
    # get the indices of the pair of points
    pair <- indices[j,]
    p.i1 <- pair[1]
    p.i2 <- pair[2]

    dp <- find.dist(p.i1, p.i2)
      
    if (d-0.5 <= dp & dp < d+0.5) {
      n.pairs <- n.pairs + 1
      # # square sum
      z <- values[p.i1]
      zh <- values[p.i2]
      s.sum <- s.sum + (z - zh)^2
    }
  }
  vario.df$num.pairs[i] <- n.pairs
  gamma.h <- 0.5 * s.sum / n.pairs
  # print(gamma.h)
  vario.df$vario.h[i] <- gamma.h
  i <- i+1
}
plot(vario.df$vario.h, xlab="index", ylab="semivariogram", main="Empirical semivariogram")

# helper function
find.vario <- function(p.i1, p.i2) {
  if (p.i1 == p.i2) {
    return(0)
  }
  n.d <- nrow(vario.df)
  dp <- find.dist(p.i1, p.i2)
  for (i in 1:n.d) {
    d <- i
    if (d-0.5 <= dp & dp < d+0.5) {
      v <- vario.df$vario.h[i]
    }
  }
  return(v)
}

# construct kriging matrix
M.k <- matrix(0, nrow=n+1, ncol=n+1)
for (i in 1:(n)) {
  for (j in 1:(n)) {
    if (i == j) {
      M.k[i,j] <- 0 
    } else {
      # find semivariogram by searching in |vario.df|
      M.k[i,j] <- find.vario(i, j)
    }
  }
}
# assign values to the last row
for (j in 1:n) {
  M.k[n+1, j] <- 1
}
# assign values to the last column
for (i in 1:n) {
  M.k[i, n+1] <- 1
}

rankMatrix(M.k)
det(M.k)

# the vector or variogram betwee point s7 and other points 
b.v <- matrix(0, nrow=n+1, ncol=1)
# the index of the target 
target <- 7
for (i in 1:n) {
  b.v[i,1] <- find.vario(i,target)
}
b.v[n+1,1] <- 1
x.v <- solve(M.k, b.v)
pred <- dot(x.v[1:n], values[1:n])
sigma <- dot(b.v, x.v)
sigma

# the kriging matrix for simple kriging
M.s <- M.k
# remove the last column
M.s <- M.s[,-(n+1)]
# remove the last row
M.s <- M.s[-(n+1),]
b.s <- b.v
b.s <- b.s[-(n+1)]
x.s <- b.s
as.matrix(x.s)

c0 <- var(values[1:n])
# weighting of the mean (k)
sum.lambda <- sum(x.s)
sum.lambda
# estimation of the mean
mu <- mean(values[1:n])
mu
pred.s <- dot(x.s, values[1:n]) + (1-sum.lambda)*mu
pred.s
# variance of simple kriging
sigma.s <- (1-sum.lambda)*c0 + dot(b.s, x.s)
sigma.s


###############################################
# Fitting using CompRandFld
#
# https://cran.r-project.org/web/packages/CompRandFld/CompRandFld.pdf
# now available here: https://mran.revolutionanalytics.com/snapshot/2019-12-03/web/packages/CompRandFld/CompRandFld.pdf
###############################################

set.seed(125)
sim.p <-  grf(7, cov.pars = c(1, .25), ylim=c(0, 10), xlim=c(0, 10),
              nug=0.5, mean=2)
par(mar = c(5, 1, 2, 1))
points(sim.p, lwd=5, xlab="x", ylab="y")  
sim.p

coord.df <- as.data.frame(sim.p$coords)

# Empirical variogram
fit.emp <- EVariogram(sim.p$data, coord.df$x, coord.df$y)

corrmodel = "gauss"

mean <- 5
sill <- 1
nugget <- 0
scale <- 0.5

param <- list(mean=mean, sill=sill, nugget=nugget, scale=scale)

start <- list(scale=scale, sill=sill)
fixed <- list(mean=mean, nugget=0)

# theoretical model
fit.teo <- FitComposite(sim.p$data, coordx=sim.p$coords, corrmodel=corrmodel, 
                        likelihood="Full", type="Standard",
                        start=start, fixed=fixed, maxdist=5, taper="Wendland1")

print(fit)

xx<-seq(0, 10, 0.1)
loc_to_pred <- as.matrix(expand.grid(xx, xx))

pr<-Kri(loc=loc_to_pred, coordx=sim.p$coords, corrmodel=corrmodel,
        param=as.list(c(fit.teo$param,fit.teo$fixed)), data=sim.p$data)

par(mfrow=c(1,2))
par(mar = c(0.5, 4, 3, 3))
color <- colorRampPalette(c("#6800ff", "blue", "#00e4ff" , "#00ff93", "#68ff00", "#a6ff00", "yellow", "orange", "red"))(100)
# simple kriging map prediction
image.plot(xx, xx, matrix(pr$pred,ncol=length(xx)),col=color,
xlab="",ylab="",main="Simple Kriging")

# simple kriging map prediction variance
image.plot(xx, xx, matrix(pr$varpred,ncol=length(xx)),col=color,
xlab="",ylab="",main="Std error")


pr<-Kri(loc=loc_to_pred, coordx=sim.p$coords, corrmodel=corrmodel,
        param=as.list(c(fit.teo$param,fit.teo$fixed)), data=sim.p$data, type_krig="Ordinary")

# ordinary kriging map prediction
image.plot(xx, xx, matrix(pr$pred,ncol=length(xx)),col=color,
           xlab="",ylab="",main="Ordinary Kriging")

# ordinary kriging map prediction variance
image.plot(xx, xx, matrix(pr$varpred,ncol=length(xx)),col=color,
           xlab="",ylab="",main="Std error")


###############################
# data examples 
###############################

head(climate)
summary(climate)

X <- climate$lat.degree
Y <- climate$long.degree
Z <- climate$temp

world_map <- map_data("world")

ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = world_map, fill = "#939393", color = "#E1E1E1") +
  geom_point(aes(x =Y, y=X, color="black"), cex=Z/31) +
  scale_color_identity() +
  coord_fixed() +
  xlab("") +
  ylab("")


df.c <- data.frame(X, Y, Z)
# convert to spatial coordinates
coordinates(df.c) = ~X + Y
vario.c <- variogram(Z ~ X + Y, data=df.c)
vario.c
plot(vario.c, cex=0.8, pch=21)

# fit the empirical model with theoretical model
par(mfrow=c(2,2))
model.vario.c <- vgm(model="Exp", nugget=0.01)
fit.vario.c <- fit.variogram(vario.c, model=model.vario.c)
plot(vario.c, model=fit.vario.c, main="Exponential model")

temp <- brick("tg_ens_mean_0.1deg_reg_v24.0e.nc")
plot(temp)

par(mar=c(2.5, 2.5, 1, 2))

temp50.winter <- temp$X1950.01.01 
temp50.spring <- temp$X1950.03.20
temp50.summer <- temp$X1950.07.21
temp50.autumn <- temp$X1950.12.21

par(mfrow=c(2, 2))
plot(temp50.winter)
plot(temp50.spring)
plot(temp50.summer)
plot(temp50.autumn)


# use the layer of 1950 Jan.1 as example for spatial data 
par(mfrow=c(1,1))
plot(temp50)

# get mean
cellStats(temp50, 'mean')
# prepare necessary data from extracting value
x.min <- xmin(temp50)
x.max <- xmax(temp50)
y.min <- ymin(temp50)
y.max <- yman(temp50)

##################################################
### generate the coordinates
##################################################

set.seed(132)
sim1 <- grf(20, cov.pars = c(1, .25))
# a display of simulated locations and values
points(sim1) 

coords <- sim1$coords
################################################################
###
### Simulation of a Gaussian random field.
### with double exponential correlation.
### One spatio-temporal replication.
### Using function RFsim
###
### https://cran.r-project.org/web/packages/CompRandFld/CompRandFld.pdf
###############################################################

# Define the spatial-coordinates of the points:
x <- coords[,1]
y <- coords[,2]
# Define the temporal-coordinates:
instances <- 6
times <- seq(1, instances, 1)

set.seed(140)
# Simulation of a spatial Gaussian random field:
sim <- RFsim(x, y, times, corrmodel="gneiting", grid=FALSE,
             param=list(mean=0,scale_s=0.4,scale_t=1,sill=1,
                        nugget=0,power_s=1,power_t=1,sep=0.5))

# Spatial simulated data at first temporal instant
values <- sim$data
sizes <- values/apply(values,1,max)

sim.df <- as.data.frame(coords)
# color palette
rbPal <- colorRampPalette(c('#ececec','black'))

par(mfrow=c(2, 3))
par(mar = c(5, 4, 3, 2))

for (i in 1:instances) {
  #This adds a column of color values
  # based on the y values
  colors <- rbPal(10)[as.numeric(cut(sizes[i,],breaks = 10))]
  plot(coords, cex = sizes[i,]*5.5, col=colors, pch=20, lwd=2, xlab="", ylab="")
  mtext(paste("t =", i), side=1, line=3, cex=1.5, col="black") 
}

fit <- EVariogram(sim$data, x, y, times, maxtime=6, maxdist=4)

# Results: Marginal spatial empirical variogram
par(mfrow=c(2,2))
par(mar = c(6, 4, 3, 2))
plot(fit$centers, fit$variograms, xlab='h', ylab=expression(gamma(h)),
     ylim=c(0, max(fit$variograms)), xlim=c(0, max(fit$centers)),
     pch=20, main=" ",cex.axis=.8,
     cex.lab=0.9)
mtext("(a) Marginal spatial Variogram", side=1, line=3, cex=1, col="black") 

plot(fit$bint, fit$variogramt, xlab='t', ylab=expression(gamma(t)),
     ylim=c(0, max(fit$variograms)),xlim=c(0,max(fit$bint)),
     pch=20,main=" ",cex.axis=.8, 
     cex.lab=0.9)
mtext("(b) Marginal temporal Variogram", side=1, line=3, cex=1, col="black") 


# Building space-time variogram
st.vario <- matrix(fit$variogramst,length(fit$centers),length(fit$bint))
st.vario <- cbind(c(0,fit$variograms), rbind(fit$variogramt,st.vario))
# Results: 3d Spatio-temporal variogram

st.grid <- expand.grid(c(0,fit$centers),c(0,fit$bint))
scatterplot3d(st.grid[,1], st.grid[,2], c(st.vario),
              highlight.3d=TRUE, xlab="h",ylab="t",
              zlab=expression(gamma(h,t)), pch=20,
              main="",cex.axis=.7,
              mar=c(6, 3, 1.5, 2), mgp=c(0,0,0),
              cex.lab=.8,
              cex=1.2)
mtext("(c) Space-time variogram", side=1, line=3, cex=1, col="black")
# A smoothed version
par(mar=c(5, 0.5, 1.2, 2),mgp=c(1,.3, 0))
persp(c(0,fit$centers), c(0,fit$bint), st.vario,
      xlab="h", ylab="u", zlab=expression(gamma(h,u)),
      ltheta=90, shade=0.75, ticktype="detailed", phi=30,
      theta=30,main="",cex.axis=.8,
      cex.lab=.9)
mtext("(d) Space-time variogram", side=1, line=3, cex=1, col="black")


fixed <- list(mean=1,nugget=0)
start <- list(scale_s=0.5,scale_t=1.2,sill=1.5)
fit.teo <- FitComposite(sim$data, coordx=coords, coordt=times,
                    corrmodel="gneiting", likelihood='Marginal',
                    type='Pairwise',start=start,fixed=fixed,
                    maxdist=0.5,maxtime=3)

par(mar = c(2, 4, 2, 2))
Covariogram(fit.teo, vario=fit, fix.lagt=1, fix.lags=1, show.vario=TRUE, pch=20)

################################################################
#### Spatio-temporal kriging
#### Using function Kri
################################################################

xx<-seq(0,1,0.01)
# xx<-seq(0,1,0.02)
loc_to_pred<-as.matrix(expand.grid(xx,xx))

times_to_pred<-7:9

param<-as.list(c(fit.teo$param,fit.teo$fixed))

corrmodel<-"gneiting"

pr<-Kri(loc=loc_to_pred, time=times_to_pred, coordx=coords, coordt=times,
        corrmodel=corrmodel, param=param, data=sim$data, type_krig="ordinary")

par(mfrow=c(3, 2))

# visualization of the kriging result and standard error
color <- colorRampPalette(c("#6800ff", "blue", "#00e4ff" , "#00ff93", "#68ff00", "#a6ff00", "yellow", "orange", "red"))(100)

par(mar = c(2, 1, 2, 1))

for(i in 1:3){
  image.plot(xx, xx, matrix(pr$pred[i,],ncol=length(xx)),col=color,
             main = paste("Kriging Time=" , i+6), xlab="", ylab="")
  image.plot(xx, xx, matrix(pr$varpred[i,],ncol=length(xx)),col=color,
             main = paste("Std error Time=" , i+6), xlab="", ylab="")
}



## Simple Kriging
# set the mean
param$mean <- 0
pr<-Kri(loc=loc_to_pred, time=times_to_pred, coordx=coords, coordt=times,
        corrmodel=corrmodel, param=param, data=sim$data, type_krig="simple")

par(mfrow=c(3, 2))

# visualization of the kriging result and standard error
color <- colorRampPalette(c("#6800ff", "blue", "#00e4ff" , "#00ff93", "#68ff00", "#a6ff00", "yellow", "orange", "red"))(100)

par(mar = c(2, 1, 2, 1))

for(i in 1:3){
  image.plot(xx, xx, matrix(pr$pred[i,],ncol=length(xx)),col=color,
             main = paste("Kriging Time=" , i+6), xlab="", ylab="")
  image.plot(xx, xx, matrix(pr$varpred[i,],ncol=length(xx)),col=color,
             main = paste("Std error Time=" , i+6), xlab="", ylab="")
}


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
df_T <- readRDS("DATA/CHMU_Output/df_T_Hodnota_V1.32.rds")
nrow(df_T)
df_T_n <- na.omit(df_T)

# check the last column
df_T[, ncol(df_T)]

num_cols_T <- ncol(df_T)
mat_T <- as.matrix(df_T[, c(6: (num_cols_T-5))])
min_T <- min(mat_T, na.rm=TRUE)
max_T <- max(mat_T, na.rm=TRUE)

# lst_T <- unlist(df_T[, c(6: (num_cols_T-5))])
# # all the values from the average temperature table without Na
# lst_T <- lst_T[!is.na(lst_T)]
# 
# plot(density(lst_T))
# quantile(lst_T, probs=0.25)
# quantile(lst_T, probs=0.75)

# reindex
# (coords.cz.re) <- NULL


df_info_indices <- c('station_id', 'station_name', 'longitude', 'latitude', 'altitude')

coords.cz.re <- coords.cz[seq(1, nrow(coords.cz), 1), ]

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

stations_T <- na.omit(df_T)[, 1:5]
colnames(stations_T) <- df_info_indices 

write.csv(stations_T, "stations_T.csv", row.names = FALSE)

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


long_max <- max(coords.cz$x)
long_min <- min(coords.cz$x)
la_max <- max(coords.cz$y)
la_min <- min(coords.cz$y)

#' Function makes a Raster object for ploting
#'
get_r_data <- function(df_pred){
  r_obj <- raster(xmn=long_min, xmx=long_max, ymn=la_min, ymx=la_max, ncol=300, nrow=250)
  r_data <- rasterize(x=df_pred[, 1:2], # lon-lat data
                      y=r_obj, # raster object
                      field=df_pred$value, # vals to fill raster with
                      fun=mean) # aggregate function
  return(r_data)
}


# color template
myPalette <- colorRampPalette(rev(brewer.pal(n=9, "Spectral")))
sc <- scale_fill_gradientn(colours = myPalette(50), limits=c(-5, 22))

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

plot_T <- function() {
  g1 <- ggplot() +
    geom_raster(data=r_data_spring_df, aes(x=x, y=y, fill=layer)) +
    sc +
    geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
    ggtitle("3-20") +
    theme(plot.title=element_text(hjust = 0.5)) +
    xlab("longitude") + ylab("latitude") +
    labs(fill = "")
  
  g2 <- ggplot() +
    geom_raster(data=r_data_summer_df, aes(x=x, y=y, fill=layer)) +
    sc +
    geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) + 
    ggtitle("6-20") +
    theme(plot.title=element_text(hjust = 0.5)) +
    xlab("longitude") + ylab("latitude") +
    labs(fill = "")
  
  g3 <- ggplot() +
    geom_raster(data=r_data_autumn_df, aes(x=x, y=y, fill=layer)) +
    sc +
    geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) + 
    ggtitle("9-22") +
    theme(plot.title=element_text(hjust = 0.5)) +
    xlab("longitude") + ylab("latitude") +
    labs(fill = "")
  
  g4 <- ggplot() +
    geom_raster(data=r_data_winter_df, aes(x=x, y=y, fill=layer)) +
    sc +
    geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
    ggtitle("12-22") +
    theme(plot.title=element_text(hjust = 0.5)) +
    xlab("longitude") + ylab("latitude") +
    labs(fill = "")
  
  grid.arrange(g1, g2, g3, g4, nrow = 2, ncol=2)
}

plot_T()

# plot the model info of gam_T_spring
par(mfrow=c(2, 2))
par(mar = c(2, 1.5, 2, 1))
plot(gam_T_spring, all.terms=TRUE)

# one year ahead prediction
gam_year_T <- gam_model(df_T, year=2020, seasonal=FALSE)
df_pred_year_T <- gam_prediction('01-01', gam_year_T)

r_data_year_T <- get_r_data(df_pred_year_T)

r_data_year_pts <- rasterToPoints(r_data_year_T, spatial = TRUE)
r_data_year_df  <- data.frame(r_data_year_pts)

ggplot() +
  geom_raster(data=r_data_year_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("1-1-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")

saveRDS(gam_year_T, file = "models/gam_year_T.rds")
saveRDS(df_pred_year_T, file = "models/df_pred_year_T.rds")



# seasonal again
# for comparing with the one-day head version
gam_T_day <- gam_model(df_T, datetime_cz='01-01')
df_pred_gam_day <- gam_prediction('01-01', gam_T_day)



r_data_day_T <- get_r_data(df_pred_gam_day)

r_data_day_pts <- rasterToPoints(r_data_day_T, spatial = TRUE)
r_data_day_df  <- data.frame(r_data_day_pts)

ggplot() +
  geom_raster(data=r_data_day_df, aes(x=x, y=y, fill=layer)) +
  sc +
  geom_point(aes(stations_T$longitude, stations_T$latitude), , size=1.8, color="#575757", pch=19) +
  ggtitle("1-1-2021") +
  theme(plot.title=element_text(hjust = 0.5)) +
  xlab("longitude") + ylab("latitude") +
  labs(fill = "")

# plot the model info of gam_T_spring
par(mfrow=c(2, 2))
par(mar = c(3, 2, 2, 1))
plot(gam_T_day, all.terms=TRUE)

saveRDS(gam_T_day, file = "models/gam_T_day.rds")


# some analysis
df_T_na <- na.omit(df_T)
df_T_1 <- df_T_na[, grep('01-01', names(df_T_na))]


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


datetime_cz <- '3-20'
df <- df_T_n

# number of stations with observations
nrow(df_T_n)

seasonal <- TRUE

# get the subset for a concrete season
if (seasonal) {
  df_sub<- df[, grep(datetime_cz, names(df))]
} else { # take data of one year
  df_sub<- df[, grep(year, names(df))]
}

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

model_T <- vgmST("productSum",
                 space=vgm(39, "Sph", 343, 0),
                 time=vgm(36, "Exp", 3, 0),
                 k=15)

fitted_vgm <- fit.StVariogram(object=sample_vgm, model=model_T, fit.method=1)

# the data on which interpolation is conducted
df_info_new <- coords.cz.re




