# this file contains the code for the theoretical chapters in the thesis

# check all the installed packages
ip = as.data.frame(installed.packages()[, c(1, 3:4)])
ip = ip[is.na(ip$Priority), 1:2, drop = FALSE]
ip

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


##########################################
# Spatial variogram
##########################################

# resource: http://gsp.humboldt.edu/OLM/R/04_01_Variograms.html
data(meuse)
meuse
# create SpatialPointsDataFrame
# resource: https://cengel.github.io/R-spatial/intro.html
coordinates(meuse) = ~ x + y
plot(
  meuse$x,
  meuse$y,
  xlab = "x",
  ylab = "y",
  pch = 21,
  bg = 20,
  cex = 0.8
)

vario.meuse <- variogram(zinc ~ x + y, data = meuse)
vario.meuse
plot(vario.meuse, cex = 0.8, pch = 21)

# fit the empirical model with theoretical model
model.vario.meuse <- vgm(model = "Sph", nugget = 0.01)
fit.vario.meuse <-
  fit.variogram(vario.meuse, model = model.vario.meuse)
fit.vario.meuse
plot(vario.meuse, model = fit.vario.meuse)


###########################################
# Point patterns generation
###########################################

# uniform distribution - random
pp <- runifpoint(n = 500, win = owin(c(0, 100), c(0, 100)))
par(mar = c(0.5, 0.5, 1, 0.5))

# gaussian distribution - clustered
x <- rnorm(200, mean = 40, sd = 10)
y <- rnorm(200, mean = 40, sd = 10)
X <- ppp(x, y, window = owin(c(0, 100), c(0, 100)))
X

# possion dis. along x, uniform along y
s.p <- seq(5, 95, 10)
x.poi <- rpois(200, lambda = 70)
y.poi <- runif(200, min = 0, max = 100)
X.poi <- ppp(x.poi, y.poi, window = owin(c(0, 100), c(0, 100)))
X.poi

# regular
s <- seq(from = 5, to = 95, by = 10)
x.reg <- rep(s, each = 10)
y.reg <- rep(s, times = 10)
X.reg <- ppp(x.reg, y.reg, window = owin(c(0, 100), c(0, 100)))
X.reg

par(mfrow = c(2, 2))
plot(
  pp,
  main = "(a) Random",
  pch = 20,
  cex = 1,
  lwd = 0.1
)
plot(
  X,
  main = "(b) Gaussian",
  pch = 20,
  cex = 1,
  lwd = 0.1
)
plot(
  X.poi,
  main = "(c) Poisson along x-axis",
  pch = 20,
  cex = 1,
  lwd = 0.1
)
plot(
  X.reg,
  main = "(d) Regular",
  pch = 20,
  cex = 1,
  lwd = 0.1
)

# Quadrat count
par(mfrow = c(1, 1))

data(swedishpines)
plot(swedishpines)

Q <- quadratcount(swedishpines, nx = 3, ny = 3)
plot(Q, add = TRUE, size = 15)

summary(swedishpines)
ei <- 71 / 9
ei
res <- sum((Q - ei) ^ 2 / ei)
# --> acording to the table of chi-square distribution, the p-value is between 0.1 and 0.9
# or directly use the function from theh package
res.q <- quadrat.test(swedishpines, 3, 3)
res.q
# --> p-value is 0.4169

# variance to mean, t-test
n.g <- 10
Q2 <- quadratcount(swedishpines, nx = n.g, ny = n.g)
plot(swedishpines, pch = 19)
plot(Q2, add = TRUE, size = 15)

n.T <- n.g * n.g # number of cells
n.bar <- 71 / (n.g * n.g)
s2 <- 1 / (n.T - 1) * sum((Q2 - n.bar) ^ 2)
ratio <- s2 / n.bar
t <- (ratio - 1) / sqrt(2 / (n.T - 1))
t
pt(t, df = n.T - 1, lower.tail = TRUE)

###########################
# Spatial Kriging
###########################

# set.seed(123)
set.seed(125)
sim.p <-  grf(
  7,
  cov.pars = c(1, .25),
  nug = 0.5,
  xlim = c(0, 10),
  ylim = c(0, 10),
  mean = 10
)
par(mar = c(5, 1, 2, 1))
points(sim.p,
       lwd = 5,
       xlab = "x",
       ylab = "y")
sim.p

coord.df <- as.data.frame(sim.p$coords)
distances <- dist(coord.df, method = "euclidean")
ma <- ceiling(max(distances))
mi <- floor(min(distances))

# automatic semivariogram calculation using |variog|
plot(sim.p)
# but this is the variogram of all the points
v.aut <-
  variog(
    sim.p,
    max.dist = 11,
    estimator.type = "classical",
    option = "bin",
    uvec = seq(1, ma, by = 1)
  )
plot(v.aut)
values <- sim.p$data

# find the coordinates of the selected data
plot(coord.df)
# target <- fhs(coord.df)

vario.df <- as.data.frame(seq(mi, ma))
names(vario.df)[1] <- "distance"

# manual calculation of the variogram (necessary to exclude the target point)
n <- nrow(coord.df) - 1
# Exclude point s7, which is the point we are going to predict
indices <- t(combn(1:n, m = 2))

# helper function
find.dist <- function(p.i1, p.i2) {
  # the points are rows in coord.df
  p1 <- coord.df[p.i1, ]
  p2 <- coord.df[p.i2, ]
  
  p1.x <- p1[1]
  p1.y <- p1[2]
  p2.x <- p2[1]
  p2.y <- p2[2]
  dp <- sqrt((p1.x - p2.x) ^ 2 + (p1.y - p2.y) ^ 2)
  
  return(dp)
}

i <- 1
for (d in vario.df$distance) {
  n.indices <- nrow(indices)
  n.pairs <- 0
  s.sum <- 0
  for (j in 1:n.indices) {
    # get the indices of the pair of points
    pair <- indices[j, ]
    p.i1 <- pair[1]
    p.i2 <- pair[2]
    
    dp <- find.dist(p.i1, p.i2)
    
    if (d - 0.5 <= dp & dp < d + 0.5) {
      n.pairs <- n.pairs + 1
      # # square sum
      z <- values[p.i1]
      zh <- values[p.i2]
      s.sum <- s.sum + (z - zh) ^ 2
    }
  }
  vario.df$num.pairs[i] <- n.pairs
  gamma.h <- 0.5 * s.sum / n.pairs
  # print(gamma.h)
  vario.df$vario.h[i] <- gamma.h
  i <- i + 1
}
plot(vario.df$vario.h,
     xlab = "index",
     ylab = "semivariogram",
     main = "Empirical semivariogram")

# helper function
find.vario <- function(p.i1, p.i2) {
  if (p.i1 == p.i2) {
    return(0)
  }
  n.d <- nrow(vario.df)
  dp <- find.dist(p.i1, p.i2)
  for (i in 1:n.d) {
    d <- i
    if (d - 0.5 <= dp & dp < d + 0.5) {
      v <- vario.df$vario.h[i]
    }
  }
  return(v)
}

# construct kriging matrix
M.k <- matrix(0, nrow = n + 1, ncol = n + 1)
for (i in 1:(n)) {
  for (j in 1:(n)) {
    if (i == j) {
      M.k[i, j] <- 0
    } else {
      # find semivariogram by searching in |vario.df|
      M.k[i, j] <- find.vario(i, j)
    }
  }
}
# assign values to the last row
for (j in 1:n) {
  M.k[n + 1, j] <- 1
}
# assign values to the last column
for (i in 1:n) {
  M.k[i, n + 1] <- 1
}

rankMatrix(M.k)
det(M.k)

# the vector or variogram betwee point s7 and other points
b.v <- matrix(0, nrow = n + 1, ncol = 1)
# the index of the target
target <- 7
for (i in 1:n) {
  b.v[i, 1] <- find.vario(i, target)
}
b.v[n + 1, 1] <- 1
x.v <- solve(M.k, b.v)
pred <- dot(x.v[1:n], values[1:n])
sigma <- dot(b.v, x.v)
sigma

# the kriging matrix for simple kriging
M.s <- M.k
# remove the last column
M.s <- M.s[, -(n + 1)]
# remove the last row
M.s <- M.s[-(n + 1), ]
b.s <- b.v
b.s <- b.s[-(n + 1)]
x.s <- b.s
as.matrix(x.s)

c0 <- var(values[1:n])
# weighting of the mean (k)
sum.lambda <- sum(x.s)
sum.lambda
# estimation of the mean
mu <- mean(values[1:n])
mu
pred.s <- dot(x.s, values[1:n]) + (1 - sum.lambda) * mu
pred.s
# variance of simple kriging
sigma.s <- (1 - sum.lambda) * c0 + dot(b.s, x.s)
sigma.s


###############################################
# Fitting using CompRandFld
#
# https://cran.r-project.org/web/packages/CompRandFld/CompRandFld.pdf
# now available here: https://mran.revolutionanalytics.com/snapshot/2019-12-03/web/packages/CompRandFld/CompRandFld.pdf
###############################################

set.seed(125)
sim.p <-  grf(
  7,
  cov.pars = c(1, .25),
  ylim = c(0, 10),
  xlim = c(0, 10),
  nug = 0.5,
  mean = 2
)
par(mar = c(5, 1, 2, 1))
points(sim.p,
       lwd = 5,
       xlab = "x",
       ylab = "y")
sim.p

coord.df <- as.data.frame(sim.p$coords)

# Empirical variogram
fit.emp <- EVariogram(sim.p$data, coord.df$x, coord.df$y)

corrmodel = "gauss"

mean <- 5
sill <- 1
nugget <- 0
scale <- 0.5

param <- list(
  mean = mean,
  sill = sill,
  nugget = nugget,
  scale = scale
)

start <- list(scale = scale, sill = sill)
fixed <- list(mean = mean, nugget = 0)

# theoretical model
fit.teo <-
  FitComposite(
    sim.p$data,
    coordx = sim.p$coords,
    corrmodel = corrmodel,
    likelihood = "Full",
    type = "Standard",
    start = start,
    fixed = fixed,
    maxdist = 5,
    taper = "Wendland1"
  )

print(fit)

xx <- seq(0, 100, 0.1)
loc_to_pred <- as.matrix(expand.grid(xx, xx))

pr <- Kri(
  loc = loc_to_pred,
  coordx = sim.p$coords,
  corrmodel = corrmodel,
  param = as.list(c(fit.teo$param, fit.teo$fixed)),
  data = sim.p$data
)

par(mfrow = c(1, 2))
par(mar = c(0.5, 4, 3, 3))
color <-
  colorRampPalette(
    c(
      "#6800ff",
      "blue",
      "#00e4ff" ,
      "#00ff93",
      "#68ff00",
      "#a6ff00",
      "yellow",
      "orange",
      "red"
    )
  )(100)
# simple kriging map prediction
image.plot(
  xx,
  xx,
  matrix(pr$pred, ncol = length(xx)),
  col = color,
  xlab = "",
  ylab = "",
  main = "Simple Kriging"
)

# simple kriging map prediction variance
image.plot(
  xx,
  xx,
  matrix(pr$varpred, ncol = length(xx)),
  col = color,
  xlab = "",
  ylab = "",
  main = "Std error"
)


pr <- Kri(
  loc = loc_to_pred,
  coordx = sim.p$coords,
  corrmodel = corrmodel,
  param = as.list(c(fit.teo$param, fit.teo$fixed)),
  data = sim.p$data,
  type_krig = "Ordinary"
)

# ordinary kriging map prediction
image.plot(
  xx,
  xx,
  matrix(pr$pred, ncol = length(xx)),
  col = color,
  xlab = "",
  ylab = "",
  main = "Ordinary Kriging"
)

# ordinary kriging map prediction variance
image.plot(
  xx,
  xx,
  matrix(pr$varpred, ncol = length(xx)),
  col = color,
  xlab = "",
  ylab = "",
  main = "Std error"
)


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
  geom_polygon(
    aes(x = long, y = lat, group = group),
    data = world_map,
    fill = "#939393",
    color = "#E1E1E1"
  ) +
  geom_point(aes(x = Y, y = X, color = "black"), cex = Z / 31) +
  scale_color_identity() +
  coord_fixed() +
  xlab("") +
  ylab("")


df.c <- data.frame(X, Y, Z)
# convert to spatial coordinates
coordinates(df.c) = ~ X + Y
vario.c <- variogram(Z ~ X + Y, data = df.c)
vario.c
plot(vario.c, cex = 0.8, pch = 21)

# fit the empirical model with theoretical model
par(mfrow = c(2, 2))
model.vario.c <- vgm(model = "Exp", nugget = 0.01)
fit.vario.c <- fit.variogram(vario.c, model = model.vario.c)
plot(vario.c, model = fit.vario.c, main = "Exponential model")

temp <- brick("tg_ens_mean_0.1deg_reg_v24.0e.nc")
plot(temp)

par(mar = c(2.5, 2.5, 1, 2))

temp50.winter <- temp$X1950.01.01
temp50.spring <- temp$X1950.03.20
temp50.summer <- temp$X1950.07.21
temp50.autumn <- temp$X1950.12.21

par(mfrow = c(2, 2))
plot(temp50.winter)
plot(temp50.spring)
plot(temp50.summer)
plot(temp50.autumn)


# use the layer of 1950 Jan.1 as example for spatial data
par(mfrow = c(1, 1))
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
x <- coords[, 1]
y <- coords[, 2]
# Define the temporal-coordinates:
instances <- 6
times <- seq(1, instances, 1)

set.seed(140)
# Simulation of a spatial Gaussian random field:
sim <- RFsim(
  x,
  y,
  times,
  corrmodel = "gneiting",
  grid = FALSE,
  param = list(
    mean = 0,
    scale_s = 0.4,
    scale_t = 1,
    sill = 1,
    nugget = 0,
    power_s = 1,
    power_t = 1,
    sep = 0.5
  )
)

# Spatial simulated data at first temporal instant
values <- sim$data
sizes <- values / apply(values, 1, max)

sim.df <- as.data.frame(coords)
# color palette
rbPal <- colorRampPalette(c('#ececec', 'black'))

par(mfrow = c(2, 3))
par(mar = c(5, 4, 3, 2))

for (i in 1:instances) {
  #This adds a column of color values
  # based on the y values
  colors <- rbPal(10)[as.numeric(cut(sizes[i, ], breaks = 10))]
  plot(
    coords,
    cex = sizes[i, ] * 5.5,
    col = colors,
    pch = 20,
    lwd = 2,
    xlab = "",
    ylab = ""
  )
  mtext(
    paste("t =", i),
    side = 1,
    line = 3,
    cex = 1.5,
    col = "black"
  )
}

fit <- EVariogram(sim$data, x, y, times, maxtime = 6, maxdist = 4)

# Results: Marginal spatial empirical variogram
par(mfrow = c(2, 2))
par(mar = c(6, 4, 3, 2))
plot(
  fit$centers,
  fit$variograms,
  xlab = 'h',
  ylab = expression(gamma(h)),
  ylim = c(0, max(fit$variograms)),
  xlim = c(0, max(fit$centers)),
  pch = 20,
  main = " ",
  cex.axis = .8,
  cex.lab = 0.9
)
mtext(
  "(a) Marginal spatial Variogram",
  side = 1,
  line = 3,
  cex = 1,
  col = "black"
)

plot(
  fit$bint,
  fit$variogramt,
  xlab = 't',
  ylab = expression(gamma(t)),
  ylim = c(0, max(fit$variograms)),
  xlim = c(0, max(fit$bint)),
  pch = 20,
  main = " ",
  cex.axis = .8,
  cex.lab = 0.9
)
mtext(
  "(b) Marginal temporal Variogram",
  side = 1,
  line = 3,
  cex = 1,
  col = "black"
)


# Building space-time variogram
st.vario <-
  matrix(fit$variogramst, length(fit$centers), length(fit$bint))
st.vario <-
  cbind(c(0, fit$variograms), rbind(fit$variogramt, st.vario))
# Results: 3d Spatio-temporal variogram

st.grid <- expand.grid(c(0, fit$centers), c(0, fit$bint))
scatterplot3d(
  st.grid[, 1],
  st.grid[, 2],
  c(st.vario),
  highlight.3d = TRUE,
  xlab = "h",
  ylab = "t",
  zlab = expression(gamma(h, t)),
  pch = 20,
  main = "",
  cex.axis = .7,
  mar = c(6, 3, 1.5, 2),
  mgp = c(0, 0, 0),
  cex.lab = .8,
  cex = 1.2
)
mtext(
  "(c) Space-time variogram",
  side = 1,
  line = 3,
  cex = 1,
  col = "black"
)
# A smoothed version
par(mar = c(5, 0.5, 1.2, 2), mgp = c(1, .3, 0))
persp(
  c(0, fit$centers),
  c(0, fit$bint),
  st.vario,
  xlab = "h",
  ylab = "u",
  zlab = expression(gamma(h, u)),
  ltheta = 90,
  shade = 0.75,
  ticktype = "detailed",
  phi = 30,
  theta = 30,
  main = "",
  cex.axis = .8,
  cex.lab = .9
)
mtext(
  "(d) Space-time variogram",
  side = 1,
  line = 3,
  cex = 1,
  col = "black"
)


fixed <- list(mean = 1, nugget = 0)
start <- list(scale_s = 0.5,
              scale_t = 1.2,
              sill = 1.5)
fit.teo <- FitComposite(
  sim$data,
  coordx = coords,
  coordt = times,
  corrmodel = "gneiting",
  likelihood = 'Marginal',
  type = 'Pairwise',
  start = start,
  fixed = fixed,
  maxdist = 0.5,
  maxtime = 3
)

par(mar = c(2, 4, 2, 2))
Covariogram(
  fit.teo,
  vario = fit,
  fix.lagt = 1,
  fix.lags = 1,
  show.vario = TRUE,
  pch = 20
)

################################################################
#### Spatio-temporal kriging
#### Using function Kri
################################################################

xx <- seq(0, 1, 0.01)
# xx<-seq(0,1,0.02)
loc_to_pred <- as.matrix(expand.grid(xx, xx))

times_to_pred <- 7:9

param <- as.list(c(fit.teo$param, fit.teo$fixed))

corrmodel <- "gneiting"

pr <-
  Kri(
    loc = loc_to_pred,
    time = times_to_pred,
    coordx = coords,
    coordt = times,
    corrmodel = corrmodel,
    param = param,
    data = sim$data,
    type_krig = "ordinary"
  )

par(mfrow = c(3, 2))

# visualization of the kriging result and standard error
color <-
  colorRampPalette(
    c(
      "#6800ff",
      "blue",
      "#00e4ff" ,
      "#00ff93",
      "#68ff00",
      "#a6ff00",
      "yellow",
      "orange",
      "red"
    )
  )(100)

par(mar = c(2, 1, 2, 1))

for (i in 1:3) {
  image.plot(
    xx,
    xx,
    matrix(pr$pred[i, ], ncol = length(xx)),
    col = color,
    main = paste("Kriging Time=" , i + 6),
    xlab = "",
    ylab = ""
  )
  image.plot(
    xx,
    xx,
    matrix(pr$varpred[i, ], ncol = length(xx)),
    col = color,
    main = paste("Std error Time=" , i + 6),
    xlab = "",
    ylab = ""
  )
}



## Simple Kriging
# set the mean
param$mean <- 0
pr <-
  Kri(
    loc = loc_to_pred,
    time = times_to_pred,
    coordx = coords,
    coordt = times,
    corrmodel = corrmodel,
    param = param,
    data = sim$data,
    type_krig = "simple"
  )

par(mfrow = c(3, 2))

# visualization of the kriging result and standard error
color <-
  colorRampPalette(
    c(
      "#6800ff",
      "blue",
      "#00e4ff" ,
      "#00ff93",
      "#68ff00",
      "#a6ff00",
      "yellow",
      "orange",
      "red"
    )
  )(100)

par(mar = c(2, 1, 2, 1))

for (i in 1:3) {
  image.plot(
    xx,
    xx,
    matrix(pr$pred[i, ], ncol = length(xx)),
    col = color,
    main = paste("Kriging Time=" , i + 6),
    xlab = "",
    ylab = ""
  )
  image.plot(
    xx,
    xx,
    matrix(pr$varpred[i, ], ncol = length(xx)),
    col = color,
    main = paste("Std error Time=" , i + 6),
    xlab = "",
    ylab = ""
  )
}
