library(ncdf4)
library(zoo)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(tidyr)
library(nlme)
library(pracma)
library(ggplot2)
library(MARSS)
library(car)
library(FactoMineR)
library(ggpubr)
library(mgcv)

# using monthly NCEP/NCAR!!

# latest SLP data downloaded from http://apdrc.soest.hawaii.edu/las/v6/constrain?var=16882
# this version (through May 2019) is also in the toyGOA folder on the team drive

nc.slp <- nc_open("/Users/MikeLitzow/Documents/R/FATE2/toyGOA/36BF63C72F8E7093345001B0869D4DF4_ferret_listing.nc")

# now process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "TIME")  # seconds since 1-1-1970

d <- dates(raw, origin = c(1,1,0001))
yr <- years(d)

x <- ncvar_get(nc.slp, "LON53_101")
y <- ncvar_get(nc.slp, "LAT45_69")

SLP <- ncvar_get(nc.slp, "SLP", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SLP <- aperm(SLP, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check
z <- colMeans(SLP)   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)

# load pdo
download.file("http://jisao.washington.edu/pdo/PDO.latest", "~pdo")
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# load npgo
download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

# load SST
# load ERSSTv5 data
download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1948-01-01):1:(2019-05-01)][(0.0):1:(0.0)][(30):1:(66)][(150):1:(250)]", "~updated.sst")

nc <- nc_open("~updated.sst")

# get lat/long
x.t <- ncvar_get(nc, "longitude")
y.t <- ncvar_get(nc, "latitude")
lat.t <- rep(y.t, length(x.t))   # Vector of latitudes
lon.t <- rep(x.t, each = length(y.t))   # Vector of longitudes

# assign dates
raw <- ncvar_get(nc, "time") # seconds since January 1, 1970
h <- raw/(24*60*60)
d.t <- dates(h, origin = c(1,1,1970))

# year for processing later
m <- as.numeric(months(d.t))
yr <- as.numeric(as.character(years(d.t)))

# get required sst data
SST <- ncvar_get(nc, "sst")

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix

# plot to check
z <- colMeans(SST)   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image(x.t,y.t,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x.t,y.t,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

# set names
dimnames(SST) <- list(as.character(d.t), paste("N", lat.t, "E", lon.t, sep=""))

