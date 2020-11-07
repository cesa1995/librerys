# An R script to estimate MQ gas sensors correlation curve and compute Ro, min and max Rs/Ro
#
# Copyright (c) Davide Gironi, 2016
#
# Released under GPLv3.
# Please refer to LICENSE file for licensing information.

# How to use this script:
# 1) set limits as datasheet curve ("xlim" and "ylim")
#    ex.
#      xlim = c(10, 1000)
#      ylim = c(0.1, 10)
# 2) find out datasheet curve points, and write it out (to "pointsdata")
#    each line it's a point on cartesian coordinate system
#    the useful WebPlotDigitizer app can help you extract points from the graph
#    ex.
#      pointsdata = "
#        10.052112405371744, 2.283698378106183
#        20.171602728600178, 1.8052797165878915
#        30.099224396434586, 1.5715748803154423
#        50.09267987761949, 1.3195287228519417
#        80.38812026903305, 1.1281218760133969
#        90.12665922665023, 1.0815121769656304
#        100.52112405371739, 1.0430967861855598
#        199.62996638292853, 0.8000946404902397
#      "
# 3) optional for Ro estimation: measure the sensor resistance (set it to "mres" ohm value) at a know amount of gas
#    set it to 0 if you do not need the Ro estimation
#    ex.
#      mres = 26954
# 4) optional for Ro estimation: set the know amount of gas for the resistance measure of the previous step (to "mppm")
#    set it to 0 if you do not need the Ro estimation
#    ex.
#      mppm = 392
# 5) optional for min-max Rs/Ro estimation: set the minand max amount of gas the sensor will react to (as "minppm" and "maxppm")
#    set it to 0 if you do not need the min-max Rs/Ro estimation
#    ex.
#      minppm = 10
#      maxppm = 200
library(data.table)

#remove old variables
rm(list=ls())

#set input values
xlim = c(10, 1000)
ylim = c(0.1, 10)
minppm = 10
maxppm = 2000
mres = 26954
mppm = 392
pointsdata = "
10,2.286565641924796
14.747453906689044,2.004660422642884
20.057865279983282,1.8077686769634338
30.062952762213477,1.5700676570541772
39.90749625364179,1.4158603118397255
50.057915311791895,1.325711365590109
59.81380161769859,1.253023406427101
69.75658919429256,1.1955039638327578
80.0459183516984,1.1299551984243776
89.64959468322745,1.0780851432948584
100.40549215207336,1.0481131341546852
112.45185089705903,0.9999999999999994
199.7686067770452,0.8056055490172105
"

#load points using fread
setnames(points <- fread(pointsdata, sep=",", sep2="\n"), c("x","y"))

#set named list of points, and swapped list of points
#points will be used to plot and compute values as datasheet figure
#pointsrev will be used to plot and compute values for the correlation function, it's the datasheet figure with swapped axis
x <- as.vector(points[,x])
y <- as.vector(points[,y])
points = list(x=x, y=y)
pointsrev = list(x=y, y=x)

#the nls (Nonlinear Least Squares) it's used to perform the power regression on points
#in order to work, nls needs an estimation of staring values
#we use log-log slope estimation to find intitial values

#estimate fit curve initial values
xfirst = head(points$x, n=1)
xlast = tail(points$x, n=1)
yfirst = head(points$y, n=1)
ylast = tail(points$y, n=1)
bstart= log(ylast/yfirst)/log(xlast/xfirst)
astart = yfirst/(xfirst^bstart)
#perform the fit
fit <- nls("y~a*x^b", start=list(a=astart,b=bstart), data=points)

#estimate fitref curve initial values
xfirstrev = head(pointsrev$x, n=1)
xlastrev = tail(pointsrev$x, n=1)
yfirstrev = head(pointsrev$y, n=1)
ylastrev = tail(pointsrev$y, n=1)
bstartrev = log(ylastrev/yfirstrev)/log(xlastrev/xfirstrev)
astartrev = yfirstrev/(xfirstrev^bstartrev)
fitrev <- nls("y~a*x^b", start=list(a=astartrev,b=bstartrev), data=pointsrev)

#plot fit curve (log-log scale)
fiteq = function(x){coef(fit)["a"]*x^(coef(fit)["b"])}
plot(points, log="xy", col="blue", xlab="ppm", ylab="Rs/Ro", xlim=xlim, ylim=ylim, panel.first=grid(equilogs=FALSE))
curve(fiteq, col="red", add=TRUE)

#plot fitrev curve (log-log scale)
fiteqrev = function(x){coef(fitrev)["a"]*x^(coef(fitrev)["b"])}
plot(pointsrev, log="xy", col="blue", xlab="Rs/Ro", ylab="ppm", xlim=ylim, ylim=xlim, panel.first=grid(equilogs=FALSE))
curve(fiteqrev, col="red", add=TRUE)

#plot fit curve (linear scale)
fiteq = function(x){coef(fit)["a"]*x^(coef(fit)["b"])}
plot(points, col="blue", xlab="ppm", ylab="Rs/Ro", panel.first=grid(equilogs=FALSE))
curve(fiteq, col="red", add=TRUE)

#plot fitrev curve (linear scale)
fiteqrev = function(x){coef(fitrev)["a"]*x^(coef(fitrev)["b"])}
plot(pointsrev, col="blue", xlab="Rs/Ro", ylab="ppm", panel.first=grid(equilogs=FALSE))
curve(fiteqrev, col="red", add=TRUE)

#estimate min Rs/Ro
cat("\nCorrelation function coefficients")
cat("\nEstimated a\n")
cat("  ")
cat(coef(fitrev)["a"])
cat("\nEstimated b\n")
cat("  ")
cat(coef(fitrev)["b"])
cat("\n")

#estimate min Rs/Ro
if (minppm != 0) {
    minRsRo = (maxppm/coef(fitrev)["a"])^(1/coef(fitrev)["b"])
    cat("\nEstimated min Rs/Ro\n")
    cat("  ")
    cat(minRsRo)
    cat("\n")
}

#estimate max Rs/Ro
if (maxppm != 0) {
    maxRsRo = (minppm/coef(fitrev)["a"])^(1/coef(fitrev)["b"])
    cat("\nEstimated max Rs/Ro\n")
    cat("  ")
    cat(maxRsRo)
    cat("\n")
}

#estimate Ro
if (mppm != 0 && mres != 0) {
    Ro = mres*(coef(fitrev)["a"]/mppm)^(1/coef(fitrev)["b"])
    cat("\nEstimated Ro\n")
    cat("  ")
    cat(Ro)
    cat("\n")
}

# MQ-135 output 
# Correlation function coefficients
#Estimated a
#  110.9327
#Estimated b
#  -2.762331

#Estimated min Rs/Ro
#  0.3510122

#Estimated max Rs/Ro
#  2.3896

#Estimated Ro
#  42568.51
              