install.packages("extRemes")
library(extRemes)
install.packages("in2extRemes")
library(in2extRemes)
install.packages("texmex")
library(texmex)
install.packages(VGAM)
library(VGAM)
install.packages("ismev")
library(ismev)
install.packages("evmix")
library(evmix)
install.packages("in2extRemes")
library(extRemes)
install.packages("evd")
library(evd)

ds <- read.csv(file = 'extremes.csv',header = TRUE,sep = ",")
x <- cbind(ds$mttc_neg,ds$deltav)
look <- as.data.frame(x)
look2 <- as.in2extRemesDataObject(look)
# BGEV
f1 <- fbvevd(x, model = c("log"), std.err = TRUE, method = 'SANN')
summary(f1)

loc1 = f1[["param"]][["loc1"]]
loc2 = f1[["param"]][["loc2"]]
sca1 = f1[["param"]][["scale1"]]
sca2 = f1[["param"]][["scale2"]]
shp1 = f1[["param"]][["shape1"]]
shp2 = f1[["param"]][["shape2"]]
dep = f1[["param"]][["dep"]]

se_loc1 = f1[["std.err"]][["loc1"]]
se_loc2 = f1[["std.err"]][["loc2"]]
se_sca1 = f1[["std.err"]][["scale1"]]
se_sca2 = f1[["std.err"]][["scale2"]]
se_shp1 = f1[["std.err"]][["shape1"]]
se_shp2 = f1[["std.err"]][["shape2"]]
dep = f1[["std.err"]][["dep"]]

# confidence interval
sc1 <- rnorm(100, sca1, (se_sca1*(sqrt(nrow(ds))))) # SE/sqrt of n # rnorm (numbers required, mean, sd)
lc1 <- rnorm(100, loc1, (se_loc1*(sqrt(nrow(ds)))))
sh1 <- rnorm(100, shp1, (se_shp1*(sqrt(nrow(ds)))))
sh1 = abs(sh1)
sc2 <- rnorm(100, sca2, (se_sca2*(sqrt(nrow(ds))))) # SE/sqrt of n # rnorm (numbers required, mean, sd)
lc2 <- rnorm(100, loc2, (se_loc2*(sqrt(nrow(ds)))))
sh2 <- rnorm(100, shp2, (se_shp2*(sqrt(nrow(ds)))))
sh2 = abs(sh2)


T = 365*12*5 #daylight hours in one year 
t = 24 #conflict observation period
N <- 0

for (i in 1:100) {
  a1 = 1/sh1[i]
  b1 = sh1[i]*(lc1[i] /sc1[i])
  c1 = 1-b1
  d1 = c1^a1
  e11 = -1/dep
  var1 = d1^e11
  a2 = 1/sh2[i]
  b2 = sh2[i]*(lc2[i] /sc2[i])
  c2 = 1-b2
  d2 = c2^a2
  var2 = d2^e11
  var3 = (var1+var2)^dep
  var4 = exp(-var3)
  R = 1-var4
  N[i] = (T/t)*R
  print(i)
}
x2 <- na.omit(N)
quantile(x2, probs = c(0.025, 0.5, 0.975))

