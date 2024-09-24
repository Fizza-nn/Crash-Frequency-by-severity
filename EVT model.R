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

ds <- read.csv(file = 'Block maxima extremes.csv',header = TRUE,sep = ",")
x <- cbind(ds$mttc_neg,ds$deltav)
look <- as.data.frame(x)
look2 <- as.in2extRemesDataObject(look)
# BGEV
#ds <- read.csv(file = 'Logan_TTC_samp_Time2.csv',header = TRUE,sep = ',') # For TTC
#x <- cbind(ds$r_TTC,ds$DRAC)
f1 <- fbvevd(x, model = c("log"), std.err = TRUE, method = 'SANN')
#summary(f1)
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
  #e1 = exp(-d1)
  e11 = -1/dep
  var1 = d1^e11
  a2 = 1/sh2[i]
  b2 = sh2[i]*(lc2[i] /sc2[i])
  c2 = 1-b2
  d2 = c2^a2
  var2 = d2^e11
  #e2 = exp(-d2)
  var3 = (var1+var2)^dep
  var4 = exp(-var3)
  #R1 = e1^e11; R1 = e2^e11; R3 = (R1 + R2)^dep
  #R4 = -1*R3; R5 = exp(R4)
  #R = 1- R5
  R = 1-var4
  N[i] = (T/t)*R
  print(i)
}
x2 <- na.omit(N)
quantile(x2, probs = c(0.025, 0.5, 0.975))




# BGPD
# Fixed quantile rule, where a higher upper quantile (e.g. 5%) is used as the threshold
ds1 <- read.csv(file = 'Book2.csv',header = TRUE,sep = ',') # For TTC
u <- apply(x, 2, quantile, probs = 0.95)# automatic threshold selection

x <- cbind(ds1$TTC,ds1$DRAC)

bvtcplot(x, spectral = FALSE) #bivariate threshold choice plot

fbvpot(x, threshold=c(-0.504,5.262), model = c("log"), likelihood =
         c("censored"), method = 'Nelder-Mead')

abvevd(x = 0.5, 0.75, model = c("log"),
       rev = FALSE, plot = TRUE, col = 1,
       xlim = c(0,1), ylim = c(0.5,1), xlab = "t",
       ylab = "A(t)")


# Confidence interval
sca1 = 0.3956; se_sca1 =  0.3 #0.051
shp1 = -0.7918; se_shp1 = 0.7 #0.017
sca2 = 3.552; se_sca2 = 0.285
shp2 = 0.996; se_shp2 = 1.049
dep = 0.693

zz = nrow(ds1[(ds1[,1]>-0.504),])
zzz = nrow(ds1[(ds1[,2]>5.262),])

sc1 <- rnorm(100, sca1, (se_sca1*(sqrt(zz)))) # SE/sqrt of n # rnorm (numbers required, mean, sd)
sc1 = abs(sc1)
sh1 <- rnorm(100, shp1, (se_shp1*(sqrt(zz))))
sh1 = abs(sh1)
sc2 <- rnorm(100, sca2, (se_sca2*(sqrt(zzz)))) # SE/sqrt of n # rnorm (numbers required, mean, sd)
sc2 = abs(sc2)
sh2 <- rnorm(100, shp2, (se_shp2*(sqrt(zzz))))
sh2 = abs(sh2)


T = 365*12*5 #daylight hours in one year 
t = 24 #conflict observation period
N <- 0

for (i in 1:100) {
  a1 = sc1[i]*-0.504
  b1 = a1/sh1[i]
  c11 = (1+b1); c12 = (-1/sc1[i])
  c1 = c11^c12
  dx1 = (sc1[i]*-0.504)/sh1[i];  dx2 = (1 + dx1)^(-1/sc1[i]);  dx3 = 1 - dx2; 
  d1 = dx3*c1
  e1 = 1-d1
  f1 = log(e1)
  g1 = f1^-1
  h1 = 1*g1
  a2 = sc2[i]*5.262
  b2 = a2/sh2[i]
  c21 = 1+b2; c22 = (-1/sc2[i])
  c2 = c21^c22
  dx21 = (sc2[i]*5.262)/sh2[i];  dx22 = (1 + dx21)^(-1/sc2[i]);  dx23 = 1 - dx22; 
  d2 = dx23*c2
  e2 = 1-d2
  f2 = log(e2)
  g2 = f2^-1
  h2 = -1*g2
  e11 = -1/dep
  var1 = h1^e11
  var2 = h2^e11
  var3 = (var1+var2)^dep
  var4 = exp(-var3)
  R = 1-var4
  N[i] = (T/t)*R
  print(i)
}
x2 <- na.omit(N)
quantile(x2, probs = c(0.025, 0.5, 0.975))

############### Bivariate modelling ########################
x = cbind(negMTT,CIF)

f1 = fbvevd(x, model = c("log"), sym = FALSE, 
            nsscale1 =  FlowTF, nsloc2 = DensityTF, 
            cloc = FALSE, std.err = TRUE,  corr = FALSE, method = "BFGS",
            warn.inf = TRUE)
dep.summary(f1)
f1[["estimate"]]
f1[["std.err"]]

###############################