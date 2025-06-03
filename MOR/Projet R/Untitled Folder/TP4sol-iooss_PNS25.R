
#####################################################################
# TP4 - Analyse de sensibilite, Criblage
# Bertrand Iooss
# Polytech Nice Sophia

# S O L U T I O N S

#####################################################################


##############################
# 1 M?thode de Morris sur la wing weight function

rm(list=ls())
graphics.off()

wingweight <- function(xx)
{
  # y  = wing weight
  # xx = c(Sw, Wfw, A, Lam, q, lam, tc, Nz, Wdg, Wp)
  
  # on met les entrees dans leur domaine de variation
  Sw      <- xx[,1]*50+150
  Wfw     <- xx[,2]*80+220
  A       <- xx[,3]*4+6
  Lam <- (xx[,4]*20-10)*pi/180 # transf en rad
  q       <- xx[,5]*29+16
  lam     <- xx[,6]*0.5+0.5
  tc      <- xx[,7]*0.1+0.08
  Nz      <- xx[,8]*3.5+2.5
  Wdg     <- xx[,9]*800+1700
  Wp      <- xx[,10]*0.055+0.025
  fact1 <- 0.036 * Sw^0.758 * Wfw^0.0035
  fact2 <- (A / ((cos(Lam))^2))^0.6
  fact3 <- q^0.006 * lam^0.04
  fact4 <- (100*tc / cos(Lam))^(-0.3)
  fact5 <- (Nz*Wdg)^0.49
  term1 <- Sw * Wp
  y <- fact1*fact2*fact3*fact4*fact5 + term1
  return(y)
}

namesWW <- c("Sw", "Wfw", "A", "Lam", "q", "lam", "tc", "Nz", "Wdg", "Wp")

# Methode de Morris

library(sensitivity)
?morris

mor1 <- morris(model=wingweight, factors=10, r=10,
           design=list(type = "oat", levels = 4, grid.jump = 1))
x11()
plot(mor1,xlim=c(0,120))

mor1 <- morris(model=wingweight, factors=namesWW, r=10,
              design=list(type = "oat", levels = 4, grid.jump = 1))
x11()
plot(mor1,xlim=c(0,120))


################################################
# 2 Outils graphiques d'analyse de sensibilite

# 2a

library(mlbench)

data(BostonHousing2)
?BostonHousing2

a <- BostonHousing2
b <- a[,-c(1,2,3,4,5,10,18)] # on enleve les infos spatiales, chas, b, et la sortie inutile
yBH <- b$cmedv # sortie
xBH <- b[,-1] # entr?es (on enleve cmedv)
dataBH <- data.frame(xBH, yBH)

#2b

panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.1) )
  h <- hist(x, plot = FALSE, breaks=20)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = rgb(0,124/255,124/255,0.8))
  #lines(density(x),col=1, lwd=2)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  cols=colorRampPalette(c(brewer.pal(n=9,"Blues")[8], "darkgrey", brewer.pal(n=14,"YlOrRd")[9]))(200)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="pairwise.complete.obs")
  col=cols[floor(r*100)+100]
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r)+1, col=col)
}

library(RColorBrewer)
x11()
pairs(dataBH, panel=panel.smooth, diag.panel=panel.hist, upper.panel=panel.cor, horOdd=F, verOdd=T,
      cex=.7, pch=".", bg="light blue", cex.labels=1.5, font.labels=2, cex.axis=1)

################################################
# 3 Calcul d'indices de sensibilite (mesures d'importance basees sur la variance)

set.seed(123456)

#3 a

n <- 100
d <- 10
x <- matrix(runif(n=n*d), nrow=n, ncol=d)
colnames(x) = namesWW
y = wingweight(x)

library(sensitivity)
help(package="sensitivity")
?src

srcWW <- src(x, y, nboot = 100)
print(srcWW)
x11() ; plot(srcWW) ; abline(h=0)
print(srcWW$SRC$original^2)
sum(srcWW$SRC$original^2)

# validation du modele lineaire
formule <- as.formula( y ~ Sw + Wfw + A + Lam + q + lam + tc + Nz + Wdg + Wp)
modele <- lm(formule,data=data.frame(x,y))
summary(modele)

# 3 b

srcBH <- src(xBH, yBH, nboot = 100)
print(srcBH)
x11() ; plot(srcBH) ; abline(h=0)
print(srcBH$SRC$original^2)
sum(srcBH$SRC$original^2)

# validation du modele lineaire
formule <- as.formula(yBH ~ crim + zn + indus + nox + rm + age + dis + rad + tax + ptratio + lstat)
modele <- lm(formule,data=data.frame(xBH,yBH))
summary(modele)$r.squared

# 3 c

library(car)
?vif
vif(modele)

library(sensitivity)
?lmg
lmgBH <- lmg(xBH, yBH)
print(lmgBH)
print(cbind(lmgBH$lmg,srcBH$SRC$original^2))
x11() ; plot(lmgBH)

# 3 d

# modele lineaire avec les 3 entrees influentes par SRC
formule <- as.formula(yBH ~ dis + rad + lstat)
modele1 <- lm(formule,data=data.frame(xBH,yBH))
summary(modele1)$r.squared

# modele lineaire avec les 3 entrees influentes par LMG
formule <- as.formula(yBH ~ rm + ptratio + lstat)
modele2 <- lm(formule,data=data.frame(xBH,yBH))
summary(modele2)$r.squared
