####################################################################
# Curso R. UPO-EBD, Nov 2014
# Pedro Jordano.
#-------------------------------------------------------------------
# CODIGO R usado en el curso. ANALISIS ESPACIAL
####################################################################
# Librerias que vamos a necesitar.
library(ecespa); library(spatial); library(spatstat); library(gstat)
library(sgeostat); library(MASS)

#-------------------------------------------------------------------
# Analisis basico
# Representacion de tres patrones de puntos basicos.
data(fig1,fig2,fig3)
par(mfrow=c(1,3))
plot(fig1, pch=19); plot(fig2, pch =19); plot(fig3, pch=19)
plot(fig1, pch=19, main= "Random"); plot(fig2, pch =19, 
     main= "Clumped"); plot(fig3, pch=19, main= "Regular")
par (mfrow=c(1,1))

#---------------------------------------------------------------
# Objetos de puntos
data(swedishpines) 
X <- swedishpines
par(mfrow=c(1,3))
plot(X)
X
summary(X)
plot(density(X, 10), axes= TRUE)
contour(density(X, 10), axes= FALSE)

#-------------------------------------------------------------------
# Patrones de puntos multi-tipos
data(lansing)
lansing
head(lansing)
summary(lansing) 
plot(lansing)
plot(split(lansing))
hick <- split(lansing)$hickory 
plot(hick) 

#-------------------------------------------------------------------
# Puntos marcados. Crear ppp con puntos marcados.
fig4<-fig1[1:40,]   # Tomo 40 coordenadas de Fig 1
m<- round(rnorm(40, mean=100, sd=60),1)  # Hago 40 marcas al azar
fig4<-cbind(fig4,m) # Lo junto
fig4.ppp <- ppp(fig4$x, fig4$y, c(0,100), c(0,100), marks=fig4$m)
plot(fig4.ppp)
# Plot con transparencia
plot(fig4.ppp, bg = rgb(0, 0, 1, 0.2), pch = 16, lty=0, frame.plot=T)

#---------- EJEMPLO
# Distancia entre dos puntos
plot(fig1, pch=19)
dist<-function(x1, y1, x2, y2)  sqrt((x2 - x1)^2 + (y2 - y1)^2)

# Function to estimate nearest neighbor for n points and draw lines
# after plotting. Substitue 100 by no. of points 
# Usamos fig1, que tiene 87 puntos.
r<-numeric(87)
nn<-numeric(87)
d<-numeric(87)
for (i in 1: 87) {
d<-0
for (k in 1: 87) d[k]<-dist(fig1$x[i],fig1$y[i],fig1$x[k],fig1$y[k])
r[i]<-min(d[-i])
nn[i]<-which(d==min(d[-i]))						
}
# Dibujamos las lineas
for (i in 1: 87) lines(c(fig1$x[i],fig1$x[nn[i]]),c(fig1$y[i],
                       fig1$y[nn[i]]))   # draw the lines

#-------------------------------------------------------------------
# Patrones interactivos
X <- clickppp(10)
plot(X)
#-------------------------------------------------------------------
# Patrones ppp
x <- c(1.94,0.32,1.74,0.64,0.12,1.44,0.29,0.74,0.32,1.35,1.23,0.53,0.98,0.96,0.91,1.28,1.24,0.14,1.75,0.24,0.45,0.94,1.22,1.60,0.62)
y <- c(0.40,0.70,0.91,0.92,0.13,0.92,0.72,0.15,0.78,0.59,0.02,0.70,0.75,0.33,0.52,0.75,0.19,0.32,0.87,0.13,0.63,0.08,0.72,0.67,0.96)
m<- c(9.2,3.2,14.4,12.3,2.5,6.1,2.7,10.4,10.2,0.4,20.9,10.4,25.7,7.7,13.7,10.4,8.1,9.7,0.3,0.2,1.9,11.5,6.1,16.8,36.2)

Q <- ppp(x, y, c(0, 2), c(0, 1), marks = m) 
Q 
plot(Q) 

m <- sample(c("Yes", "No"), 25, replace = TRUE) 
m <- factor(m) 
YN <- ppp(x, y, c(0, 2), c(0, 1), marks = m) 
YN 
plot(YN)

#-------------------------------------------------------------------
# Objetos ppp
help (ppp)
fig1.ppp = haz.ppp(fig1)
cosa1=Kest(fig1.ppp)
fig2.ppp = haz.ppp(fig2)
cosa2=Kest(fig2.ppp)
fig3.ppp = haz.ppp(fig3)
cosa3=Kest(fig3.ppp)

# Creamos un objeto llamado "cosa1" en el que se almacena el resultado de estimar la funcion K del patron de fig1 con Kest. Para conocer los detalles del resultado, opciones de analisis, etc, posdemos teclear help(Kest). Para representarla graficamente, tecleamos:
plot(cosa1)
plot(cosa1, col=c(1,0,0,1), lwd=c(2,2,2,2), lty=c(1,1,1,2), main="")

# Los argumentos col, lwd y lty hacen referencia al color, anchura y tipo de las lineas; main indica que escribir como rotulo del grafico ). Lo anterior nos representa la estimacion de la funcion K con la correccion isotropica de Ripley (en trazo continuo, lty=1) y la funcion K teorica de un proceso CSR (en trazo discontinuo, lty=2).

# La funcion L(r) [= (K(r) /π )1/2)-r], se puede representar con:

plot(cosa1, sqrt(./pi)-r~r, col=c(2,1,1,3), lwd=c(2,0.2,0.2,2), lty=c(1,1,1,2), main="",ylab="L(r)")
par(mfrow=c(1,3))
plot(fig2)
plot(cosa2, col=c(1,0,0,1), lwd=c(2,2,2,2), lty=c(1,1,1,2), main="")
plot(cosa2, sqrt(./pi)-r~r, col=c(2,1,1,3), lwd=c(2,0.2,0.2,2), lty=c(1,1,1,2), main="",ylab="L(r)")
par(mfrow=c(1,3))
plot(fig3)
plot(cosa3, col=c(1,0,0,1), lwd=c(2,2,2,2), lty=c(1,1,1,2), main="")
plot(cosa3, sqrt(./pi)-r~r, col=c(2,1,1,3), lwd=c(2,0.2,0.2,2), lty=c(1,1,1,2), main="",ylab="L(r)")

# En la aproximación de Ripley para Lm, en áreas rectangulares, los valores críticos aproximados del 1 % son ± 1.68 √A / N, siendo A la superficie del área de estudio y N el número de puntos. Podemos comparar las envueltas pointwise y las envueltas simultáneas del patrón fig3 con las siguentes instrucciones:

cosa1.env=envelope(fig1.ppp,Kest)
cosa2.env=envelope(fig2.ppp,Kest)
cosa3.env=envelope(fig3.ppp,Kest)
par(mfrow=c(1,3))
plot(cosa1.env,sqrt(./pi)-r~r, lwd=c(2,1,1,1), lty=c(1,1,3,3), col=c(1,1,2,2), xlab="r", ylab="L(r)",main="")
plot(cosa2.env,sqrt(./pi)-r~r, lwd=c(2,1,1,1), lty=c(1,1,3,3), col=c(1,1,2,2), xlab="r", ylab="L(r)",main="")
plot(cosa3.env,sqrt(./pi)-r~r, lwd=c(2,1,1,1), lty=c(1,1,3,3), col=c(1,1,2,2), xlab="r", ylab="L(r)",main="")

#-------------------------------------------------------------------
# Funcion K de Ripley
data(cells)
K<-Kest(cells)
par(mfrow=c(1,3))
plot(cells)
plot(Kest(cells))
L<-Lest(cells) 
plot(L,main="Lfunction")

# Pair correlation
plot(pcf(cells))

# Simulaciones de Monte Carlo para estima de CI
data(cells) 
E <- envelope(cells, Kest, nsim = 39, rank = 1) 
E
plot(E,main="pointwise envelopes")

#-------------------------------
# Puntos inhomogeneos para los arces (maple)
data(lansing)
X <- unmark(split(lansing)$maple)
plot(X,main="Acer", cex = 0.5, pch = "+")

# (1) Funcion de intensidad estimada ajustando un modelo
# Ajuste de la tendencia espacial: polynomial para las coordenadas
# x and y
fit <- ppm(X, ~ polynom(x,y,2), Poisson())

# (a) valores de intensidad esperados en los puntos, para
# obtener un valor de lambda (vector de valores)
lambda <- predict(fit, locations=X, type="trend")

# Inhomogeneous K function
Ki <- Kinhom(X, lambda)
plot(Ki)
  
# Como hacer 'envueltas' de simulacion:
#      Example shows method (2)
smo <- density.ppp(X, sigma=0.1)
Ken <- envelope(X, Kinhom, nsim=99,
                simulate=expression(rpoispp(smo)),
                sigma=0.1, correction="trans")
par(mfrow=c(1,2))
plot(smo, main="Densidad de arboles")
plot(Ken,main="Ken")   
#-------------------------------------------------------------------
# Quadrat counts
Q<-quadratcount(X,nx=4,ny=3)
plot(X) 
plot(Q,add=TRUE,cex=2)
K<-Kest(X) 
plot(K)  # mas tarde explicamos mas sobre la K de Ripley

#-------------------------------------------------------------------
# Tests
data(bei) 
summary(bei) 

par(mfrow=c(1,2))
quadratcount(bei, nx = 6, ny = 3) 
Q <- quadratcount(bei, nx = 6, ny = 3) 
plot(bei, cex = 0.5, pch = "+") 
plot(Q, add = TRUE, cex = 2)

M <- quadrat.test(bei, nx = 6, ny = 3)
M 
plot(Q) 
plot(M, add = TRUE, cex = 1)

#-------------------------------------------------------------------
# Otros plots de densidad y tendencia

den <- density(bei, sigma = 70) 
plot(den) 
plot(bei, add = TRUE, cex = 0.5)

persp(den)
contour(den)

plot(rpoispp(100))

#-------------------------------------------------------------------
# Analisis de cuadriculas
data(nztrees) 
nztrees 
plot(nztrees)
contour(density(nztrees,10),axes=FALSE)

# "Recortamos" el mapa a fin de excluir una parte de la muestra
chopped<-owin(c(0,148),c(0,95)) 
# O bien podemos decir:
win<-nztrees$window 
chopped<-trim.rectangle(win,xmargin=c(0,5),ymargin=0) 
chopped 
nzchop<-nztrees[chopped] 
summary(nzchop) 
plot(density(nzchop,10)) 
plot(nzchop,add=TRUE)

# Cuadriculas
M <- quadrat.test(nzchop, nx = 3, ny = 2)
M 
plot(nzchop) 
plot(M, add = TRUE, cex = 2)

# KS Test
KS <- kstest(nzchop, function(x, y) {x}) 
plot(KS) 
pval <- KS$p.value

#-------------------------------------------------------------------
# Test de modelos de distribucion. Funcion ppm
# Ajusta un modelo Poisson estacionario marcado
# con diferente intensidad para cada especie.
plot(lansing)
ppm(lansing, ~ marks, Poisson())
# Se interpreta de forma analoga a un glm.

#-------------------------------------------------------------------
# Test de modelos con covariantes
data(bei)
grad <- bei.extra$grad
par(mfrow=c(1,2))
quadratcount(bei, nx = 6, ny = 3) 
Q <- quadratcount(bei, nx = 6, ny = 3) 
plot(bei, cex = 0.5, pch = "+") 
plot(Q, add = TRUE, cex = 2)
plot(grad)
fit <- ppm(bei, ~slope, covariates = list(slope = grad)) 
fitnull <- update(fit, ~1) 
anova(fitnull, fit, test = "Chi")
#-------------------------------------------------------------------
# Variogramas, triangulaciones de Delaunay

library(sgeostat); library(MASS)
data(lansing)
hick <- split(lansing)$hickory 
plot(hick) 

# Conteos de arboles por parcelas
par(mfrow=c(3,1))
Q<-quadratcount(hick, nx = 20, ny = 20) # Quadrats
plot(hick, cex = 0.5, pch = "+") 
plot(Q, cex = 0.5)
# Genero las coordenadas de las parcelas
hick.cnt<-as.data.frame(cbind(rep(1:20, each=20),rep(1:20, 20),as.vector(Q)))
colnames(hick.cnt)<-c("x","y","z")
str(hick.cnt)

# Con libreria MASS
topo.kr <- surf.ls(2, hick.cnt)
correlogram(topo.kr, 30) 
d <- seq(0, 30, 0.1)
lines(d, expcov(d,0.7))

#--------------------------------
# Delaunay triangulations
library(deldir)
hick.tr<-deldir(hick$x,hick$y)
par(mfrow=c(1,3))
plot(hick)
plot(delaunay(hick))
plot(hick.tr)
#-------------------------------------------------------------------


