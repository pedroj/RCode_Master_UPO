####################################################################
# Curso R. UPO-EBD, Nov 2012
# Pedro Jordano.
#-------------------------------------------------------------------
# CODIGO R usado en el curso. ANALISIS MULTIVARIANTE
####################################################################
library(MASS)
library(datasets)
library(HH)
library(mvoutlier)
library(mvnormtest)
library(psych)
library(nFactors)

#-------------------------------------------------------------------
# Outliers multivariantes
# La funcion aq.plot( ) en la libreria mvoutlier permite identificar outliers multivariantes por medio de distancias de Mahalanobis (plots de las "ordered squared robust Mahalanobis distances" de las observaciones a la funcion empirica.

# Detect Outliers in the iris Data
library(mvoutlier)
iris<-read.table("iris.txt",header=TRUE,sep="\t",dec=".",na.strings="NA")
# Vars: sep_len,sep_wid,pet_len,pet_wid,species
outliers <- 
aq.plot(c("sep_len","sep_wid","pet_len","pet_wid"))
outliers # show list of outliers

#-------------------------------------------------------------------
# Normalidad multivariante, etc...
# MANOVA y otros analisis multivariantes asumen normalidad multivariante. La funcion mshapiro.test( ) en la libreria mvnormtest da un test Shapiro-Wilk test para normalidad multivariante. La entrada al test debe ser una matriz numerica.

# Test Multivariate Normality 
mshapiro.test(as.matrix(iris[,1:4]))
# Da error de singularidad. Para ver otro ejemplo:
data(EuStockMarkets)
C <- t(EuStockMarkets[15:29,1:4])
mshapiro.test(C)

# Podemos ver graficamente la normalidad multivariante con un Q-Q plot. Las distancias de Mahalanobis al centroide se distribuyen como chi-sq con p grados de libertad. El plot Q-Q de ambos debe dar una linea recta 1:1.
x <- as.matrix(iris[,1:4])       # n x p matriz numerica
center <- colMeans(x)            # centroide
n <- nrow(x); p <- ncol(x); cov <- cov(x); 
d <- mahalanobis(x,center,cov)   # distancias al centroide 
qqplot(qchisq(ppoints(n),df=p),d,
  main="QQ Plot Assessing Multivariate Normality",
  ylab="Mahalanobis D2")
abline(a=0,b=1)

# Homogeneidad de varianzas
# La funcion bartlett.test( ) es un test parametrico K-sample de la igualdad de varianzas. La funcion fligner.test( ) es no-parametrico. 
# Bartlett Test de homogeneidad de varianzas
bartlett.test(sep_len~species, data=iris)

# Figner-Killeen Test de homogeneidad de varianzas
fligner.test(sep_len~species, data=iris)

# La funcion plot.hov( ) en la libreria HH da un test grafico de homogeneidad de varianzas basado en Brown-Forsyth.

# Plot de homogeneidad de varianzas. El plot es de trellis, con tres paneles. El test resulta de un ANOVA sobre los datos del tercer panel (las desviaciones absolutas respecto a las medianas de cada grupo).
library(HH)
hov(sep_len~species, data=iris)
plot.hov(sep_len~species, data=iris)

#-------------------------------------------------------------------
# Analisis basico - PCA
library(MASS)
data(crabs)
# Loading required package: grDevices
head(crabs)
# FL: size of frontal lobe
# RW: rear width
# CL: shell length
# CW: shell width
# BD: body depth

    sp sex index   FL   RW   CL   CW   BD
1    B   M     1  8.1  6.7 16.1 19.0  7.0
2    B   M     2  8.8  7.7 18.1 20.8  7.4
3    B   M     3  9.2  7.8 19.0 22.4  7.7
4    B   M     4  9.6  7.9 20.1 23.1  8.2
5    B   M     5  9.8  8.0 20.3 23.0  8.2

# Log-transformations
lcrabs<-log(crabs[,4:8])
head(lcrabs)
        FL     RW     CL     CW     BD
1   2.0919 1.9021 2.7788 2.9444 1.9459
2   2.1748 2.0412 2.8959 3.0350 2.0015
3   2.2192 2.0541 2.9444 3.1091 2.0412

crabs.grp<-factor(c("B","b","O","o")[rep(1:4,each=50)])
lcrabs.pca<-princomp(lcrabs)
lcrabs.pca

Call:
princomp(x = lcrabs)

Standard deviations:
   Comp.1    Comp.2    Comp.3    Comp.4    Comp.5 
0.5166405 0.0746536 0.0479144 0.0248040 0.0090522 

 5  variables and  200 observations.

summary(lcrabs.pca)

Importance of components:
                        Comp.1   Comp.2    Comp.3    Comp.4     Comp.5
Standard deviation     0.51664 0.074654 0.0479144 0.0248040 0.00905219
Proportion of Variance 0.96891 0.020230 0.0083337 0.0022333 0.00029745
Cumulative Proportion  0.96891 0.989136 0.9974692 0.9997026 1.00000000

plot(lcrabs.pca)
loadings(lcrabs.pca)
biplot(lcrabs.pca)

Loadings:
   Comp.1 Comp.2 Comp.3 Comp.4 Comp.5
FL -0.452 -0.157  0.438  0.752  0.114
RW -0.387  0.911                     
CL -0.453 -0.204 -0.371        -0.784
CW -0.440        -0.672         0.591
BD -0.497 -0.315  0.458 -0.652  0.136

               Comp.1 Comp.2 Comp.3 Comp.4 Comp.5
SS loadings       1.0    1.0    1.0    1.0    1.0
Proportion Var    0.2    0.2    0.2    0.2    0.2
Cumulative Var    0.2    0.4    0.6    0.8    1.0

fit$scores # los componentes principales
lcrabs.pca$scores

splom(~lcrabs.pca$scores[,1:3],groups=crabs.grp,panel=panel.superpose,
	key=list(text=list(c("Blue male","Blue female","Orange male",
	"Orange female")),
	points= Rows(trellis.par.get("superpose.symbol"),1:4),
	columns=4))
#-------------------------------------------------------------------
# PCA con rotacion de ejes
# La funcion principal( ) en la libreria psych puede usarse para extraer y rotar los PCAs.
# Varimax Rotated Principal Components
# Retenemos 5 componentes. Los datos de entrada pueden ser una matriz original de datos (dataset) o bien una matriz de correlacion o bien una matriz de covarianza.
library(psych)
lcrabs.rpca <- principal(lcrabs, nfactors=5, rotate="varimax")
lcrabs.rpca    # Ver los resultados

#-------------------------------------------------------------------
# Como determinar cuantos factores retener...
# La libreria nFactors tiene funciones para determinar objetivamente 
# cuantos factores o componentes retener. Detalles pueden 
# encontrarse aqui: 
#  http://www.er.uqam.ca/nobel/r17165/RECHERCHE/COMMUNICATIONS/2006/IMPS/IMPS_2006.ppt#1

library(nFactors)
ev <- eigen(cor(lcrabs)) # get eigenvalues
ap <- parallel(subject=nrow(lcrabs),var=ncol(lcrabs),
  rep=100,cent=.05)
nS <- nScree(ev$values, ap$eigen$qevpea)
plotnScree(nS)


