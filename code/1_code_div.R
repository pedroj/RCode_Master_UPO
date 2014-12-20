#######################################################################
# Curso R. UPO-EBD, Nov 2013
# Pedro Jordano.
#----------------------------------------------------------------------
# CODIGO R usado en el curso. OTROS ANALISIS
#######################################################################
library(MASS)
library(datasets)
library(ade4)
library(vegan)
library(FD)
install.packages("vegan")
#----------------------------------------------------------------------
# Diversidad
library(vegan)
# BCI: A data frame with 50 plots (rows) of 1 hectare with counts of 
# trees on each plot with total of 225 species (columns). 
# Full Latin names are used for tree species.
data(BCI)
dim(BCI)
head(BCI)
str(BCI)

# Shannon–Weaver Simpson
H <- diversity(BCI)  # H'= -Sum pi*log pi
H <- diversity(BCI,index = "simpson")  # D2= 1-Sum pi^2

# Pielou’s evenness J = H′/ log(S):
J <- H/log(specnumber(BCI))

# Riqueza especifica
## Species richness (S):
S <- specnumber(BCI) ## rowSums(BCI > 0) hace lo mismo...

# Renyi diversities de orden a:
# Los indices de diversidad son casos especiales de la diversidad de 
# Renyi:
# Ha = 1/(1− a) log(Sum pi^a) ; donde a es un parametro de escala
# Los numeros de Hill correpondientes son: 
# Na=exp(Ha). Hill sugirio usar los numeros 0, 1, 2 e infinito para los
# indices: S, H', inverso de Simpson y 1/max(pi), respectivamente. 
# O sea, los indices de diversidad mas usados son casos especiales de 
# numeros de Hill:
# N0 = S
# N1 = exp(H′)
# N2 = D2 (inverso del indice de Simpson): 1/Sum pi^2
# N∞ = 1/(maxpi). 
# Las diversidades de Renyi correspondientes son:
# H0 = log(S) 
# H1 = H′
# H2 = − log(Sum pi^2)
# H∞ = − log(max pi). 

# Una forma de comparar diversidades es examinar los indices de
# Renyi para un rango de valores de a, tipicamente 0, 1, 2, 4... inf
data(BCI)

# Seleccionamos 9 plots para estimar su diversidad Renyi:
i <- sample(nrow(BCI), 9)
mod <- renyi(BCI[i,])
exp(mod[,1]) # Numero de especies en cada plot
plot(mod)
mod <- renyiaccum(BCI[i,])
plot(mod, as.table=TRUE, col = c(1, 2, 2))

#----------------------------------------------------------------------
# Rarefaccion
# Riqueza esperada de especies en submuestras aleatorias de la 
# comunidad con el tamaño dado. 
# Hacer:
rowSums(BCI); specnumber(BCI)    # Esta S relacionado con Nind?
plot(rowSums(BCI),specnumber(BCI))
Srar <- rarefy(BCI, min(rowSums(BCI)))

## Rarefaction
(raremax <- min(rowSums(BCI)))
Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)

# Curvas de rarefaccion para cada plot (fila del dataset). 
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)

# Acumulacion de especies y estimadores asintoticos
sac <- specaccum(BCI)
plot(sac, ci.type = "polygon", ci.col = "yellow")
sac2<-specaccum(BCI,"random")
plot(sac, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sac2, col="yellow", add=TRUE, pch="+") # Solo para "random"
sac2

# Otra forma- BiodiversityR
BCI$site.totals <- apply(BCI,1,sum)   # Totals by plot
Accum.1<- accumresult(BCI, y= BCI, scale='site.totals', method='exact')
Accum.1
accumplot(Accum.1, ci.type="poly", col="blue",
          lwd=2, ci.lty=0, ci.col="lightblue",type="n") # With CI intervals
          
accumplot(Accum.1) # Simple plot

# Estimadores asintoticos
specpool(BCI)
# La funcion specpool usa una coleccion de sitios (p.ej., parcelas) pero hay metodos para estimar el numero de especies "no vistas" en un unico sitio. Estas funciones precisan conteos de individuos.
# Hacer:
#estimateR(BCI[k, ])  Donde k es cualquiera de los plots 1 ha 
estimateR(BCI[12, ])
# Chao, ACE, etc son estimadores de riqueza de especies no parametricos.

#----------------------------------------------------------------------
# Ver aplicaciones a otros ejemplos: genetica de poblaciones, 
# redes,etc., en el pdf de la presentacion
#----------------------------------------------------------------------




