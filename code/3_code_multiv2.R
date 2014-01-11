####################################################################
# Curso R. UPO-EBD, Nov 2013.
# Pedro Jordano.
#-------------------------------------------------------------------
# CODIGO R usado en el curso. ANALISIS MULTIVARIANTE (2).
####################################################################
library(MASS)
library(vegan) # Multivariante
library(ade4)  # Multivariante
library(datasets)
library(klaR) # Para funcion discriminante cuadratica qda( )
library(CCA)  # Correlaciones canonicas
library(yacca)  # Correlaciones canonicas
#-------------------------------------------------------------------
# Analisis de correspondencias
data(avimedi)

# avimedi es una lista con informacion sobre 302 sitios: 
# frecuencias de 51 especies de aves; 
# dos factores (habitats y Mediterranean origin). 
# La lista contiene los objetos:
# fau: data frame de 302 sitios - 51 especies de aves.
# plan: data frame de 302 sitios - 2 factores : reg con dos niveles Provence (Pr, S Francia) y Corcega (Co) ; str con seis niveles describiendo la vegetacion- desde un matorral bajo (1) hasta bosque maduro de encinas (6).
# nomesp: vector con 51 nombres latinos.
# Source: Blondel, J., Chessel, D., & Frochot, B. (1988) Bird species impoverishment, niche expansion, and density inflation in mediterranean island habitats. Ecology, 69, 1899–1917.

coa1 <- dudi.coa(avimedi$fau, scan = FALSE, nf = 3)
s.class(coa1$li,avimedi$plan$str:avimedi$plan$reg, sub = "Correspondences Analysis")
score(coa1)

# Waders. 
# Este data frame tiene 15 filas (sitios) y 19 columnas (especies de limicolas). Los datos son conteos en verano.
# S1 Oystercatcher
# S2 White-fronted Plover
# S3 Kitt Lutzs Plover
# S4 Three-banded Plover
# S5 Grey Plover
# S6 Ringed Plover
# S7 Bar-tailed Godwit
# S8 Whimbrel
# S9 Marsh Sandpiper
# S10 Greenshank
# S11 Common Sandpiper
# S12 Turnstone
# S13 Knot
# S14 Sanderling
# S15 Little Stint
# S16 Curlew Sandpiper
# S17 Ruff
# S18 Avocet
# S19 Black-winged Stilt

# Sites:
# A = Namibia North coast
# B = Namibia North wetland
# C = Namibia South coast
# D = Namibia South wetland
# E = Cape North coast
# F = Cape North wetland
# G = Cape West coast
# H = Cape West wetland
# I = Cape South coast
# J= Cape South wetland
# K = Cape East coast
# L = Cape East wetland
# M = Transkei coast
# N = Natal coast
# O = Natal wetland

plot(corresp(waders, nf=2))

#-------------------------------------------------------------------
# Analisis discriminante. Ejemplo completo.

# Ejemplo con el iris dataset. 
# El dataset tiene 150 medidas de flores (longitud de petalos y sepalos y anchura) de tres especies de lirios gen. Iris (50 especimenes de cada especie).
iris<-read.table("iris.txt",header=TRUE,sep="\t",dec=".",na.strings="NA")
# Vars: sep_len,sep_wid,pet_len,pet_wid,species

# Necesitamos MASS.
library(MASS)

# La sintaxis con MASS para el analisis discrimiante lineal es  lda(classvariable~., dataset). El punto '.' significa "todas las otras variables", pero podriamos indicarlas por su nombre tambien  lda(classvariable~var1+var2+ ...+varN, dataset). Escribimos los resultados a un archivo de salida que luego podemos ver:
output= lda(species~., iris)

# Obtenemos primero los discriminant analysis scores, o sea, las nuevas coordenadas de los puntos tras la rotacion de la matriz. Esto es exactamente igual que un PCA. Luego hacemos un plot de las dos primeras funciones discriminantes (equivalente a PC1, PC2) and coloreamos los puntos por especie. La opcion asp=1 fuerza que los dos ejes esten en la misma escala y eso ayuda a que identifiquemos que funcion discriminante es la mas efectiva separando los grupos:
scores<-predict(output, iris)$x
plot(scores, col=rainbow(3)[iris$species], asp=1)

# Plot de puntos identificados sobre las dos primeras funciones discrim. 
plot(output) # fit from lda

# Histogramas y plots de densidad para las observaciones de cada grupo sobre la primera funcion discriminante. 
plot(output, dimen=1, type="both") # fit from lda

# Vemos el archivo de salida para confirmar que todo lo lleva la primera funcion discriminante, que absorbe el 99% de la varianza (“Proportion of trace”). Tambien tenemos las funciones exactas o combinaciones lineales de las variables originales.

# Clasificando nuevas observaciones
# Colectamos 6 nuevos lirios que queremos clasificar. Medimos las flores y dejamos vacio el campo de especie. Y hacemos la clasificiacion igual que antes. This assumes that you still have R open and that you are finished with the previous exercise, which gives us the “output” file:
unknown<-read.table("iris_unknown.txt",header=TRUE,sep="\t",dec=".",na.strings="NA")
predict(output, unknown)$class

# Tambien podemos añadir estos puntos al plot discrimiante anterior:
plot(scores, asp=1, col=rainbow(3)[iris$species])
scores_unknown=predict(output, unknown)$x
points(scores_unknown, pch=19)

# Test de diferencias entre grupos
# Para ver si existen diferencias significativas entre grupos podemos hacer un MANOVA- Multivariate Analysis of Variance. De hecho esto es matematicamente identico al discriminante anterior. 

# Primero debemos cambiar la estructura de los datos. Partimos el dataset en la variable clasificatoria por un lado y las medidas por otro, definiendo ambas como matrices:
species=as.factor(as.matrix(iris[,5])) 
measurements=as.matrix(iris[,1:4])

# Corremos el MANOVA, y escribimos el resultado en un archivo de salida. Pedimos la lambda de Wilks, que es el test estadistico que rinde el MANOVA y odemos pedir tambien ANOVA de las variables individuales.
output=manova(measurements~species) 
summary(output, test="Wilks")
summary.aov(output)

# Tablas de clasificaciones correctas
x <- lda(species ~ sep_len+sep_wid+pet_len+pet_wid, data=iris)
y <- predict(x, iris)
# absolute numbers: 
errormatrix(iris$species, y$class)
# relative frequencies: 
errormatrix(iris$species, y$class, relative = TRUE)
# percentages: 
round(100 * errormatrix(iris$species, y$class, relative = TRUE), 0)

#-------------------------------------------------------------------
# Ejemplo con library(vegan). Analisis discriminante de correspondencias.

data(perthi02) # perthi02 es una list con 2 componentes.
# 
# tab: data frame 904 filas (proteinas de 201 especies) 20 columnas (amino acidos).
# cla: es un factor de 3 clases de proteina.
# Los niveles de perthi02$cla son cyto (cytoplasmic proteins) memb (integral membran proteins) y peri (periplasmic proteins)

# Source: Perriere, G. and Thioulouse, J. (2002) Use of Correspondence Discriminant Analysis to predict the subcellular location of bacterial proteins. Computer Methods and Programs in Biomedicine, 70, 2, 99–105.

d1<-discrimin.coa(perthi02$tab, perthi02$cla, scan = FALSE)

summary(d1)

# Salida   
# eig  autovalores
# nf   numero de ejes mantenidos
# fa   matriz con loadings: canonical weights
# li   dataframe con los scores canonicos
# va   cosenos (corr) entre variables
# cp   cosenos (corr) entre componentes
# gc   dataframe con los scores de las clases


# Retorna una lista de clase 'discrimin' con:
# nf	 un valor numerico indicando el numero de ejes extraidos
# eig	 un vector numerico  con todos los eigenvalues (autovalores)
# fa	 matriz con los pesos: los "canonical weights"
# li	 data frame con los "canonical scores"
# va	 matriz con los valores de los cosenos (corr) entre las variables y 
#        los canonical scores
# cp	 matriz con los valores de los cosenos (corr) entre los componentes y 
#        los canonical scores
# gc	 data frame con los scores de las clases (grupos)

plot(d1)     # Plot general

#-------------------------------------------------------------------
# Funcion discriminante cuadratica.

# La funcion discriminante cuadratica qda( ) se puede usar en vez de la lineal lda( ) cuando no asumimos homogeneidad de las matrices var-covar.

# Quadratic Discriminant Analysis con 3 grupos aplicando prior probabilities iguales. 
library(MASS); library(datasets); library(klaR)

# Usamos iris dataset. Podriamos usar tambien crabs.
# Vars: sep_len,sep_wid,pet_len,pet_wid,species
fit <- qda(species ~ sep_len+sep_wid+pet_len+pet_wid, data=na.omit(iris),prior=c(1,1,1)/3)

# La funcion partimat( ) en la libreria klaR muestra los resultados de clasificaciones cuadraticas a partir de plots de 2 variables a la vez. La funcion drawparti permite un control mas preciso de este tipo de plot.
# Graficas exploratorias para LDA o QDA
library(klaR)
partimat(species ~ sep_len+sep_wid+pet_len+pet_wid, data=iris, method="lda")

# Matriz de clasificaciones con pares de variables
partimat(species ~ ., data = iris, method = "lda", 
      plot.matrix = TRUE, imageplot = T)
# Scatterplot para 3 grupos
pairs(iris[c("sep_len","sep_wid","pet_len","pet_wid")], main="Iris", pch=22, bg=c("red", "yellow", "blue")[unclass(iris$species)])

# Tablas de clasificaciones correctas
data(iris)
library(MASS)
x <- qda(species ~ sep_len+sep_wid+pet_len+pet_wid, data=iris)
y <- predict(x, iris)
# absolute numbers: 
errormatrix(iris$species, y$class)
# relative frequencies: 
errormatrix(iris$species, y$class, relative = TRUE)
# percentages: 
round(100 * errormatrix(iris$species, y$class, relative = TRUE), 0)

# Ver las probabilidades de clasificacion correcta en relacion
# con las variables
plineplot(species ~ ., data = iris, method = "qda", 
      x = "pet_len", xlab = "Petal length")

#-------------------------------------------------------------------
# Analisis de Correlaciones Canonicas
library(vegan);library(yacca); library(CCA)

# El analisis de correlaciones canonicas es una forma de analisis de asociaciones entre dos conjuntos de variables medidos sobre las mismas unidades muestrales. Tecnicamente es un analisis de subespacio lineal, y conlleva la proyeccion de los dos conjuntos de variables en un subespacio comun. El objetivo del CCA es encontrar una secuencia de transformaciones de cada conjunto de variables de tal forma que las correlaciones entre estas variables transformadas queden maximizadas- con el requerimiento de que cada variable transformada (var. canonica) sea ortogonal a la(s) precedentes. Estas CV expresan la variacion compartida entre los dos conjuntos, de forma analoga a un PCA para el caso de un unico conjunto. 

data(mite)
data(mite.env)
str(mite)

out <- CCorA(mite, mite.env[,1:2])

# Salida principal
out
# Plots
biplot(out, "ob")                 # Dos plots de objetos
biplot(out, "v", cex=c(0.7,0.6))  # Dos plots de variables
biplot(out, "ov", cex=c(0.7,0.6)) # Los 4 plots
biplot(out, "b", cex=c(0.7,0.6))  # Dos biplots
biplot(out, plot.type="biplots", xlabs = NULL) 

#-----------------------------------------------

cca.fit <- cca(mite, mite.env[,1:2])

# Mas resultados...
summary(cca.fit)
plot(cca.fit)

# Cantidades que importan: corr entre las vars originales y las canonicas respectivas (structural correlations or loadings), los coeficientes de las vars originales en las CV, y las corrs entre los scores CV de un conjunto sobre los de otro (canonical correlations). Las correlaciones canonicas son una medida de concordancia pero es mejor usar las canonical redundancies. 
# Para interpretar las CCA, mirar los loadings de las variables originales. 

# El redundancy index expresa la varianza en un conjunto que esta explicada por el otro conjunto. Combina la var explicada en un conjunto con las corr canonicas. 

correl=matcor(mite, mite.env[,1:2])
img.matcor(correl)
img.matcor(correl,type=2)
res.cc<-cc(mite, mite.env[,1:2])
plt.cc(res.cc)
plt.cc(res.cc,type="v",var.label=TRUE)

# Test
# Test de las correlaciones canonicas (ver en detalle ?summary.cca)
F.test.cca(cca.fit)

# NOTA: para datos con mayor numero de vars que de casos, ver library(CCA), funcion rcc (Regularized canonical correlation analysis)
#-------------------------------------------------------------------






