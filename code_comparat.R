####################################################################
# Curso R. UPO-EBD, Nov 2012
# Pedro Jordano.
#-------------------------------------------------------------------
# CODIGO R usado en el curso. ANALISIS COMPARATIVO
####################################################################
# Librerias que nos van a hacer falta.
library(ape); library(ade4)
#-------------------------------------------------------------------
# Plot of different levels of inertia
x1 <- c(1,3,5,7,9,12,14,18,19)
y1 <- c(2,4,2,5,11,9,17,20,21)
x2 <- c(1,3,5,7,9,12,14,18,19)
y2 <- c(2,2,1,4,5,7,7,9,9)

plot(x1,y1,type="n",xlab="Tiempo",ylab="Rasgo fenotípico")
points(x1,y1)
abline(lm(y1~x1))
scatterutil.ellipse(x1, y1, rep(c(1,0), c(10,10)), cell = 1.5,
                    ax = F)
points(x2,y2,pch=19)
abline(lm(y2~x2), col="red")
scatterutil.ellipse(x2, y2, rep(c(1,0),c(10,10)), cell= 1.5,  
                    ax= F)
#-------------------------------------------------------------------
# Simular un modelo Ornstein-Uhlenbeck
# Primero el Brownian
x <- cumsum(c(0,rnorm(99)))
# Replicamos el Brownian model
x.br <- replicate(5, cumsum(c(0,rnorm(99))))  # Cinco replicas
var(x.br[100,])         # Varianza del Browniano

# Ornstein-Uhlenbeck
x.ou <- numeric(100)
for (i in 2:100)
    x.ou[i] <- -0.2 * x.ou[i-1] + rnorm(1)

# Creo una funcion para hacer la simulacion del OU
sim.ou <- function () {
	x <- numeric(100)
	for (i in 2:100)
	    x[i] <- -0.2 * x[i-1] + rnorm(1)
	x   # retorna el valor de x
}
x.ou <- replicate(5, sim.ou())    # La simulacion de cinco OU's
var(x.ou[100,])         # Varianza del O-U

# Plots
y1<-range(x.br)
par(mfrow=c(2,1))
matplot(x.br,ylim=y1,type="l",
        col=1,main="Brownian",xlab="Tiempo",ylab="Estado del caracter")
        abline(0,0,col="red",lwd=0.5)
matplot(x.ou,ylim=y1,type="l",	
        col=1,main="Ornstein-Uhlenbeck",
        xlab="Tiempo",ylab="Estado del caracter")
        abline(0,0,col="red",lwd=0.5)
#-------------------------------------------------------------------
# Sylvia data (Paradis 2006; Bohning-Gaese et al. 2003)
x <- paste("AJ5345", 26:49, sep="")   # Accesion numbers
x <- c("Z73494", x)
sylvia.seq <- read.GenBank(x) # Read sequences from GenBank
write.dna(sylvia.seq,"sylvia.fas",format="fasta")  # Write in fasta
                                                   # format
taxa.sylvia <- attr(sylvia.seq,"species")
names(taxa.sylvia) <- names(sylvia.seq)
sylvia.seq[c(1,24)]
taxa.sylvia[1]<-"Sylvia_atricapilla"
taxa.sylvia[24]<-"Sylvia_abyssinica"
#-------------------------------------------------------------------
# Contrasts
tree.primates <- read.tree(text="((((Homo:0.21,Pongo:0.21):0.28,Macaca:0.49):0.13,Ateles: 0.62):0.38,Galago:1.00);")
body <- c(4.09434,3.61092,2.37024,2.02815,-1.46968)
longevity <- c(4.74493,3.3322,3.3673,2.89037,2.30259)
names(body)<-names(longevity)<-c("Homo","Pongo","Macaca","Ateles","Galago")
pic.body <- pic(body,tree.primates)
pic.longevity <- pic(longevity,tree.primates)

plot(tree.primates)
nodelabels(round(pic.body,3),adj=c(0,0.5),frame="n")
nodelabels(round(pic.longevity,3),adj=c(0,2),frame="n")

par(mfrow=c(2,1))
plot(body,longevity)
plot(pic.body,pic.longevity,pch=19)
abline(0,1,lty=2)
#-------------------------------------------------------------------
# GLS - Generalized Least Squares
# corGrafen(value,phy,fixed=FALSE)
tree.primates <- read.tree(text="((((Homo:0.21,Pongo:0.21):0.28,Macaca:0.49):0.13,Ateles: 0.62):0.38,Galago:1.00);")
body <- c(4.09434,3.61092,2.37024,2.02815,-1.46968)
longevity <- c(4.74493,3.3322,3.3673,2.89037,2.30259)
names(body)<-names(longevity)<-c("Homo","Pongo","Macaca","Ateles","Galago")
bm.prim <- corBrownian(phy=tree.primates)
DF.prim <- data.frame(body,longevity)
m1 <- gls(longevity~body,
   correlation=bm.prim,data=DF.prim)
summary(m1)

# Ornstein-Uhlenbeck
ou.prim <- corMartins(1,tree.primates)
m2 <- gls(longevity~body,
   correlation=ou.prim,data=DF.prim)
summary(m2)

# GEE-Generalized Estimating Equations
compar.gee(longevity~body,
   phy= tree.primates)
#-------------------------------------------------------------------
# Particion de la varianza taxonomica
data(carnivora)
names(carnivora)
 [1] "Order"       "SuperFamily" "Family"      "Genus"      
 [5] "Species"     "FW"          "SW"          "FB"         
 [9] "SB"          "LS"          "GL"          "BW"         
[13] "WA"          "AI"          "LY"          "AM"         
[17] "IB"     
m <- lme(log10(SW) ~ 1, random = ~ 1|Order/SuperFamily/Family/Genus, data=carnivora)
v <- varcomp(m, TRUE, TRUE)
plot(v)
barplot(as.vector(v),beside=FALSE,
        names.arg=c("Order","SuperFamily","Family","Genus","Within"))

# Multivariate decomposition
tpg<-newick2phylog(write.tree(tree.primates))
variance.phylog(tpg,body)

# Autocorrelation methods
### Example from Rohlf's 2001 article:
W<- matrix(c(
0,1,1,2,0,0,0,0,
1,0,1,2,0,0,0,0,
1,1,0,2,0,0,0,0,
2,2,2,0,0,0,0,0,
0,0,0,0,0,1,1,2,
0,0,0,0,1,0,1,2,
0,0,0,0,1,1,0,2,
0,0,0,0,2,2,2,0
),8)
W <- 1/W
W[W == Inf] <- 0
y<-c(-0.12,0.36,-0.1,0.04,-0.15,0.29,-0.11,-0.06)
rohlf<-compar.cheverud(y,W)
rohlf

1-var(rohlf$residuals)/var(y)
 
# Autocorrelation
Moran.I(body,cophenetic(tree.primates))

# Estimate with ade4
library(ade4)
gearymoran(cophenetic(tree.primates),data.frame(body,longevity))

# Autocorrelogram
data(carnivora)
correl.carn <- correlogram.formula(SW ~ Order/SuperFamily/Family/Genus,
             data=carnivora)
correl.carn
plot(correl.carn,ylim=c(-0.15,0.7),main="Carnivora; Average body weight (kg)")

#-------------------------------------------------------------------
### The example in Lynch (1991)
cat("((((Homo:0.21,Pongo:0.21):0.28,",
    "Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);",
    file = "ex.tre", sep = "\n")
tree.primates <- read.tree("ex.tre")
unlink("ex.tre")
X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
compar.lynch(cbind(X, Y),
             G = vcv.phylo(tree.primates, cor = TRUE))
#-------------------------------------------------------------------
# Plotting trees
# example tree in NH format (a string)
data("landplants.newick")
landplants.newick

# get corresponding phylo object
tree.landplants <- read.tree(text = landplants.newick)

# plot tree
plot(tree.landplants, label.offset = 0.001)
plot(tree.landplants, label.offset = 0.001, type="radial")

# Labelling the nodes
trape<-read.tree(text="((Homo,Pan),Gorilla);")
 plot(trape,x.lim=c(-0.1,2.2))
nodelabels(c("6.4 Ma","5.4 Ma"),frame="c",bg="white")

R> bs.pars<-scan()
1: NA 76 34 54 74 100 56 91 74 60 63 100 100

R> bs.nj<-scan()
1: NA 74 48 68 75 100 NA 91 67 82 52 100 100

R> bs.ml<-scan()
1: NA 88 76 73 71 100 45 81 72 67 63 100 100

R> plot(tr,no.margin=TRUE)
R> nodelabels(bs.pars, adj= c(-0.2,-0.1), 
   frame= "n",cex= 0.8, font= 2)
R> nodelabels(bs.nj, adj= c(1.2,-0.5), 
   frame= "n",cex= 0.8, font= 3)
R> nodelabels(bs.ml, adj= c(1.2,0.5), 
   frame= "n",cex= 0.8)
add.scale.bar(length=0.01)

# Using zoom
data(bird.families)
zoom(bird.families,list(1:15,38:48),col=rep("grey",2),no.margin=TRUE,font=1,subtree=TRUE)

#--------------------------------------------
R> plot(tr,no.margin=TRUE)
R> nodelabels(thermo=bs.ml/100, col= "grey", bg="white")
#-------------------------------------------------------------------
### Drawing the character maps on the phylogeny
## Phylogenetic tables, PALMS
data(palm)
?palm
# Variables in palm$traits:
# h: height in meter (squared transform)
# vfruit: fruit volume in mm^3 (log-transform)
# vgrain: seed volume in mm^3 (log-transform)
# aire: spatial distribution area (km^2)
# alti: maximum altitude in meter (logged transform)

palm.phy <- newick2phylog(palm$tre)
radial.phylog(palm.phy,clabel.l=1.25)
tt<-cbind.data.frame(scalewt((palm$traits[,c(2,4:7)])))
table.phylog(tt, palm.phy, clabel.r = 0.7, f = 0.4)

#-------------------------------------------------------------------
# NOTES. Help.
# Input data directly from the clipboard
# From the Clipboard (MacOSX)
data<-read.table(pipe("pbpaste"),header=TRUE,dec=",")














   
