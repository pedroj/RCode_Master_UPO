Curso R. UPO-EBD, Nov 2012. Pedro Jordano. DIVERSITY
------------------------------------------------------

Here we include the `R` code used in the practices related to diversity analysis.

Diversity
---------
First, let us load the packages that will be used in this script.

```r
library(MASS)
library(datasets)
library(ade4)
library(BiodiversityR)
library(vegan)
library(FD)
```


Now let's start with the estimation of several diversity indexes. First, the Shannon–Weaver $H'$:
$H'= -\sum p_i log p_i$

The Simpson's index:
$D_2= 1-\sum {p_i}^2$

Pielou’s evenness:
$J = H'/log(S)$

and Fisher's alpha ($\alpha$ parameter for the logseries distribution)


```r
data(BCI)  # We load the BCI plot dataset
H <- diversity(BCI)  # H'
H <- diversity(BCI, index = "simpson")  # Simpson's
J <- H/log(specnumber(BCI))  # Evenness
alpha <- fisher.alpha(BCI)  # Fisher's alpha
```


Renyi diversities of order $a$:
$H_a = \frac{1}{(1−a)} log(\sum {p_i}^a)$ ; where $a$ is a scale parameter.

The corresponding Hill's numbers are: 
$Na=exp(Ha)$

The diversity indexes used most often are special cases of the Hill's numbers:
$N0 = S$
$N1 = exp(H')$
$N2 = D2$ (inverse of Simpson's index): $1/\sum {p_i}^2$
$N\infty = 1/(max p_i)$

The corresponding Renyi diversities are:
$H_0 = log(S)$ 
$H_1 = H'$
$H_2 = -log(\sum {pi}^2)$
$H_{\infty}= − log(max p_i)$

A way to compare diversities is to examine the Renyi indexes for a range of $a$ values, typically $0, 1, 2, 4 ... \infty$.

We select 9 BCI subplots to estimate their diversity:


```r
data(BCI)
i <- sample(nrow(BCI), 9)
mod <- renyi(BCI[i, ])
plot(mod)
mod <- renyiaccum(BCI[i, ])
plot(mod, as.table = TRUE, col = c(1, 2, 2))
```


Rarefaction
---------
Expected species richness in random subsamples- of a given size- of the community.
 

```r
rowSums(BCI)
specnumber(BCI)  # Esta S relacionado con Nind?
plot(rowSums(BCI), specnumber(BCI))
Srar <- rarefy(BCI, min(rowSums(BCI)))
```


### Species accumulation and asymptotic estimators


```r
sac <- specaccum(BCI)
plot(sac, ci.type = "polygon", ci.col = "yellow")
sac2 <- specaccum(BCI, "random")
plot(sac, ci.type = "poly", col = "blue", lwd = 2, ci.lty = 0, ci.col = "lightblue")
boxplot(sac2, col = "yellow", add = TRUE, pch = "+")
sac2
```


Another way, using the `BiodiversityR` package.

```r
BCI$site.totals <- apply(BCI, 1, sum)  # Totals by plot
Accum.1 <- accumresult(BCI, y = BCI, scale = "site.totals", method = "exact")
Accum.1
accumplot(Accum.1, ci.type = "poly", col = "blue", lwd = 2, ci.lty = 0, ci.col = "lightblue", 
    type = "n")  # With CI intervals

accumplot(Accum.1)  # Simple plot
```


### Asymptotic estimators

```r
specpool(BCI)
```


Function `specpool` uses a collection of sites (i.e., plots or subplots); but there are methods to estimate the number of "unseen" species in a single site. These functions require the count of individuals.

Just do: 

```r
estimateR(BCI[k, ])  # Not run
```


where $k$ is any of the 1 ha plots in the BCI dataset. Chao, ACE, etc., are non-parametric estimators of species richness.


```r
estimateR(BCI[12, ])  # For subplot no. 12
```


See the application to other examples (population genetics, complex networks, etc.) in the pdf files of the presentations.


