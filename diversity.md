Curso R. UPO-EBD, Nov 2013. Pedro Jordano. 
------------------------------------------

Here we include the `R` code used in the practices related to diversity analysis.

Diversity
=========
First, let us load the packages that will be used in this script.

```r
require(MASS)
require(datasets)
require(ade4)
require(BiodiversityR)
require(vegan)
require(FD)
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
$N0 = S$ ; 
$N1 = exp(H')$; 
$N2 = D2$ (inverse of Simpson's index): $1/\sum {p_i}^2$; and 
$N\infty = 1/(max p_i)$

The corresponding Renyi diversities are:
$H_0 = log(S)$; 
$H_1 = H'$; 
$H_2 = -log(\sum {pi}^2)$; and 
$H_{\infty}= − log(max p_i)$

A way to compare diversities is to examine the Renyi indexes for a range of $a$ values, typically $0, 1, 2, 4 ... \infty$.

We select 9 BCI subplots to estimate their diversity:


```r
data(BCI)
i <- sample(nrow(BCI), 9)
mod <- renyi(BCI[i, ])
plot(mod)
```

![plot of chunk Renyi](figure/Renyi1.png) 

```r
mod <- renyiaccum(BCI[i, ])
plot(mod, as.table = TRUE, col = c(1, 2, 2))
```

![plot of chunk Renyi](figure/Renyi2.png) 


Rarefaction
---------
Expected species richness in random subsamples- of a given size- of the community.
 

```r
rowSums(BCI)
```

```
##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18 
## 448 435 463 508 505 412 416 431 409 483 401 366 409 438 462 437 381 347 
##  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36 
## 433 429 408 418 340 392 442 407 417 387 364 475 421 459 436 447 601 430 
##  37  38  39  40  41  42  43  44  45  46  47  48  49  50 
## 435 447 424 489 402 414 407 409 444 430 425 415 427 432
```

```r
specnumber(BCI)  # Esta S relacionado con Nind?
```

```
##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18 
##  93  84  90  94 101  85  82  88  90  94  87  84  93  98  93  93  93  89 
##  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36 
## 109 100  99  91  99  95 105  91  99  85  86  97  77  88  86  92  83  92 
##  37  38  39  40  41  42  43  44  45  46  47  48  49  50 
##  88  82  84  80 102  87  86  81  81  86 102  91  91  93
```

```r
plot(rowSums(BCI), specnumber(BCI))
```

![plot of chunk Rarefaction](figure/Rarefaction.png) 

```r
Srar <- rarefy(BCI, min(rowSums(BCI)))
```


### Species accumulation and asymptotic estimators


```r
sac <- specaccum(BCI)
```

```
## Warning: el desvío estándar es cero
```

```r
plot(sac, ci.type = "polygon", ci.col = "yellow")
```

![plot of chunk Species accumulation](figure/Species accumulation1.png) 

```r
sac2 <- specaccum(BCI, "random")
plot(sac, ci.type = "poly", col = "blue", lwd = 2, ci.lty = 0, ci.col = "lightblue")
boxplot(sac2, col = "yellow", add = TRUE, pch = "+")
```

![plot of chunk Species accumulation](figure/Species accumulation2.png) 

```r
sac2
```

```
## Species Accumulation Curve
## Accumulation method: random, with 100 permutations
## Call: specaccum(comm = BCI, method = "random") 
## 
##                                                                       
## Sites     1.000   2.000   3.000   4.000   5.00   6.000   7.000   8.000
## Richness 90.730 120.760 138.320 150.080 158.39 164.890 170.720 175.220
## sd        7.326   7.555   7.391   6.862   6.05   5.699   5.063   4.766
##                                                                         
## Sites      9.000  10.000  11.000  12.000  13.000  14.000  15.000  16.000
## Richness 179.160 182.530 185.490 188.090 190.540 192.510 194.460 196.190
## sd         4.722   4.523   4.428   4.321   3.981   4.054   4.152   4.182
##                                                                        
## Sites     17.000  18.000  19.000  20.000  21.00  22.000  23.000  24.000
## Richness 197.940 199.410 200.950 202.380 203.51 204.570 205.690 206.740
## sd         4.131   3.901   4.066   3.813   3.68   3.565   3.515   3.661
##                                                                         
## Sites     25.000  26.000  27.000  28.000  29.000  30.000  31.000  32.000
## Richness 207.810 208.960 210.000 211.030 211.950 212.760 213.460 214.210
## sd         3.603   3.646   3.559   3.424   3.433   3.204   3.135   3.059
##                                                                        
## Sites     33.000  34.000  35.000  36.000  37.000  38.00  39.000  40.000
## Richness 214.940 215.770 216.500 217.230 217.890 218.54 219.260 220.030
## sd         3.021   2.912   2.725   2.585   2.412   2.41   2.448   2.153
##                                                                        
## Sites     41.000  42.000  43.000  44.000  45.00  46.000  47.000  48.000
## Richness 220.500 221.130 221.710 222.370 222.86 223.190 223.680 224.150
## sd         2.115   2.058   2.022   1.686   1.45   1.376   1.246   1.009
##                      
## Sites     49.0000  50
## Richness 224.5200 225
## sd         0.8466   0
```


Another way, using the `BiodiversityR` package.

```r
BCI$site.totals <- apply(BCI, 1, sum)  # Totals by plot
Accum.1 <- accumresult(BCI, y = BCI, scale = "site.totals", method = "exact")
Accum.1
```

```
## Species Accumulation Curve
## Accumulation method: exact
## Call: specaccum(comm = x, method = method, permutations = permutations,      conditioned = conditioned, gamma = gamma) 
## 
##                                                                      
## Sites    429.140 858.280 1287.420 1716.560 2145.700 2574.840 3003.980
## Richness  91.780 122.610  140.046  151.712  160.236  166.831  172.142
## sd         6.935   7.193    7.002    6.636    6.248    5.895    5.591
##                                                                      
## Sites    3433.120 3862.260 4291.40 4720.540 5149.680 5578.820 6007.96
## Richness  176.553  180.306  183.56  186.429  188.991  191.306  193.42
## sd          5.337    5.126    4.95    4.803    4.677    4.567    4.47
##                                                                    
## Sites    6437.100 6866.2 7295.380 7724.520 8153.660 8582.8 9011.940
## Richness  195.355  197.1  198.819  200.381  201.849  203.2  204.543
## sd          4.382    4.3    4.222    4.147    4.073    4.0    3.927
##                                                                   
## Sites    9441.080 9870.220 10299.360 10728.500 11157.640 11586.780
## Richness  205.785  206.967   208.094   209.170   210.199   211.184
## sd          3.853    3.779     3.703     3.627     3.549     3.469
##                                                                     
## Sites    12015.920 12445.060 12874.200 13303.340 13732.480 14161.620
## Richness   212.130   213.037   213.909   214.747   215.554   216.331
## sd           3.388     3.306     3.221     3.135     3.047     2.957
##                                                                     
## Sites    14590.760 15019.900 15449.040 15878.180 16307.320 16736.460
## Richness   217.079   217.800   218.495   219.165   219.811   220.434
## sd           2.864     2.769     2.671     2.571     2.466     2.358
##                                                                     
## Sites    17165.600 17594.740 18023.880 18453.020 18882.160 19311.300
## Richness   221.035   221.616   222.175   222.716   223.237   223.740
## sd           2.245     2.127     2.002     1.871     1.729     1.576
##                                                       
## Sites    19740.440 20169.580 2.060e+04 2.103e+04 21457
## Richness   224.225   224.693 2.251e+02 2.256e+02   226
## sd           1.408     1.216 9.916e-01 6.954e-01     0
```

```r
accumplot(Accum.1, ci.type = "poly", col = "blue", lwd = 2, ci.lty = 0, ci.col = "lightblue", 
    type = "n")  # With CI intervals
```

![plot of chunk BiodiversityR](figure/BiodiversityR1.png) 

```r

accumplot(Accum.1)  # Simple plot
```

![plot of chunk BiodiversityR](figure/BiodiversityR2.png) 


### Asymptotic estimators

```r
specpool(BCI)
```

```
##     Species  chao chao.se jack1 jack1.se jack2  boot boot.se  n
## All     226 237.6   6.659 246.6    5.651 248.9 236.7   3.469 50
```


Function `specpool` uses a collection of sites (i.e., plots or subplots); but there are methods to estimate the number of "unseen" species in a single site. These functions require the count of individuals.

Just do: 

```r
estimateR(BCI[k, ])  # Not run
```

```
## Error: objeto 'k' no encontrado
```


where $k$ is any of the 1 ha plots in the BCI dataset. Chao, ACE, etc., are non-parametric estimators of species richness.


```r
estimateR(BCI[12, ])  # For subplot no. 12
```

```
##               12
## S.obs     85.000
## S.chao1  111.714
## se.chao1  12.997
## S.ACE    129.001
## se.ACE     6.843
```


See the application to other examples (population genetics, complex networks, etc.) in the pdf files of the presentations.

