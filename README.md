# R package SLIDE for Structural Learning and Integrative Decomposition of Multi-View Data

SLIDE (Structural Learning and Integrative DEcomposition) is an R package that implements methods described in 

Gaynanova, Irina, and Gen Li. "Structural learning and integrative decomposition of multi‐view data." Biometrics 75.4 (2019): 1121-1132.

  * [Gaynanova, Irina, and Gen Li. *"Structural learning and integrative decomposition of multi‐view data."* Biometrics 75.4 (2019): 1121-1132.](https://onlinelibrary.wiley.com/doi/full/10.1111/biom.13108)

To install from Github:

``` install
devtools::install_github("irinagain/SLIDE")
```

The main functions are `slide`(... To be filled ...) Each function has a documentation with a simple example which can be accessed using standard ? commands in R (i.e. `?slide`).

Please feel free to contact me at irinag [at] stat [dot] tamu [dot] edu if you have any questions or experience problems with the package.

Example
-------
```{r}
n = 100
p1 = 25
p2 = 25
data = generateModel1(n = n, pvec = c(p1, p2))
out_slide = slide(X = data$X, pvec = c(p1,p2))
out_slide$S
```
