# R package SLIDE for Structural Learning and Integrative Decomposition of Multi-View Data

More info coming soon.



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