# Convergent Cross Mapping algorithm in MATLAB
This repository contains the implementation of the algortihm for Convergent Cross Mapping (CCM) used in the article:

Mønster, D., Fusaroli, R., Tylén, K., Roepstorff, A., & Sherson, J. F. (2017). Causal inference from noisy time-series data—Testing the Convergent Cross-Mapping algorithm in the presence of noise and external influence. *Future Generation Computer Systems*, 73, 52-62. DOI: [10.1016/j.future.2016.12.009](https://doi.org/10.1016/j.future.2016.12.009)

## How to use
In order to use the CCM function `xmap()` on data from two time series `X`  and
`Y` these have to be phase space embedded, using the time delayed coordinates method. This can be done using the function `psembed()`. An example of how to use these functions can be seen in the file `example.m`, which should produce a plot
that looks roughly like the figure below. This is the result of CCM on coupled logistic maps with unidirectional coupling. The correlation coefficient between the observed *X* and the cross-mapped estimates of *X* (in blue) are observed to
converge to high values as the library size *L* increases, indicating a causal influence of *X* on *Y*. The corresponding
correlation coefficient between observed *Y* and cross-mapped estimates of *Y* are consistently low, and do not display any convergence.

![Output from example.m](/example_plot.png)

## Documentation
### Phase space embedding
The function `psembed` returns the phase space embedding *MX* of a time series *X* in
*m* dimensions using a time delay *tau*. The time delay is measured in dimensionless
units, i.e., in terms of the vector index.
```
MX = psembed(X, m, tau);
```

### Cross mapping
The function `xmap` takes as input two time series vectors *X* and *Y* and their
phase space embeddings *MX* and *MY*, the embedding parameters *m* and
 *tau*, the library length *L* and the sampling method (either 'linear' or 'random')
 and returns the cross-mapped estimates *X_MY* and *Y_MX* along with the corresponding
 points *X1* and *Y1* from the original time series.
```
[ X_MY, Y_MX, X1, Y1] = xmap( X, Y, MX, MY, m, tau, L,'random');
```

For further details, see the MATLAB implementation source code which is documented,
as well as the article. A [preprint](https://arxiv.org/abs/1603.01155) of an earlier version of the article is also available
from arXiv.
