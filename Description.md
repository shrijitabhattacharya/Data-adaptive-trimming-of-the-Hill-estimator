
This code implements a framework for _extreme outlier detection_ and _adaptive robust tail inference_ based on the paper 
[Bhattacharya et al 2019](https://arxiv.org/abs/1705.03088 "Shrijita's paper"). The data are merely 
assumed to be heavy-tailed, where the largest k0 sample points may be contaminated. Such outliers can 
severely bias the classic tail-inference (see the tab on the trimmed Hill and classic Hill plots). The _trimmed Hill_
estimator introduced in [Bhattacharya et al 2019](https://arxiv.org/abs/1705.03088 "Shrijita's paper") provides
asymptotically unbiased and efficient estimators of the tail index that _do not_ depend on the top-k0 sample values.
In this sense, the trimmed Hill estimator is robust and optimal among the estimators with _strict upper breakdown points_.

In applications, the value of k0 is unknown and should be estmiated from the data. [Bhattacharya et al 2019](https://arxiv.org/abs/1705.03088 "Shrijita's paper") introduced a data-driven method for the  selection of the trimming parameter k0 demostrated in the tab on the _Trimming Diagnostics_ plot.  The resulting _adaptive trimmed Hill estimator_ becomes 
a robust estimator, which adapts to the unknown degree of contamination in the extremes. All these methods are implemented in the script below.

The third tab shows in the _Trimmed Hill plot_ (red/solid), which correctly accounts or the bias if the top-k0 observations are missing. Note that the classic Hill plot
(black/dotted) can be severely biased if outliers are present.  The naive Hill plot based on (incorrectly) ignoring the missingness 
of the top-k0 observations 
is also severely biased (black/dashed).

In the side bar, one can choose from three different categories of data. The first of them is the Simulated Data set where observations can be drawn from either Pareto, Burr, Frechet and |T| distribution. The outliers are injected according to Relation 4.4 in [Bhattacharya et al 2019](https://arxiv.org/abs/1705.03088 "Shrijita's paper"). One can also control for the tail index (xi), number of injected outliers (k0_outliers) and magnitude of the outliers (L).

The second one in the sidebar is the default Calcium Data set [NIS61072.TXT -- Condroz data. Variables: Calcium, pH value](https://lstat.kuleuven.be/Wiley/). From this data set, we analyze the calcium content levels between pH levels 7-7.5. In the last option of Uploaded Data set, the user can upload any data set of their choice provide  as long as it is in  a .txt or .csv format.
