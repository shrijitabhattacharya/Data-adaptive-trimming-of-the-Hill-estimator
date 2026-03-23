This code implements a framework for extreme outlier detection and adaptive robust tail inference based on the paper Bhattacharya et al 2019. The data are merely assumed to be heavy-tailed, where the largest k0 sample points may be contaminated. Such outliers can severely bias the classic tail-inference. The trimmed Hill estimator introduced in Bhattacharya et al 2019 provides asymptotically unbiased and efficient estimators of the tail index that do not depend on the top-k0 sample values. In this sense, the trimmed Hill estimator is robust and optimal among the estimators with strict upper breakdown points.

In applications, the value of k0 is unknown and should be estmiated from the data. Bhattacharya et al 2019 introduced a data-driven method for the selection of the trimming parameter k0. The resulting adaptive trimmed Hill estimator becomes a robust estimator, which adapts to the unknown degree of contamination in the extremes.

k0_estimation.R has the code for implementing the adaptive trimmed Hill estimator.

app.R gives a shiny implementation of the algorithm.
