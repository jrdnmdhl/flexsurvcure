# flexsurvcure 1.3.3
- Fixes issue with quantile function generation that led to incorrect results when generating quantiles through summary.flexsurvreg

# flexsurvcure 1.3.2
- Fixes issue with incompatible vector length in quantile functions

# flexsurvcure 1.3.1
- Fixes bug where wrong function was used in quantile calculations for summary.flexsurvreg
- Fixes issue with vectorization of quantile functions

# flexsurvcure 1.3.0
- Updated to include vectorized versions of mean_mixsurv and mean_nmixsurv

# flexsurvcure 1.2.0
- Changed p function to satisfy convention that p function is P[X <= x] when lower.tail=TRUE, rather than P[X < x]
- Adds random sampling function to object returned by flexsurvcure

# flexsurvcure 1.1.0
- Added probit link option

# flexsurvcure 1.0.0
- Fixes and performance improvements to quantile & random generation functions

# flexsurvcure 0.0.2
- Fixes to cumulative hazard and RMST functions

# flexsurvcure 0.0.1
- Initial release
