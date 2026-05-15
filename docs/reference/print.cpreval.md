# Print estimation results

Simple print function for the four estimators (i.e. mle, conditional
mle, moment based and survey sample).

## Usage

``` r
# S3 method for class 'cpreval'
print(x, ...)
```

## Arguments

- x:

  A `cpreval` object

- ...:

  Further arguments passed to or from other methods

## Value

Prints object

## Author

Stephane Guerrier

## Examples

``` r
X = sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, seed = 18)
survey_mle(X$R, X$n)
#> Method: Survey MLE
#> 
#> Estimated proportion: 3.4667%
#> Standard error      : 0.4723%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 2.5409% - 4.3924%
#> Clopper-Pearson    : 2.5997% - 4.5214%
#> 
#> Assumed measurement error: alpha = 0%, beta = 0% 
#> Sampling: Random
```
