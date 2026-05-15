# Print (simulated) sample

Simple print function for the `cpreval_sim` objects

## Usage

``` r
# S3 method for class 'cpreval_sim'
print(x, ...)
```

## Arguments

- x:

  A `cpreval_sim` object

- ...:

  Further arguments passed to or from other methods

## Value

Prints object

## Author

Stephane Guerrier

## Examples

``` r
# Samples without measurement error
sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, seed = 18)
#> Data: R = 52, R0 = 19, n = 1500
#>       R1 = 19, R2 = 0, R3 = 33, R4 = 1448
#> 
#> Assumed measurement error: alpha = 0%, alpha0 = 0%, beta  = 0% 
#> 
#> False negative rate of the official procedure: beta0 = 66.67%

# With measurement error
sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, alpha0 = 0,
alpha = 0.01, beta = 0.05, seed = 18)
#> Data: R = 56, R0 = 19, n = 1500
#>       R1 = 18, R2 = 1, R3 = 38, R4 = 1443
#> 
#> Assumed measurement error: alpha = 1%, alpha0 = 0%, beta  = 5% 
#> 
#> False negative rate of the official procedure: beta0 = 66.67%
```
