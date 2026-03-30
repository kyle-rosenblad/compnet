# Inverse logit function

Inverse logit function

## Usage

``` r
expit(x)
```

## Arguments

- x:

  A numeric value, or a vector or array of numeric values.

## Value

The inverse-logit(s) of the supplied numeric value(s).

## Examples

``` r
expit(-1)
#> [1] 0.2689414

exvec <- c(-2:2)
expit(exvec)
#> [1] 0.1192029 0.2689414 0.5000000 0.7310586 0.8807971

exmat <- matrix(c(1:12), nrow=3)
expit(exmat)
#>           [,1]      [,2]      [,3]      [,4]
#> [1,] 0.7310586 0.9820138 0.9990889 0.9999546
#> [2,] 0.8807971 0.9933071 0.9996646 0.9999833
#> [3,] 0.9525741 0.9975274 0.9998766 0.9999939

```
