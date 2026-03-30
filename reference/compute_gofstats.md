# Compute summary statistics quantifying row/column-level and third order dependencies in a symmetric matrix with NA diagonals

Compute summary statistics quantifying row/column-level and third order
dependencies in a symmetric matrix with NA diagonals

## Usage

``` r
compute_gofstats(Y)
```

## Arguments

- Y:

  A symmetric matrix with NA diagonals

## Value

A named vector containing: 1- the standard deviation of row means, and
2- The triadic dependency metric used by Hoff, Fosdick, & Volfovsky's
"amen" package.

## Details

Internal helper function.
