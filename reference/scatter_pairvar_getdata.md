# Deprecated. Please use plot_pairvar_getdata().

Deprecated. Please use plot_pairvar_getdata().

## Usage

``` r
scatter_pairvar_getdata(
  mod,
  xvar,
  orig.scale = TRUE,
  ci_width = 0.95,
  grid_size = 100,
  thin = TRUE,
  thin_to = 100
)
```

## Arguments

- mod:

  An object of class "compnet" created by the buildcompnet() function.

- xvar:

  Character string for the name of the trait to be used. Must match the
  trait name in the input data used to build the model.

- orig.scale:

  Logical value indicating whether to back-transform trait data to the
  original scale (TRUE) or leave them with mean zero and unit variance
  (FALSE).

- ci_width:

  A real number (0,1) describing the desired widths of credible band.
  Defaults to 0.95.

- grid_size:

  A positive integer defining the number of discrete steps to use in
  approximating the shape of mean prediction curves and credible bands.
  Defaults to 100.

- thin:

  Logical value determining whether to use a random subsample of the
  full posterior sample.

- thin_to:

  Integer value determining how many random samples to draw from the
  full posterior sample if thin=TRUE.
