# Get data for a plot of the effect of a pair-level variable like phylogenetic distance

Get data for a plot of the effect of a pair-level variable like
phylogenetic distance

## Usage

``` r
plot_pairvar_getdata(
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

## Value

A ggplot2 graphic.

## Examples

``` r
data(ex_presabs)
data(ex_phylo)

# Quick demo run. Will prompt warnings.
# Run with default warmup and iter for good posterior sampling.
ex_compnet_phylo <- buildcompnet(presabs=ex_presabs,
pairvars=ex_phylo, warmup=10, iter=20, family='binomial')
#> 
#> SAMPLING FOR MODEL 'srm_binomial' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 8.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.85 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: No variance estimation is
#> Chain 1:          performed for num_warmup < 20
#> Chain 1: 
#> Chain 1: Iteration:  1 / 20 [  5%]  (Warmup)
#> Chain 1: Iteration:  2 / 20 [ 10%]  (Warmup)
#> Chain 1: Iteration:  4 / 20 [ 20%]  (Warmup)
#> Chain 1: Iteration:  6 / 20 [ 30%]  (Warmup)
#> Chain 1: Iteration:  8 / 20 [ 40%]  (Warmup)
#> Chain 1: Iteration: 10 / 20 [ 50%]  (Warmup)
#> Chain 1: Iteration: 11 / 20 [ 55%]  (Sampling)
#> Chain 1: Iteration: 12 / 20 [ 60%]  (Sampling)
#> Chain 1: Iteration: 14 / 20 [ 70%]  (Sampling)
#> Chain 1: Iteration: 16 / 20 [ 80%]  (Sampling)
#> Chain 1: Iteration: 18 / 20 [ 90%]  (Sampling)
#> Chain 1: Iteration: 20 / 20 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.017 seconds (Warm-up)
#> Chain 1:                0.019 seconds (Sampling)
#> Chain 1:                0.036 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 2.12, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> [1] "compnet uses Stan under the hood. You may see warnings from Stan alongside, this message. To deal with any warnings Stan might issue, Please see the links provided in Stan's output, as well as the compnet website:https://kyle-rosenblad.github.io/compnet/"

plot_pairvar_getdata(ex_compnet_phylo, xvar="phylodist")
#>             x        qlow      means      qhigh
#> 1   0.4909603 0.008685444 0.01428367 0.02610390
#> 2   0.5563972 0.008710683 0.01422740 0.02584512
#> 3   0.6218340 0.008735995 0.01417175 0.02558886
#> 4   0.6872709 0.008761380 0.01411670 0.02533510
#> 5   0.7527078 0.008786837 0.01406226 0.02508381
#> 6   0.8181446 0.008812369 0.01400841 0.02483498
#> 7   0.8835815 0.008837973 0.01395515 0.02458858
#> 8   0.9490184 0.008863652 0.01390248 0.02434459
#> 9   1.0144552 0.008889404 0.01385039 0.02410298
#> 10  1.0798921 0.008915231 0.01379887 0.02386373
#> 11  1.1453290 0.008941132 0.01374792 0.02362682
#> 12  1.2107659 0.008967107 0.01369754 0.02339223
#> 13  1.2762027 0.008993158 0.01364772 0.02315993
#> 14  1.3416396 0.009019283 0.01359846 0.02292991
#> 15  1.4070765 0.009045484 0.01354974 0.02270214
#> 16  1.4725133 0.009071759 0.01350157 0.02247660
#> 17  1.5379502 0.009098111 0.01345394 0.02225328
#> 18  1.6033871 0.009124538 0.01340684 0.02203214
#> 19  1.6688239 0.009151042 0.01336028 0.02181317
#> 20  1.7342608 0.009177622 0.01331424 0.02159635
#> 21  1.7996977 0.009204278 0.01326873 0.02138166
#> 22  1.8651345 0.009231011 0.01322373 0.02116908
#> 23  1.9305714 0.009257821 0.01317924 0.02095858
#> 24  1.9960083 0.009284708 0.01313526 0.02075015
#> 25  2.0614451 0.009311672 0.01309178 0.02054377
#> 26  2.1268820 0.009338714 0.01304881 0.02033942
#> 27  2.1923189 0.009365834 0.01300633 0.02013708
#> 28  2.2577558 0.009393031 0.01296433 0.01993673
#> 29  2.3231926 0.009420307 0.01292283 0.01973835
#> 30  2.3886295 0.009447662 0.01288180 0.01956558
#> 31  2.4540664 0.009475095 0.01284126 0.01940563
#> 32  2.5195032 0.009502607 0.01280118 0.01924712
#> 33  2.5849401 0.009530198 0.01276158 0.01909005
#> 34  2.6503770 0.009557869 0.01272244 0.01893441
#> 35  2.7158138 0.009585619 0.01268376 0.01878018
#> 36  2.7812507 0.009613449 0.01264554 0.01862735
#> 37  2.8466876 0.009641358 0.01260777 0.01847590
#> 38  2.9121244 0.009669349 0.01257045 0.01832583
#> 39  2.9775613 0.009697419 0.01253357 0.01817712
#> 40  3.0429982 0.009725570 0.01249714 0.01802976
#> 41  3.1084350 0.009753803 0.01246114 0.01788373
#> 42  3.1738719 0.009782116 0.01242558 0.01773904
#> 43  3.2393088 0.009810511 0.01239045 0.01759566
#> 44  3.3047457 0.009838987 0.01235574 0.01745358
#> 45  3.3701825 0.009867545 0.01232145 0.01731279
#> 46  3.4356194 0.009896185 0.01228759 0.01717328
#> 47  3.5010563 0.009924908 0.01225414 0.01703503
#> 48  3.5664931 0.009953713 0.01222110 0.01689805
#> 49  3.6319300 0.009982601 0.01218847 0.01676231
#> 50  3.6973669 0.010011572 0.01215624 0.01662780
#> 51  3.7628037 0.010040626 0.01212441 0.01649451
#> 52  3.8282406 0.010069763 0.01209298 0.01636244
#> 53  3.8936775 0.010098985 0.01206195 0.01623157
#> 54  3.9591143 0.010128290 0.01203131 0.01610188
#> 55  4.0245512 0.010151561 0.01200105 0.01597338
#> 56  4.0899881 0.010164827 0.01197118 0.01584605
#> 57  4.1554249 0.010178191 0.01194169 0.01571987
#> 58  4.2208618 0.010182456 0.01191258 0.01559484
#> 59  4.2862987 0.010176566 0.01188384 0.01547095
#> 60  4.3517356 0.010171082 0.01185547 0.01534818
#> 61  4.4171724 0.010105297 0.01182747 0.01522653
#> 62  4.4826093 0.010017268 0.01179984 0.01510599
#> 63  4.5480462 0.009930394 0.01177257 0.01498654
#> 64  4.6134830 0.009839844 0.01174565 0.01486818
#> 65  4.6789199 0.009739729 0.01171910 0.01475090
#> 66  4.7443568 0.009640742 0.01169289 0.01463468
#> 67  4.8097936 0.009542872 0.01166704 0.01451952
#> 68  4.8752305 0.009446105 0.01164153 0.01440541
#> 69  4.9406674 0.009350428 0.01161637 0.01429234
#> 70  5.0061042 0.009255830 0.01159155 0.01418030
#> 71  5.0715411 0.009162296 0.01156707 0.01406927
#> 72  5.1369780 0.009069816 0.01154292 0.01395926
#> 73  5.2024148 0.008978377 0.01151911 0.01385024
#> 74  5.2678517 0.008887967 0.01149563 0.01377774
#> 75  5.3332886 0.008798574 0.01147247 0.01373095
#> 76  5.3987255 0.008710187 0.01144964 0.01368448
#> 77  5.4641623 0.008622793 0.01142713 0.01363831
#> 78  5.5295992 0.008536381 0.01140495 0.01359244
#> 79  5.5950361 0.008450940 0.01138307 0.01356316
#> 80  5.6604729 0.008366459 0.01136152 0.01355468
#> 81  5.7259098 0.008282926 0.01134027 0.01354626
#> 82  5.7913467 0.008200331 0.01131934 0.01353788
#> 83  5.8567835 0.008118662 0.01129871 0.01352956
#> 84  5.9222204 0.008037910 0.01127838 0.01352129
#> 85  5.9876573 0.007958062 0.01125836 0.01351307
#> 86  6.0530941 0.007879110 0.01123864 0.01350490
#> 87  6.1185310 0.007801042 0.01121921 0.01349678
#> 88  6.1839679 0.007723848 0.01120008 0.01349636
#> 89  6.2494047 0.007647518 0.01118124 0.01352167
#> 90  6.3148416 0.007572042 0.01116269 0.01354707
#> 91  6.3802785 0.007497411 0.01114443 0.01357257
#> 92  6.4457154 0.007423614 0.01112645 0.01359815
#> 93  6.5111522 0.007350641 0.01110876 0.01362384
#> 94  6.5765891 0.007278484 0.01109135 0.01364961
#> 95  6.6420260 0.007207132 0.01107422 0.01367549
#> 96  6.7074628 0.007136577 0.01105736 0.01370145
#> 97  6.7728997 0.007066808 0.01104078 0.01372751
#> 98  6.8383366 0.006997818 0.01102447 0.01375367
#> 99  6.9037734 0.006929596 0.01100843 0.01377992
#> 100 6.9692103 0.006862135 0.01099266 0.01380627
```
