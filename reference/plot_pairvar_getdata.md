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
#> Chain 1: Gradient evaluation took 7.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.78 seconds.
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
#> Chain 1:  Elapsed Time: 0.012 seconds (Warm-up)
#> Chain 1:                0.012 seconds (Sampling)
#> Chain 1:                0.024 seconds (Total)
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
#> [1] "compnet uses Stan under the hood. You may see warnings from Stan alongside, \nthis message. To deal with any warnings Stan might issue, \nPlease see the links provided in Stan's output, as well as the compnet website:\nhttps://kyle-rosenblad.github.io/compnet/"

plot_pairvar_getdata(ex_compnet_phylo, xvar="phylodist")
#>             x       qlow      means      qhigh
#> 1   0.4909603 0.01031995 0.01764614 0.03339400
#> 2   0.5563972 0.01039384 0.01760982 0.03318396
#> 3   0.6218340 0.01046827 0.01757400 0.03297520
#> 4   0.6872709 0.01054322 0.01753868 0.03276771
#> 5   0.7527078 0.01061871 0.01750386 0.03256148
#> 6   0.8181446 0.01069474 0.01746953 0.03235650
#> 7   0.8835815 0.01077131 0.01743570 0.03215277
#> 8   0.9490184 0.01084842 0.01740236 0.03195028
#> 9   1.0144552 0.01092609 0.01736952 0.03174902
#> 10  1.0798921 0.01100431 0.01733717 0.03154899
#> 11  1.1453290 0.01108309 0.01730531 0.03135018
#> 12  1.2107659 0.01116243 0.01727394 0.03115258
#> 13  1.2762027 0.01124234 0.01724306 0.03095618
#> 14  1.3416396 0.01132282 0.01721267 0.03076099
#> 15  1.4070765 0.01140387 0.01718277 0.03056698
#> 16  1.4725133 0.01148550 0.01715336 0.03037416
#> 17  1.5379502 0.01156771 0.01712443 0.03018252
#> 18  1.6033871 0.01165051 0.01709600 0.02999205
#> 19  1.6688239 0.01173389 0.01706804 0.02980274
#> 20  1.7342608 0.01181787 0.01704058 0.02961459
#> 21  1.7996977 0.01190245 0.01701359 0.02942760
#> 22  1.8651345 0.01198764 0.01698709 0.02924175
#> 23  1.9305714 0.01207343 0.01696107 0.02905703
#> 24  1.9960083 0.01215983 0.01693554 0.02887345
#> 25  2.0614451 0.01224684 0.01691049 0.02869099
#> 26  2.1268820 0.01233448 0.01688591 0.02850966
#> 27  2.1923189 0.01242274 0.01686182 0.02832943
#> 28  2.2577558 0.01249620 0.01683821 0.02815031
#> 29  2.3231926 0.01255751 0.01681507 0.02797229
#> 30  2.3886295 0.01261926 0.01679242 0.02779536
#> 31  2.4540664 0.01268145 0.01677024 0.02761952
#> 32  2.5195032 0.01274410 0.01674854 0.02744477
#> 33  2.5849401 0.01280720 0.01672732 0.02727108
#> 34  2.6503770 0.01284205 0.01670658 0.02709846
#> 35  2.7158138 0.01284838 0.01668631 0.02692691
#> 36  2.7812507 0.01285486 0.01666652 0.02675641
#> 37  2.8466876 0.01286149 0.01664720 0.02658696
#> 38  2.9121244 0.01286827 0.01662835 0.02645214
#> 39  2.9775613 0.01287521 0.01660999 0.02632654
#> 40  3.0429982 0.01288229 0.01659209 0.02620186
#> 41  3.1084350 0.01288954 0.01657467 0.02607807
#> 42  3.1738719 0.01289693 0.01655772 0.02595519
#> 43  3.2393088 0.01290449 0.01654125 0.02583319
#> 44  3.3047457 0.01290634 0.01652525 0.02571209
#> 45  3.3701825 0.01290126 0.01650972 0.02559187
#> 46  3.4356194 0.01289621 0.01649466 0.02547254
#> 47  3.5010563 0.01289121 0.01648008 0.02535408
#> 48  3.5664931 0.01288624 0.01646596 0.02523648
#> 49  3.6319300 0.01288132 0.01645232 0.02511976
#> 50  3.6973669 0.01287643 0.01643915 0.02500390
#> 51  3.7628037 0.01287158 0.01642645 0.02488889
#> 52  3.8282406 0.01286677 0.01641422 0.02477474
#> 53  3.8936775 0.01286200 0.01640246 0.02466144
#> 54  3.9591143 0.01285523 0.01639117 0.02454898
#> 55  4.0245512 0.01282316 0.01638035 0.02443736
#> 56  4.0899881 0.01279123 0.01637001 0.02432658
#> 57  4.1554249 0.01275943 0.01636013 0.02421663
#> 58  4.2208618 0.01272776 0.01635072 0.02410750
#> 59  4.2862987 0.01269621 0.01634178 0.02399920
#> 60  4.3517356 0.01266479 0.01633332 0.02389172
#> 61  4.4171724 0.01263350 0.01632532 0.02378505
#> 62  4.4826093 0.01260233 0.01631779 0.02367919
#> 63  4.5480462 0.01257129 0.01631073 0.02357414
#> 64  4.6134830 0.01254038 0.01630414 0.02346988
#> 65  4.6789199 0.01250958 0.01629803 0.02336643
#> 66  4.7443568 0.01247891 0.01629238 0.02326377
#> 67  4.8097936 0.01244836 0.01628720 0.02316189
#> 68  4.8752305 0.01241794 0.01628249 0.02306081
#> 69  4.9406674 0.01238763 0.01627826 0.02296050
#> 70  5.0061042 0.01235745 0.01627449 0.02286097
#> 71  5.0715411 0.01232739 0.01627119 0.02276221
#> 72  5.1369780 0.01229744 0.01626837 0.02266423
#> 73  5.2024148 0.01226762 0.01626601 0.02256701
#> 74  5.2678517 0.01223791 0.01626413 0.02247054
#> 75  5.3332886 0.01220832 0.01626271 0.02237484
#> 76  5.3987255 0.01217885 0.01626177 0.02227989
#> 77  5.4641623 0.01214949 0.01626130 0.02218569
#> 78  5.5295992 0.01212025 0.01626131 0.02209223
#> 79  5.5950361 0.01209113 0.01626178 0.02199952
#> 80  5.6604729 0.01206212 0.01626273 0.02190754
#> 81  5.7259098 0.01200928 0.01626415 0.02181630
#> 82  5.7913467 0.01194979 0.01626604 0.02172579
#> 83  5.8567835 0.01189064 0.01626841 0.02163600
#> 84  5.9222204 0.01183183 0.01627125 0.02154694
#> 85  5.9876573 0.01177334 0.01627457 0.02155207
#> 86  6.0530941 0.01171519 0.01627836 0.02158163
#> 87  6.1185310 0.01165736 0.01628262 0.02161159
#> 88  6.1839679 0.01159987 0.01628736 0.02170969
#> 89  6.2494047 0.01154269 0.01629258 0.02180988
#> 90  6.3148416 0.01148584 0.01629827 0.02191060
#> 91  6.3802785 0.01142931 0.01630444 0.02201186
#> 92  6.4457154 0.01137310 0.01631109 0.02211365
#> 93  6.5111522 0.01131721 0.01631822 0.02221598
#> 94  6.5765891 0.01126163 0.01632582 0.02237064
#> 95  6.6420260 0.01120636 0.01633391 0.02253137
#> 96  6.7074628 0.01115141 0.01634247 0.02269332
#> 97  6.7728997 0.01109677 0.01635152 0.02285648
#> 98  6.8383366 0.01104243 0.01636104 0.02302087
#> 99  6.9037734 0.01098841 0.01637105 0.02318649
#> 100 6.9692103 0.01093468 0.01638154 0.02335336
```
