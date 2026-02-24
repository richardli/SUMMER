# Internal functions for population simulation

Functions for calculating valuable quantities and for drawing from
important distributions for population simulation.

## Usage

``` r
getExpectedNperEA(
  easpa,
  pop.mat,
  level = c("grid", "EA"),
  pixel.index.mat = NULL
)

getSortIndices(
  i,
  urban = TRUE,
  pop.mat,
  stratify.by.urban = TRUE,
  validation.pixel.I = NULL
)

rStratifiedMultnomial(n, pop.mat, easpa, stratify.by.urban = TRUE)

rStratifiedMultnomialBySubarea(
  n,
  pop.mat,
  easpa,
  stratify.by.urban = TRUE,
  poppsub = NULL,
  min1.per.subarea = TRUE,
  min.sample = 1
)

rMyMultinomial(
  n,
  i,
  stratify.by.urban = TRUE,
  urban = TRUE,
  pop.mat = NULL,
  easpa = NULL,
  min1.per.subarea = FALSE,
  method = c("mult1", "mult", "indepMH"),
  min.sample = 1
)

rMyMultinomialSubarea(
  n,
  i,
  easpsub,
  stratify.by.urban = TRUE,
  urban = TRUE,
  pop.mat = NULL
)

rmultinom1(
  n = 1,
  size,
  prob,
  max.size = 8000 * 8000,
  method = c("mult1", "mult", "indepMH"),
  verbose = FALSE,
  min.sample = 100,
  max.expected.size.before.switch = 1000 * 1e+07,
  init = NULL,
  burnin = floor(n/4),
  filter.every = 10,
  zero.prob.zero.samples = TRUE,
  allow.size.less.than.K = FALSE
)

sampleNMultilevelMultinomial(
  ndraws = ncol(pixel.index.mat),
  pixel.index.mat = NULL,
  urban.mat = NULL,
  area.mat = NULL,
  easpa.list,
  pop.mat,
  stratify.by.urban = TRUE,
  verbose = TRUE,
  return.EA.info = FALSE,
  min.HH.per.EA = 25,
  fix.HH.per.EA = NULL,
  fix.pop.per.HH = NULL
)

sampleNMultilevelMultinomialFixed(
  clusters.per.pixel,
  ndraws = ncol(pixel.indices),
  pixel.indices = NULL,
  urbanVals = NULL,
  areaVals = NULL,
  easpa,
  pop.mat,
  stratify.by.urban = TRUE,
  verbose = TRUE
)
```

## Arguments

- easpa:

  Census frame. See
  [`simPopCustom`](https://richardli.github.io/SUMMER/reference/simPop.md)
  for details

- pop.mat:

  data.frame of pixellated grid of population densities. See
  [`simPopCustom`](https://richardli.github.io/SUMMER/reference/simPop.md)
  for details

- level:

  Whether to calculate results at the integration grid or EA level

- pixel.index.mat:

  Matrix of pixel indices associated with each EA and draw. Not required
  by getExpectedNperEA unless level == "EA"

- i:

  Index

- urban:

  If TRUE, calculate only for urban part of the area. If FALSE, for only
  rural part

- stratify.by.urban:

  whether or not to stratify calculations by urban/rural classification

- validation.pixel.I:

  CURRENTLY FOR TESTING PURPOSES ONLY a set of indices of pixels for
  which we want to simulate populations (used for pixel level
  validation)

- n:

  Number of samples

- poppsub:

  Population per subarea. See
  [`simPopCustom`](https://richardli.github.io/SUMMER/reference/simPop.md)
  for details

- min1.per.subarea:

  Whether or not to ensure there is at least 1 EA per subarea. See
  [`simPopCustom`](https://richardli.github.io/SUMMER/reference/simPop.md)
  for details

- min.sample:

  The minimum number of samples per `chunk` of samples for truncated
  multinomial sampling. Defaults to 1

- method:

  If min1.per.subarea is TRUE, the sampling method for the truncated
  multinomial to use with rmulitnom1. rmultinom1 automatically switches
  between them depending on the number of expected samples. The methods
  are:

  mult1

  :   rejection sampling from multinomial plus 1 in each category

  mult

  :   rejection sampling from multinomial if any category has zero count

  indepMH

  :   independent Metropolis-Hastings using multinomial plus 1
      distribution as proposal

- easpsub:

  This could either be total EAs per subarea, or subarea crossed with
  urban or rural if stratify.by.urban is TRUE

- size:

  Multinomial size parameter. See
  [`rmultinom`](https://rdrr.io/r/stats/Multinom.html)

- prob:

  Multinomial probability vector parameter. See
  [`rmultinom`](https://rdrr.io/r/stats/Multinom.html)

- max.size:

  The maximum number of elements in a matrix drawn from the proposal
  distribution per sample chunk.

- verbose:

  Whether to print progress as the function proceeds

- max.expected.size.before.switch:

  Max expected number of samples / k, the number of categories, before
  switching method

- init:

  Initial sample if method is `indepMH`

- burnin:

  Number of initial samples before samples are collected if method is
  `indepMH`

- filter.every:

  Store only every filter.every samples if method is i`indepMH`

- zero.prob.zero.samples:

  If TRUE, set samples for parts of prob vector that are zero to zero.
  Otherwise they are set to one.

- allow.size.less.than.K:

  If TRUE, then if size \< the number of categories (k), returns matrix
  where each column is vector of size ones and k - size zeros. If FALSE,
  throws an error if size \< k

- ndraws:

  Number of draws

- urban.mat:

  Matrix of urbanicities associated with each EA and draw

- area.mat:

  Matrix of areas associated with each EA and draw

- easpa.list:

  A list of length n with each element being of the format of easpa
  giving the number of households and EAs per stratum. It is assumed
  that the number of EAs per stratum is the same in each list element.
  If easpa.list is a data frame, number of households per stratum is
  assumed constant

- return.EA.info:

  Whether a data frame at the EA level is desired

- min.HH.per.EA:

  The minimum number of households per EA (defaults to 25, since that is
  the number of households sampled per DHS cluster)

- fix.HH.per.EA:

  If not NULL, the fixed number of households per EA

- fix.pop.per.HH:

  If not NULL, the fixed target population per household

- clusters.per.pixel:

  CURRENTLY FOR TESTING PURPOSES ONLY a vector of length
  nIntegrationPoints specifying the number of clusters per pixel if they
  are fixed

- pixel.indices:

  A nEA x n matrix of pixel indices associated with each EA per
  simulation/draw

- urbanVals:

  A nEA x n matrix of urbanicities associated with each EA per
  simulation/draw

- areaVals:

  A nEA x n matrix of area names associated with each EA per
  simulation/draw

## Details

**\[experimental\]**

## Functions

- `getExpectedNperEA()`: Calculates expected denominator per enumeration
  area.

- `getSortIndices()`: For recombining separate multinomials into the
  draws over all grid points

- `rStratifiedMultnomial()`: Gives nIntegrationPoints x n matrix of
  draws from the stratified multinomial with values corresponding to the
  value of \|C^g\| for each pixel, g (the number of EAs/pixel)

- `rStratifiedMultnomialBySubarea()`: Gives nIntegrationPoints x n
  matrix of draws from the stratified multinomial with values

- `rMyMultinomial()`:

- `rMyMultinomialSubarea()`:

- `rmultinom1()`: Random (truncated) multinomial draws conditional on
  the number of each type being at least one

- `sampleNMultilevelMultinomial()`: Take multilevel multinomial draws
  first from joint distribution of number of households per EA given the
  total per stratum, and then from the joint distribution of the total
  target population per household given the total per stratum

- `sampleNMultilevelMultinomialFixed()`: Same as
  sampleNMultilevelMultinomial, except the number of EAs per pixel is
  fixed
