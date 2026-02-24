# Helper function of [`pixelPopToArea`](https://richardli.github.io/SUMMER/reference/aggPop.md)

Aggregates population from the pixel level to the level of the area of
interest.

## Usage

``` r
aggPixelPreds(
  Zg,
  Ng,
  areas,
  urban = target.pop.mat$urban,
  target.pop.mat = NULL,
  use.density = FALSE,
  stratify.by.urban = TRUE,
  normalize = use.density
)
```

## Arguments

- Zg:

  nIntegrationPoint x nsim matrix of simulated response (population
  numerators) for each pixel and sample

- Ng:

  nIntegrationPoint x nsim matrix of simulated counts (population
  denominators) for each pixel and sample

- areas:

  nIntegrationPoint length character vector of areas (or subareas)

- urban:

  nIntegrationPoint length vector of indicators specifying whether or
  not pixels are urban or rural

- target.pop.mat:

  same as in
  [`simPopCustom`](https://richardli.github.io/SUMMER/reference/simPop.md)

- use.density:

  whether to use population density as aggregation weights.

- stratify.by.urban:

  whether or not to stratify simulations by urban/rural classification

- normalize:

  if TRUE, pixel level aggregation weights within specified area are
  normalized to sum to 1. This produces an average of the values in Zg
  rather than a sum. In general, should only be set to TRUE for smooth
  integrals of risk.
