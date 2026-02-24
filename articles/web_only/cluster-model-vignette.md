# Specifying cluster-level model in SUMMER for mortality estimation

In this vignette, we will describe different ways to specify the
cluster-level model for U5MR estimation. We first load the package and
an example dataset `DemoData`, which contains model survey data provided
by DHS. Note that this data is simulated and does not represent any real
country’s data.

``` r
library(SUMMER)
library(dplyr)
data(DemoData)
```

For this simulated dataset, the strata variable is coded as region
crossed by urban/rural status. For our analysis with urban/rural
stratified model, we first construct a new strata variable that contains
only the urban/rural status, i.e., the additional stratification within
each region.

``` r
for (i in 1:length(DemoData)) {
    strata <- DemoData[[i]]$strata
    DemoData[[i]]$strata[grep("urban", strata)] <- "urban"
    DemoData[[i]]$strata[grep("rural", strata)] <- "rural"
}
counts.all <- NULL
for (i in 1:length(DemoData)) {
    vars <- c("clustid", "strata", "region", "time", "age")
    counts <- getCounts(DemoData[[i]][, c(vars, "died")], variables = "died", by = vars,
        drop = TRUE)
    counts <- counts %>%
        mutate(cluster = clustid, years = time, Y = died)
    counts$survey <- names(DemoData)[i]
    counts.all <- rbind(counts.all, counts)
}
periods <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14")
```

The `DemoData` includes $5$ surveys. We index surveys by $k = 1,...,5$.
The sampling frame that was used for survey $k$, will be denoted by
$r\lbrack k\rbrack$. For example, if the first $3$ surveys are based on
the first sampling frame and the next $2$ surveys are based on the
second sampling frame, then $r = \lbrack 1,1,1,2,2\rbrack$. When only
one survey is used, the $k$ index can be dropped in all formulas below.

We assume a discrete hazards model described in Mercer et al. (2015). We
use discrete time survival analysis to estimate age-specific monthly
probabilities of dying in user-defined age groups. We assume constant
hazards within the age bands. The default choice uses the monthly age
bands
$$\lbrack 0,1),\lbrack 1,12),\lbrack 12,24),\lbrack 24,36),\lbrack 36,48),\lbrack 48,60)$$
for U5MR and they can be easily specified by the user. The U5MR for area
$i$ and time $t$ can be calculated as,
$${\widehat{p}}_{it}^{\texttt{𝙷𝚃}} = {}_{60}{\widehat{q}}_{0}^{it} = 1 - \prod\limits_{a = 1}^{6}\left( 1 - {}_{n_{a}}{\widehat{q}}_{x_{a}}^{it} \right),$$
where $x_{a}$ and $n_{a}$ are the start and end of the $a$-th age group,
and ${}_{n_{a}}q_{x_{a}}^{it}$ is the probability of death in age group
$\lbrack x_{a},x_{a} + n_{a})$ in area $i$ and time $t$, with
${}_{n_{a}}{\widehat{q}}_{x_{a}}^{it}$ the estimate of this quantity.
This calculation follows the synthetic cohort life table approach in
which mortality probabilities for each age segments based on real cohort
mortality experience are combined.

For cluster-level modeling, we consider a beta-binomial model for the
probability (hazard) of death from month $m$ to $m + 1$ in survey $k$
and at cluster $c$ in year $t$. This model allows for overdispersion
relative to the binomial model. Assuming constant hazards within age
bands, we assume the number of deaths occurring within age band
$a\lbrack m\rbrack$, in cluster $c$, time $t$, and survey $k$ follow the
beta-binomial distribution,
$$Y_{a{\lbrack m\rbrack},k,c,t}\ |\ p_{a{\lbrack m\rbrack},k,c,t} \sim \text{BetaBinomial}\left( \ n_{a{\lbrack m\rbrack},k,c,t}\ ,\ p_{m,k,c,t}\ ,d\  \right),$$
where $p_{m,k,c,t}$ is the monthly hazard at $m$-th month of age, in
cluster $c$, time $t$, and survey $k$ and $d$ is the overdispersion
parameter.

In the likelihood model, $n_{a{\lbrack m\rbrack},k,c,t}$ is a known
constant. The default prior for $d$ is
$$\text{logit}(d) \sim N\left( 0,\sqrt{1/0.4} \right)$$ The mean and
precision of the normal prior can be modified by specifying the
`overdisp.mean` and `overdisp.prec` arguments.

## Specifying the latent logistic model

For simplicity of notation, we will consider within-region
stratification by urban/rural status in this vignette. In SUMMER, more
than two levels of stratification is automatically allowed when `strata`
variable in the input data frame contains more unique values. Under the
urban/rural stratification, the most general form of the latent logistic
model we use is, $$\begin{aligned}
{p_{m,k,c,t} =} & {\text{expit}\left( \alpha_{m,c,k,t} + \epsilon_{t} + b_{k} \right),} \\
{\alpha_{m,k,c,t} =} & {\beta_{a{\lbrack m\rbrack},r{\lbrack k\rbrack},t}I\left( s_{c} \in {\mspace{6mu}\text{rural}\mspace{6mu}} \right) + \gamma_{a{\lbrack m\rbrack},r{\lbrack k\rbrack},t}I\left( s_{c} \in {\mspace{6mu}\text{urban}\mspace{6mu}} \right)} \\
 & {\; + S_{i{\lbrack s_{c}\rbrack}} + e_{i{\lbrack s_{c}\rbrack}} + \delta_{i{\lbrack s_{c}\rbrack},t} + \text{BIAS}_{k,t}.}
\end{aligned}$$

We will break down the terms next.

### Sharing information across age groups

In the temporal model, we usually would like to share some information
across different age groups. This requires specifying a reduced group of
age bands for some component. We use $a^{*}\lbrack m\rbrack$ here to
differentiate from the age bands used in the hazard likelihood. The
default choice for U5MR in the package is
$$a^{*}\lbrack m\rbrack = \begin{cases}
1 & {{\text{if}\mspace{6mu}}m = 0,} \\
2 & {{\text{if}\mspace{6mu}}m = 1,\ldots,11,} \\
3 & {{\text{if}\mspace{6mu}}m = 12,\ldots,59.}
\end{cases}$$ The $a\lbrack m\rbrack$ are specified in `age.group`
argument, and the size of each age group is specified in `age.n`
argument. The reduced mapping $a^{*}\lbrack m\rbrack$ is specified in
`age.time.group` argument.

### Survey fixed effect

We include a survey fixed effect $b_{k}$ with the constraint
$\sum_{k}b_{k}\mathbf{1}_{r{\lbrack k\rbrack} = r} = 0$ for each
sampling frame $r$, so that the main temporal trends are identifiable
for each sampling frame. The $b_{k}$ terms are not included in the
prediction, i.e., they are set to zero. If we ignore any survey-specific
difference, we can remove the $b_{k}$ term by setting
`survey.effect = FALSE`.

### Unstructured temporal effect

The $\epsilon_{t}$ are unstructured temporal effects that allow for
perturbations over time. It is a contextual choice whether they are used
in predictions. By default, when applying
[`getSmoothed()`](https://richardli.github.io/SUMMER/reference/getSmoothed.md)
to a fitted cluster-level model, $\epsilon_{t}$ are not included in the
predictions. They can be included by specifying
`include.time.unstruct = TRUE`.

### Structured temporal effect

The temporal main effects are specified for urban and rural stratum
seperately with $\beta_{a^{*}{\lbrack m\rbrack},r{\lbrack k\rbrack},t}$
and $\gamma_{a^{*}{\lbrack m\rbrack},r{\lbrack k\rbrack},t}$. These
effects are also specific to each sampling frame $r$, as the definition
of urban and rural usually change across different sampling frames.

We do not share any information across frames currently in SUMMER, thus
all terms specific to different sampling frames are modeled as
independent a priori.

The age-frame-stratum-specific temporal effects can be modelled by a
first order random walk (RW1), second order frandom walk (RW2), or first
order autoregressive process (AR1). In all cases, we apply a sum-to-zero
constraint on the random walk or autoregressive process and include an
explicit intercept for the full age groups $a$ (instead of $a^{*}$. For
example, $$\beta_{a,r,t} = \beta_{a,r,0} + \beta_{a^{*},r,t}^{*}$$ where
$\sum_{t}\beta_{a,r,t}^{*} = 0$. For RW1 and AR1 model, we include an
additional linear trend by default, so that
$$\beta_{a,r,t} = \beta_{a,r,0} + {\widetilde{\beta}}_{a^{*},r}t + \beta_{a^{*},r,t}^{*}$$
We will discuss the modeling of the linear trend later. First, we start
with the default RW2 model.

#### No urban/rural effect

In the simplest case, we can let $\beta_{a,r,t} = \gamma_{a,r,t}$. This
can be achieved by setting `strata` variable of the input data frame to
be all NA.

``` r
counts.all.no.strat <- counts.all
counts.all.no.strat$strata <- NA
fit1 <- smoothCluster(data = counts.all.no.strat, Amat = DemoMap$Amat, family = "betabinomial",
    year.label = c(periods, "15-19"), time.model = "rw2", survey.effect = TRUE, strata.time.effect = FALSE)
rownames(fit1$fit$summary.fixed)
```

    ## [1] "age.intercept0"     "age.intercept1-11"  "age.intercept12-23"
    ## [4] "age.intercept24-35" "age.intercept36-47" "age.intercept48-59"

We can see the fixed effects include only the six elements
$\beta_{1,r,0},...,\beta_{6,r,0}$.

#### Time-invariant urban/rural effect

We can also specify the random effects
$\beta_{a^{*},r,t}^{*} = \gamma_{a^{*},r,t}^{*} + \Delta$. This

This can be specified by setting `strata.time.effect = FALSE`.

``` r
fit2 <- smoothCluster(data = counts.all, Amat = DemoMap$Amat, family = "betabinomial",
    year.label = c(periods, "15-19"), survey.effect = TRUE, strata.time.effect = FALSE)
rownames(fit2$fit$summary.fixed)
```

    ##  [1] "age.intercept0:rural"     "age.intercept1-11:rural" 
    ##  [3] "age.intercept12-23:rural" "age.intercept24-35:rural"
    ##  [5] "age.intercept36-47:rural" "age.intercept48-59:rural"
    ##  [7] "age.intercept0:urban"     "age.intercept1-11:urban" 
    ##  [9] "age.intercept12-23:urban" "age.intercept24-35:urban"
    ## [11] "age.intercept36-47:urban" "age.intercept48-59:urban"
    ## [13] "strataurban"

Now the fixed effects include the six elements
$\beta_{1,r,0},...,\beta_{6,r,0}$ for the urban strata, six elements
$\gamma_{1,r,0},...,\gamma{6,r,0}$ for the rural strata and a constant
urban/rural effect for the random effect difference.

#### Time-varying urban/rural effect

A more flexible and usually reasonable model is to assume each stratum
have indepent age-specific trends, i.e., $\beta_{a^{*},r,t}^{*}$ and
$\gamma_{a^{*},r,t}^{*}$ each follow its own trends.

``` r
fit3 <- smoothCluster(data = counts.all, Amat = DemoMap$Amat, family = "betabinomial",
    year.label = c(periods, "15-19"), survey.effect = TRUE, strata.time.effect = TRUE)
rownames(fit3$fit$summary.fixed)
```

    ##  [1] "age.intercept0:rural"     "age.intercept1-11:rural" 
    ##  [3] "age.intercept12-23:rural" "age.intercept24-35:rural"
    ##  [5] "age.intercept36-47:rural" "age.intercept48-59:rural"
    ##  [7] "age.intercept0:urban"     "age.intercept1-11:urban" 
    ##  [9] "age.intercept12-23:urban" "age.intercept24-35:urban"
    ## [11] "age.intercept36-47:urban" "age.intercept48-59:urban"

We can check now that there are in total $42$ structured temporal random
effect terms, corresponding to $3$ urban time series and $3$ rural time
series.

``` r
nrow(fit3$fit$summary.random$time.struct)
```

    ## [1] 42

#### Linear trend

As mentioned before, for RW1 and AR1 model, by default we include an
additional linear trend by default, so that
$$\beta_{a,r,t} = \beta_{a,r,0} + {\widetilde{\beta}}_{a^{*},r}X_{t} + \beta_{a^{*},r,t}^{*}$$
The linear trend is modelled as a fixed effect for a rescaled time index
variable $X_{t} = (t - T/2)/(T - 1)$. That is, the coefficient
$\widetilde{\beta}$ should be interpreted as unit change on the logit
scale comparing the last to the first time period.

``` r
fit4 <- smoothCluster(data = counts.all, Amat = DemoMap$Amat, family = "betabinomial",
    year.label = c(periods, "15-19"), time.model = "ar1", survey.effect = TRUE, strata.time.effect = TRUE)
rownames(fit4$fit$summary.fixed)
```

    ##  [1] "time.slope.group1"        "time.slope.group2"       
    ##  [3] "time.slope.group3"        "time.slope.group4"       
    ##  [5] "time.slope.group5"        "time.slope.group6"       
    ##  [7] "age.intercept0:rural"     "age.intercept1-11:rural" 
    ##  [9] "age.intercept12-23:rural" "age.intercept24-35:rural"
    ## [11] "age.intercept36-47:rural" "age.intercept48-59:rural"
    ## [13] "age.intercept0:urban"     "age.intercept1-11:urban" 
    ## [15] "age.intercept12-23:urban" "age.intercept24-35:urban"
    ## [17] "age.intercept36-47:urban" "age.intercept48-59:urban"

The $6$ additional linear trends are included in the fixed effect. To
see which slope correspond to which age-stratum pairs, we can inspect
the `slope.fixed.output` field in the returned object.

``` r
fit4$slope.fixed.output
```

    ## $time.slope.group1
    ## [1] "0:rural"
    ## 
    ## $time.slope.group2
    ## [1] "1-11:rural"
    ## 
    ## $time.slope.group3
    ## [1] "12-23:rural" "24-35:rural" "36-47:rural" "48-59:rural"
    ## 
    ## $time.slope.group4
    ## [1] "0:urban"
    ## 
    ## $time.slope.group5
    ## [1] "1-11:urban"
    ## 
    ## $time.slope.group6
    ## [1] "12-23:urban" "24-35:urban" "36-47:urban" "48-59:urban"

The linear trend term can also be set to $0$ by setting
`linear.trend = FALSE`

``` r
fit5 <- smoothCluster(data = counts.all, Amat = DemoMap$Amat, family = "betabinomial",
    year.label = c(periods, "15-19"), time.model = "ar1", survey.effect = TRUE, strata.time.effect = TRUE,
    linear.trend = FALSE)
rownames(fit5$fit$summary.fixed)
```

    ##  [1] "age.intercept0:rural"     "age.intercept1-11:rural" 
    ##  [3] "age.intercept12-23:rural" "age.intercept24-35:rural"
    ##  [5] "age.intercept36-47:rural" "age.intercept48-59:rural"
    ##  [7] "age.intercept0:urban"     "age.intercept1-11:urban" 
    ##  [9] "age.intercept12-23:urban" "age.intercept24-35:urban"
    ## [11] "age.intercept36-47:urban" "age.intercept48-59:urban"

#### Shared Linear trend

We can also enforce shared linear trends across all age groups, though
this is usually too strong an assumption, i.e.,
$$\beta_{a,r,t} = \beta_{a,r,0} + {\widetilde{\beta}}_{r}t + \beta_{a^{*},r,t}^{*}$$

``` r
fit6 <- smoothCluster(data = counts.all, Amat = DemoMap$Amat, family = "betabinomial",
    year.label = c(periods, "15-19"), time.model = "ar1", survey.effect = TRUE, strata.time.effect = TRUE,
    linear.trend = TRUE, common.trend = TRUE)
rownames(fit6$fit$summary.fixed)
```

    ##  [1] "time.slope"               "age.intercept0:rural"    
    ##  [3] "age.intercept1-11:rural"  "age.intercept12-23:rural"
    ##  [5] "age.intercept24-35:rural" "age.intercept36-47:rural"
    ##  [7] "age.intercept48-59:rural" "age.intercept0:urban"    
    ##  [9] "age.intercept1-11:urban"  "age.intercept12-23:urban"
    ## [11] "age.intercept24-35:urban" "age.intercept36-47:urban"
    ## [13] "age.intercept48-59:urban"

The $6$ additional linear trends are included in the fixed effect. To
see which slope correspond to which age-stratum pairs, we can inspect
the `slope.fixed.output` field in the returned object.

``` r
fit6$slope.fixed.output
```

    ## $time.slope.group1
    ## [1] "0:rural"
    ## 
    ## $time.slope.group2
    ## [1] "1-11:rural"
    ## 
    ## $time.slope.group3
    ## [1] "12-23:rural" "24-35:rural" "36-47:rural" "48-59:rural"
    ## 
    ## $time.slope.group4
    ## [1] "0:urban"
    ## 
    ## $time.slope.group5
    ## [1] "1-11:urban"
    ## 
    ## $time.slope.group6
    ## [1] "12-23:urban" "24-35:urban" "36-47:urban" "48-59:urban"

### Spatial effects

The spatial effects are assumed to be shared across all age groups and
strata. The $S_{i{\lbrack s_{c}\rbrack}} + e_{i{\lbrack s_{c}\rbrack}}$
terms correspond to the sum of a structured spatial effect and an
unstructured effect for each area. It is parameterized with a BYM2 model
described in Riebler et al. (2016).

### Space-time interaction effects

The space-time interaction term $\delta_{it}$ are modelled with the type
I to IV interactions of the chosen temporal model and the ICAR or
independent model in space. The four types of models are described in
Knorr-Held (2000) . The specification of type I to IV interactions is
through the argument `type.st`.

We can further include random slopes in the space-time interaction
terms, i.e., $$\delta_{it} = \rho_{i}X_{t} + \delta_{it}^{*}$$ where
$\delta^{*}$ is the usual type I to IV interactions and the time
covariate $X_{t}$ is rescaled to be -0.5 to 0.5, so that the random
slope can be interpreted as the total deviation from the main trend from
the first year to the last year to be projected, on the logit scale.

The random slopes are included when user specifies their hyperprior
`pc.st.slope.alpha` and`pc.st.slope.u`. For example, the following model
specifies a type IV interaction between AR1 in time and ICAR in space,
with region-specific random slope in time, with PC prior such that
$$Prob\left( \left| \rho_{i} \right| > 2 \right) = 0.1$$ The main
temporal random effect is set to be a second order random walk.

``` r
fit6 <- smoothCluster(data = counts.all, Amat = DemoMap$Amat, family = "betabinomial",
    year.label = c(periods, "15-19"), survey.effect = TRUE, strata.time.effect = TRUE,
    linear.trend = FALSE, time.model = "rw2", st.time.model = "ar1", pc.st.slope.u = 2,
    pc.st.slope.alpha = 0.1)
```

The random slopes can be extracted from

``` r
fit6$fit$summary.random$st.slope.id[, c(1:3)]
```

    ##   ID   mean   sd
    ## 1  1  0.051 0.30
    ## 2  2  0.113 0.30
    ## 3  3  0.120 0.31
    ## 4  4 -0.284 0.33

### Hyperpriors

By default, all hyperpriors are assumed to be PC priors. The hyperprior
arguments start with `pc.`. Their definition can be easily found in the
help file.

### References

Knorr-Held, Leonhard. 2000. “Bayesian Modelling of Inseparable
Space-Time Variation in Disease Risk.” *Statistics in Medicine* 19:
2555–67.

Mercer, Laina D, Jon Wakefield, Athena Pantazis, Angelina M Lutambi,
Honorati Masanja, and Samuel Clark. 2015. “Small Area Estimation of
Childhood Mortality in the Absence of Vital Registration.” *Annals of
Applied Statistics* 9: 1889–1905.

Riebler, Andrea, Sigrunn H Sørbye, Daniel Simpson, and Håvard Rue. 2016.
“An Intuitive Bayesian Spatial Model for Disease Mapping That Accounts
for Scaling.” *Statistical Methods in Medical Research* 25: 1145–65.
