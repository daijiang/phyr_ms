---
title: "`phyr`: An R package for phylogenetic species-distribution modelling in ecological communities"
author: "Daijiang Li, Russell Dinnage, Lucas Nell, Anthony Ives"
date: "20 September, 2019"
output:
  bookdown::html_document2:
    number_sections: no
    theme: sandstone
    toc: yes
    css: style.css
    keep_md: yes
  bookdown::word_document2: 
    toc: no
  bookdown::tufte_html2:
    number_sections: no
    toc: yes
fontsize: 12pt
link-citations: yes
always_allow_html: yes
csl: https://raw.githubusercontent.com/citation-style-language/styles/master/methods-in-ecology-and-evolution.csl
bibliography: ref.bib
---



Running Title: Model-based phylogenetic analyses

*Summary*

1. Phylogenetic relationships among species are both a challenge and an opportunity for community ecology and comparative biology. Model-based phylogenetic analysis has advantages ... Currently, R functions to conduct these analyses were distributed across different packages and were mostly written in R, which can be very slow with large datasets.
2. To help remedy this situation, we created `phyr`, an R package that collects and updates functions to conduct model-based phylogenetic analyses. By collecting these functions into one package and rewrite core functions in c++, this package can facilitate the use of model-based phylogenetic analyses.
3. This paper presents the main functions of the package and provides an example using simulated data.
4. The `phyr` package provides a unified software environment with improved performance to conduct model-based phylogenetic analyses in Ecology. 

_Keywords_: phylogenetic diversity, phylogenetic generalized linear mixed models, functional trait, trait correlation, model-based methods


# Introduction

Ecological communities are collections of species that occur within the same geographical area. Which species occur within communities depends on the dispersal ability of species to enter the community, the environmental conditions that they find there, and the interactions that they have with other species in the community. These three processes -- dispersal, environmental tolerance, and species interactions -- depend on the traits that species possess and hence reflect evolutionary history and biogeographic processes [@warren2014mistaking; @gerhold2018deep]. For example, an insect species might only occur in a lake if its adult stage has long-distance flight capabilities, if it can tolerate the low pH of the lake, and if it can avoid the predators that are common. Because traits play a central role in the composition of species that make up a community, community composition will likely reflect, at least in part, phylogenetic relationships among species. For example, two closely related insects might have similar dispersal capability, pH tolerance, and predator avoidance behavior, making them more likely to occur in the same lake. The recognition that phylogenetic relationships can increase our understanding of communities has led to a growing number of statistical methods for analyzing phylogenetic community composition.

Just as the distributions of two species might reflect their proximity on a phylogenetic tree, the species occurring at two sites might reflect the sites’ geographical proximity. The most immediate possible cause of spatial correlations in species distributions is dispersal, if nearby sites are more likely to be colonized by a species. Spatial proximity may also be a surrogate for environmental variables that are unknown or unmeasured. For example, an insect species might occur in two nearby lakes because they both have low pH, yet pH has not been measured. Just as phylogenetic relationships among species can generate correlations between species in which sites they occupy, so too can spatial proximity generate correlations between sites in the species they contain.

How species respond to environmental factors, and how they respond to each other, depend on their traits. Therefore, relationships among functional traits can provide insights about the evolutionary history that has shaped species trait so they can occupy the same sites. For example, two insect species that occur in the same lake might share both long-range flight abilities and tolerance to low pH. Is the positive correlation between these two traits caused by correlated selective forces? A challenge to answering this question is that phylogenetic correlations between trait values of closely related species can give the appearance of positive correlations between traits: two species might have both long-range flight abilities and tolerance to low pH only because they are phylogenetically closely related. Therefore, it is important to account for species’ evolutionary history when studying correlations among functional traits. 

Statistical models for phylogenetic community composition provide flexible tools for exploring the many possible factors underlying the distribution of species and the composition of communities. The models can describe complex relationships in the data, such as how phylogenetically related species might respond similarly to the same environmental gradient, or how phylogenetically related species might exclude each other from the same communities. They also give a firm statistical basis to test these patterns, the ability to simulate data sets from the fitted model, and the ability to predict the composition of unsurveyed communities. These benefits of phylogenetic community composition models come with costs: building models can be intricate and fitting them computationally slow. The R package `phyr` is designed to overcome many of the costs with a user-friendly interface, flexibility to build a rich collection of models, and good computational performance. Below, we first give a brief overview of the structure and syntax of two key functions `pglmm()` and `cor_phylo()`. We then compare them to methods and programs that are currently available. Finally, we give applications of `pglmm()` and `cor_phylo()` to simulated data to illustrate their implementation and output.

# Overview of `phyr` 

`Phyr` contains three groups of functions (Table \@ref(tab:mainFunc)): phylogenetic generalized linear mixed models (`pglmm()`), phylogenetic correlations between functional traits (`cor_phylo()`), and phylogenetic diversity metrics (e.g. `psv()`, `pse()`). The workhorse functions of all groups are written in C++ to increase computational speed. Here, we will focus on the first two groups of functions because they are more complicated than phylogenetic diversity metrics. 

(ref:maunFuncCaption) List of main functions in `phyr` package.

<table class="table table-striped table-hover" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:mainFunc)(ref:maunFuncCaption)</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Group </th>
   <th style="text-align:left;"> Main Functions </th>
   <th style="text-align:left;"> Brief Description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> A </td>
   <td style="text-align:left;width: 25%; "> `psv`; `pse`; `psr`; `psc`; `psd` </td>
   <td style="text-align:left;"> Phylogenetic alpha diversity of communities </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;width: 25%; "> `pcd` </td>
   <td style="text-align:left;"> Pairwise phylogenetic beta diversity of communities </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;width: 25%; "> `vcv2` </td>
   <td style="text-align:left;"> Convert a phylogeny to a var-cov matrix, a faster version of `ape::vcv` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> B </td>
   <td style="text-align:left;width: 25%; "> `cor_phylo` </td>
   <td style="text-align:left;"> Correlations among multiple traits with phylogenetic signal </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;width: 25%; "> `binaryPGLMM` </td>
   <td style="text-align:left;"> Phylogenetic Generalized Linear Mixed Model for binary data; each species only allowed to have one value </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C </td>
   <td style="text-align:left;width: 25%; "> `pglmm` </td>
   <td style="text-align:left;"> Phylogenetic Generalised Linear Mixed Model for ecological community data (e.g., species composition across sites; bipartite interactions) </td>
  </tr>
</tbody>
</table>

## _`pglmm()`_

`pglmm()` constructs and fits generalized linear mixed models that incorporate covariance matrices containing the phylogenetic relationships among species. The syntax for `pglmm()` resembles that used in the R package `lme4`, and indeed `pglmm()` will fit most of the models that can be fit with `lmer()` and `glmer()`. `pglmm()` goes beyond `lmer()`and `glmer()` by allowing the specification of covariance matrices, which could be a phylogenetic covariance matrix or any other covariance matrix that the user defines (e.g., spatial or temporal autocorrelation matrix). `pglmm()` can also fit models with “nested” covariance structures (e.g. species phylogenetic covariance matrix nested within site covariance matrix). `pglmm()` can operate in both frequentist mode, with the distribution of species among communities being Gaussian, Binary, Binomial or Poisson, and Bayesian mode with the addition of zero-inflated Binomial and Poisson distributions. Finally, it is our hope that the formula syntax of `pglmm()` can be used to fit similar models with other programs such as Stan (e.g. via R package `brm`).

A general example of the syntax for `pglmm()` is


```r
pglmm(
  Y ~ trait * env +
    (1 | sp__) +
    (1 | site__) +
    (trait | site) +
    (env | sp__) +
    (1 | sp__@site),
  data = data,
  cov_ranef = list(sp = phy.sp, site = V.space),
  family = 'binomial',
  bayes = FALSE,
  REML = TRUE
) 
```

Here, Y is a binary (Bernoulli) dependent variable which takes values of either 0 or 1. The specification `family = 'binomial'` allows binary data and also binomial data for which Y is a matrix containing columns for successes and failures. The independent variables `trait` and `env` take on different values for each species and site, respectively. Sites and species (`sp`) are treated as random effects: thus, `(1|site)` implies that a value from a Gaussian random variable is picked for each site, thereby representing unmeasured differences among sites. For the case of species, the double underscore in `(1|sp__)` implies that, in addition to a random effect for species, there is a second random effect which contains the phylogenetic relationships among species. The phylogenetic random effect assumes that values for each species are picked from a multivariate Gaussian distribution with phylogenetic covariance matrix 𝚺. A covariance matrix 𝚺 is specified by `cov_ranef = list(sp = phy.sp, site = V.space)`. The covariance matrix `phy.sp` associated with species can be a `phylo` object from the R package `ape`. It is also possible to specify an explicit covariance matrix, such as `site = V.space`. To construct 𝚺 from a 'phylo' object, `pglmm()` assumes that the residual variation associated with species follows a Brownian motion model of evolution, so that the covariance between two species is proportional to the distance between them on the phylogenetic tree. 

The syntax `(1|sp__)` or  `(1|site__)` generates two random effects, one without and one with phylogenetic or spatial covariances; in contrast, `(1|sp)` would generate only a single random effect that is independent among species. `pglmm()` forces in a term for `(1|sp)` whenever `(1|sp__)` is specified, because otherwise any difference among species would be captured by the diagonal elements in 𝚺 even in the absence of covariances among phylogenetically related species which are specified by the off-diagonal elements of 𝚺. Therefore, if `(1|sp)` were not included, this could lead to the identification of phylogenetic signal in the abundances of species even in its absence. To account for differences among sites in how they select for species with different traits, `(trait|site)` allows the slope of Y against `trait` to be a Gaussian random variable. Similarly, to account for the differences among species for how they respond to `env`, `(env|sp__)` allows the relationship of Y against `env` to be given by two slopes, the first slope that is picked from a Gaussian random variable in which species are independent and the second slope that is picked from a multivariate Gaussian with covariance matrix 𝚺. Finally, `(1|sp__@site)` generates a nested term: within a site, the residual variation in Y shows phylogenetic relatedness, with phylogenetically related species more likely to occur in the same site. Note that `(1|sp__)` differs from `(1|sp__@site)` because `(1|sp__)` generates differences in the mean value of Y for species across all sites, whereas `(1|sp__@site)` is local, giving the covariances among species only within sites. Other forms of a nested term are available in `pglmm()`, which can be used to study more complicated questions such as bipartite networks.

With `bayes = FALSE`, `pglmm()` is fit using a frequentist approach. ML or REML is used for fitting, with `REML = TRUE` as the default. For a non-Gaussian model (e.g., `family = 'binomial'`), an iterated quasi-likelihood method is used for model fitting which gives the approximate likelihood; p-values for the fixed effects are given by a Wald test and for the random effects by profile likelihood, although we recommend bootstrap-based tests when computationally feasible. Note that `REML = TRUE` is an option for non-Gaussian models (in contrast to `glmer`) due to the algorithm used. With `bayes = TRUE`, a Bayesian approach is implemented using INLA [@rue2009approximate], which gives parameter estimates and credible intervals. For large problems with the number of species-site combinations exceeding 2000, the Bayesian computations are considerably faster than the frequentist computations. Finally, a key to interpreting the results from a model is understanding the structure of the covariance matrices associated with the random effects. Therefore, `pglmm()` has associated plotting functions `pglmm.plot.ranef()` that present the design matrices for the random effects (Fig. \@ref(fig:designPlot)). 

(ref:designPlotCap) The structures of design matrices of random terms in a phylogenetic generalized linear mixed model with 5 species and 10 sites.


```r
knitr::include_graphics(normalize_path("designPlot.pdf"))
```

<div class="figure" style="text-align: center">
<img src="/Users/dli/Dropbox/UFL/phyr_ms/designPlot.pdf" alt="(ref:designPlotCap)" width="100%" />
<p class="caption">(\#fig:designPlot)(ref:designPlotCap)</p>
</div>

## *`cor_phylo()`*
 
`cor_phylo()` makes it possible to compare suites of traits among species, accounting for their phylogenetic relatedness. To identify suites of traits under joint selection, such as traits that together make up adaptive syndromes, it is necessary to perform a correlation analysis in which phylogenetic relatedness is factored out. `cor_phylo()` does this. It can also include within-species variation (e.g., measurement error) which should better-expose the underlying correlations in traits among species. Whereas `pglmm()` can be used to identify the composition of communities within a region, `cor_phylo()` can be used to assess patterns of traits among species that make up the regional species pool.

The syntax for `cor_phylo()` is


```r
cor_phylo(
  variates = ~ trait1 + trait2,
  species = ~ sp,
  phy = phy.sp,
  covariates = list(trait1 ~ env),
  meas_errors = list(trait1 ~ me1, trait2 ~ me2),
  data = data,
  boot = 2000
)
```

In this example, the correlation between two traits, `trait1` and `trait2`, is assessed, and the column named `sp` in `data` identifies the species. The object `phy.sp` specifies the phylogenetic covariance matrix as a 'phylo' object from the `ape` package. `cor_phylo()` estimates the phylogenetic signal for each trait by assuming that trait evolution is given by a Ornstein-Uhlenbeck process. The term `covariates = list(trait1 ~ env)` includes the independent variable `env` for `trait1`, to remove possible confounding effects; only an intercept is estimated if no covariate is provided for a trait. Within-species variation (measurement error) is specified by `meas_errors = list(trait1 ~ me1, trait2 ~ me2)`, where `me1` and `me2` are the standard errors for `trait1` and  `trait2`, respectively, of values at the tips of the phylogenetic tree. If within-species standard errors are not provided for a given trait, the trait values are assumed to be known without error. Finally, `cor_phylo()` can perform parametric bootstrapping to give confidence intervals for all parameter estimates: correlations, phylogenetic signals, covariate coefficients, and coefficient covariances.

# Relationships to other methods and software

`pglmm()` and `cor_phylo()` both appear in existing R packages (`pez` and `ape`, respectively), although the versions in `phyr` represent considerable improvements in ease-of-use, computational speed, and flexibility. Both have new syntax that makes them more intuitive to use. `pglmm()` also has new associated functions that plot the design of the covariance matrices (Fig. \@ref(fig:designPlot)), making model interpretation easier. Both are now coded in C++ (for key functions), which speeds computation time by 5-10X. `pglmm()` now supports several non-Gaussian distributions and allows Bayesian analyses using INLA that is particularly useful for large datasets. Finally, both include more output; for example, both now include facilities to perform likelihood ratio tests and compute AIC and BIC values for model comparisons.

## *`pglmm()`*

`pglmm()` is syntactically modeled after `lmer()` and `glmer()` in `lme4`, although it allows the specification of phylogenetic covariance matrices. `pglmm()` also allows "nested" models which arise when phylogenetic covariances only act within single communities, rather than among communities. Such nested models make it possible to assess whether phylogenetic relatedness affects the abundance of species within the same communities, such as whether competition between closely related species excludes one of the competitors from communities where the other is present. Nested models are structurally incompatible with the architecture of `lme4`.

There are alternative programs to `pglmm()`, although they have limitations that `pglmm()` overcomes. @hadfield2013tale use the R package `MCMCglmm` [@mcmcglmm2010] to perform phylogenetic community analyses, although they also use ASReml because its penalized quasi-likelihood (PQL) approach is computationally much faster. Hierarchical Modelling of Species Communities (HMSC-R) [@tikhonov2019joint] performs community analyses using Bayesian MCMC approaches, although it does not include nested terms. It is also possible to code up specific phylogenetic community models using flexible Bayesian platforms such as WinBugs, Stan, and JAGS, although this will involve considerable programming and expertise.

### pglmm as a Joint Species Distribution Model (JSDM)

pglmm is a joint species distribution model where the (residual) dependencies among species are modelled as a simple function of phylogeny [@wilkinson2019comparison]. 


## *`cor_phylo()`*

The R package `mvMORPH` [@clavel2015mvmorph] can fit a broad range of models, of which `cor_phylo()` can be formulated as a special case. While `cor_phylo()` does not have the flexibility of `mvMORPH`, it is correspondingly simpler to use. Also, `cor_phylo()` has built-in bootstrapping capabilities that are necessary to give confidence in the parameter estimates and p-values. The function `evolvcv.lite()` in the R package `phytools` [@revell2012phytools] will compute phylogenetic correlations, and changes in phylogenetic correlations through time [see also @caetano2018estimating], although the phylogenetic covariance matrix is derived under the assumption of Brownian motion evolution. This contrasts `cor_phylo()` in which the strength of phylogenetic signal is computed at the same time as the correlation. It is also possible to code the `cor_phylo()` model using platforms such as WinBugs, Stan, and JAGS; but again, this will require considerable programming and expertise.   

# Example usage

We simulated datasets to demonstrate how to use `pglmm()` and `cor_phylo()`. Details about simulations can be found in the Appendix. Our goal of this section is to provide some general ideas about the input and output of these two functions instead of testing their statistical performances or interpreting the ecological meanings of model results. For those purposes, please see the package vignettes and @ives2018book.

## _`pglmm()`_

We fitted a PGLMM that examined how a hypothetical functional trait, environmental gradient, and their interaction affect distributions of 5 species across 10 sites (we focused on abundance and used the default family of data distribution [Gaussian] but other distributions can also be specified by setting the `family` argument) while accounting for species’ phylogenetic relationships and site spatial autocorrelations (by setting `cov_ranef = list(species = phy, site = Vspace)` where `species` and `site` are group variables of random terms, `phy` can be a phylogeny with class `phylo` or a phylogenetic variance-covariance matrix, `Vspace` is a matrix of spatial distances among sites). This model can also be fitted with Bayesian framework by setting `bayes = TRUE`.


```r
z <- pglmm(
  abund ~ 1 + envi + trait + envi:trait +
    (1 | species__)  + (1 | site__) + 
    (envi | species__) + (trait | site) +
    (1 | species__@site),
  data = dat, 
  cov_ranef = list(species = phy, site = Vspace)
)
summary(z)
```

```
## Linear mixed model fit by restricted maximum likelihood
## 
## Call:abund ~ 1 + envi + trait + envi:trait
## 
## logLik    AIC    BIC 
##  -75.4  176.9  186.8 
## 
## Random effects:
##                  Variance Std.Dev
## 1|species        2.61e-01 0.51097
## 1|species__      5.15e-05 0.00718
## 1|site           1.15e-01 0.33946
## 1|site__         1.11e-01 0.33363
## envi|species     6.01e-02 0.24508
## envi|species__   1.95e-05 0.00442
## trait|site       3.96e+01 6.29179
## 1|species__@site 1.59e-01 0.39882
## residual         9.32e-01 0.96531
## 
## Fixed effects:
##              Value Std.Error Zscore  Pvalue    
## (Intercept)  1.427     0.412   3.46 0.00053 ***
## envi         0.982     0.298   3.29 0.00099 ***
## trait       -9.260     4.059  -2.28 0.02254 *  
## envi:trait  -3.323     3.246  -1.02 0.30589    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The summary of model results includes the model fitting method (maximum likelihood or bayesian), the model formula, log likelihood and other statistics based on it (e.g., AIC, BIC, DIC), estimates of variances of random terms, coefficients of fixed terms and their uncertainties (Standard errors, Z scores, and P values if fitted by maximum likelihood or lower and upper boundaries of credible intervals).

## *`cor_phylo()`*

Here, we simulated two hypothetical functional traits (`trait_1` and `trait_2`) of 50 species. We set the true correlation between these two traits to be 0.7 and their phylogenetic signals to be 0.3 and 0.95, respectively.




```r
z2 <- cor_phylo(variates = ~ trait_1 + trait_2, 
                covariates = list(trait_2 ~ cov_trait_2),
                species = ~ sp, phy = phy, 
                meas_errors = list(trait_1 ~ se_trait_1, trait_2 ~ se_trait_2), 
                data = traits)
z2
```

```
## 
## Call to cor_phylo:
## cor_phylo(variates = ~trait_1 + trait_2, species = ~sp, phy = phy, covariates = list(trait_2 ~ cov_trait_2), meas_errors = list(trait_1 ~ se_trait_1, trait_2 ~ se_trait_2), data = traits) 
## 
## logLik    AIC    BIC 
##  -33.0   82.0   88.1 
## 
## Correlation matrix:
##         trait_1 trait_2
## trait_1   1.000   0.848
## trait_2   0.848   1.000
## 
## Phylogenetic signal (OU process):
##             d
## trait_1 0.176
## trait_2 0.857
## 
## Coefficients:
##                     Estimate      SE Z-score P-value    
## trait_1_0            -0.1234  0.2063   -0.60    0.55    
## trait_2_0            -0.3551  0.9341   -0.38    0.70    
## trait_2_cov_trait_2   1.0102  0.0159   63.67  <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The output of `cor_phylo()` includes log likelihood values, AIC, BIC, estimated correlation matrix of traits, estimated phylogenetic signals of traits, estimated coefficients and their uncertainties (standard errors, Z scores, and p values). If bootstrapping was enabled by setting the `boot` argument, the lower and upper boundaries of correlations, phylogenetic signal values, and coefficients will be appended. 

# Closing remarks

In recent years, there has been increasing effort to apply model-based approaches in community ecology. One particular example is the joint species distribution models. Despite the well-known importance of phylogenetic relationship in structuring species distributions and community composition, relatively few studies have incorporated phylogenetic relationships in model-based analyses of species distributions and community ecology. One potential reason is the lack of easy-to-use tools to facilitate the usage of phylogenetic species-distribution modelling in ecological communities. The package `phyr` fills this gap by providing implementations of phylogenetic species-distribution models with flexible model formula syntax (`pglmm()`). It also includes other model-based functions that are useful for ecological studies such as calculating community phylogenetic diversity (e.g. `psv()`) and estimating correlations among functional traits while accounting for their evolutionary history (`cor_phylo()`) (Table \@ref(tab:mainFunc)). 

The model formula of `pglmm()` is general and can be applied to other tools to fit phylogenetic species-distribution models. In the future, we plan to develop interfaces to fit such models with more bayesian programs such as Stan (based on the R package `brm`) [@brmsR2018]. It is our hope that the `phyr` package and the proposed model formula for phylogenetic species-distribution models will facilitate the usage of such models by other researchers.


<!-- # Example analysis

To demonstrate the usage of main functions in `phyr`, we simulated a dataset with 50 species and 30 communities. We first simulated a coalescent phylogeny of 50 species with function `ape::rcoal`. For each species, we then simulated one continuous functional trait along the phylogeny. We also simulated one environmental variable with all 30 communities located evenly along the gradient. Environmental variable, functional trait, and their interaction all determine the abundance of species over sites following the model below. 

$$y=(\beta_{0}+b_{0})+(\beta_{1}+b_{1})envi+(\beta_{2}+b_{2})trait+(\beta_{3}+b_{3})envi*trait+e$$
$$b_{0}\sim Gaussian(0,\sigma_{0}^{2}\Sigma_{spp})$$
$$b_{1}\sim Gaussian(0,\sigma_{1}^{2}\Sigma_{spp})$$
$$b_{2}\sim Gaussian(0,\sigma_{2}^{2}\Sigma_{site})$$
$$e\sim Gaussian(0,\sigma_{e}^2)$$

Here, we set all coefficients (beta_0 to beta_3) to be 1; we also set variances of all random terms to be 1. $\Sigma_{spp}$ is a variance-covariance matrix converted from the phylogeny. $\Sigma_{site}$ is an identity matrix since we treat sites as independent replications. In the above model, different species have different overall abundance (intercept) but closely related species have similar overall abundance. Abundance of different species also change differently along the environmental gradient; however, closely related species again change similarly in response to environmental change. We also included a random term for functional trait; therefore the relationship between functional traits and species abundance varied across sites independently to account for unmeasured environmental variables. For demonstration purpose, we only simulated one functional trait and one environmental variable; in real datasets, multiple environmental variables and multiple functional traits can be included in the model. However, more random terms will come with higher computational burdens and the program will take longer to finish.  -->
  


# Acknowledgements

Funding for this work was provided by the National Science Foundation (NSF ...).

# Authors' contributions

D.L and A.R.I conceived the idea. All authors wrote the software and package documentations. All authors wrote the manuscript.

# Data Accessibility

No data were used in this study. `phyr` is available at Github (https://github.com/daijiang/phyr) and will be submitted to CRAN in the near future.

# References
