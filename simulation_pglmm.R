# devtools::install_github("daijiang/phyr")

library(ape)

#######################################################
############## Set up simulation dataset ##############
#######################################################

###### Set up species and site variables
nspp <- 5
nsite <- 10

# construct phylogeny from 'ape'
phy <- rcoal(nspp)
phy <- compute.brlen(phy, method = "Grafen", power = 0.5)
# set up phylogenetic attraction covariance matrix which requires a nested covariance matrix for phylogenetic attraction
Vphy <- vcv(phy)
Vphy <- Vphy / max(Vphy)
Vphy <- Vphy / exp(determinant(Vphy)$modulus[1] / nspp)

# Perform a Cholesky decomposition of Vphy. This is used to generate phylogenetic
# signal: a vector of independent normal random variables, when multiplied by the
# transpose of the Cholesky deposition of Vphy will have covariance matrix
# equal to Vphy. Another way is to use mvtnorm package.
cholVnested = t(chol(kronecker(diag(1, nrow = nsite), Vphy))) # nested component


###### Set up covaraince matrices for simulations of spatial and nested covariances
# create spatial coordinates for sites
# (to make a visually more clear spatial covariance matrix, sites are placed on a line)
x.coord <- rep(0, nsite)
y.coord <- sort(runif(n = nsite))

# Set up spatial correlation covariance matrix
Dist <- matrix(0, nrow = nsite, ncol = nsite)
for (i in 1:nsite)
  for (j in 1:nsite)
    Dist[i, j] <-
  ((x.coord[i] - x.coord[j]) ^ 2 + (y.coord[i] - y.coord[j]) ^ 2) ^ .5

range <- .25
sd.space <- 1
Vspace <- sd.space ^ 2 * exp(-Dist / range)
rownames(Vspace) <- 1:nsite
colnames(Vspace) <- 1:nsite

cholVspace = t(chol(kronecker(Vspace, diag(1, nrow = nspp)))) # nested spatial correlation component

# generate environmental site variable
envi <- matrix(1:nsite, nrow = 1, ncol = nsite)
envi <- (envi - mean(envi)) / sd(envi)

trait <- rTraitCont(phy, sigma = .1)

# create data.frame with species, sites, trait, and environmental variable
dat = data.frame(
  species = rep(phy$tip.label, nsite),
  site = rep(1:nsite, each = nspp),
  x.coord = rep(x.coord, each = nspp),
  y.coord = rep(y.coord, each = nspp),
  envi = rep(envi, each = nspp),
  trait = rep(trait, times = nsite)
)

###### Set up model parameters for simulations

# parameters
beta0 <- beta1 <- beta2 <- beta3 <- 1 # fixed terms
sd.b0 <- sd.b1 <- sd.b2 <- sd.nested <- sd.resid <- 1 # sd of random terms
sd.nested <- 1

# whether or not to include phylogenetic signal in B0 and B1
signal.b0 <- TRUE
signal.b1 <- TRUE

# set up species-specific regression coefficients as random effects

if (signal.b0) {
  b0 <- rTraitCont(phy, sigma = sd.b0 / 10)
} else {
  b0 <- rnorm(nspp, sd = sd.b0)
}
if (signal.b1) {
  b1 <- beta1 + rTraitCont(phy, sigma = sd.b1 / 10)
} else {
  b1 <- beta1 + rnorm(nspp, sd = sd.b1)
}
b2 <- beta2 + rnorm(nsite, sd = sd.b2)

b0_site = t(chol(Vspace)) %*% (rnorm(nsite, mean = 0, sd = sd.b0))


# set up data.frame for parameters
B <- data.frame(
  b0_spp = rep(b0, nsite),
  b0_site = rep(b0_site[,1], each = nspp),
  b1 = rep(b1, nsite),
  b2 = rep(b2, each = nspp),
  beta3 = beta3
)

###### Simulate data

y <- beta0 + B$b0_spp + B$b0_site +
  B$b1 * dat$env +
  B$b2 * dat$trait +
  beta3 * dat$env * dat$trait +
  cholVnested %*% rnorm(nsite * nspp, sd = sd.nested) #+
  # cholVspace %*% rnorm(nsite * nspp, sd = sd.space)

# Create both presence-absence and continuous abundance data
dat$presabs <- rbinom(n = length(y),
                      size = 1,
                      prob = exp(y) / (1 + exp(y)))
dat$abund <- y + rnorm(nspp * nsite, sd = sd.resid)
head(dat)

#######################################################
############## Analyze simulation dataset #############
#######################################################

library(phyr)

# Required components:
#	dat = data.frame containing species, sites, environmental variables, and trait values
#	phy = a phylo{ape} object for the species phylogeny, or alternatively a covariance matrix for species
#	Vspace = a covariance matrix for spatial covariances among sites

# An example with phylogenetic random effects for species-specific slopes for the effect of the environment
z <- pglmm(
  abund ~ 1 + envi + trait + envi:trait +
    (1 | species__)  + (1 | site__) +
    (envi | species__) + (trait | site) +
    (1 | species__@site),
  data = dat,
  cov_ranef = list(species = phy, site = Vspace)
)
# 
# z.bayes <- pglmm(
#   abund ~ 1 + envi + trait + envi:trait +
#     (1 | species__)  + (1 | site__) + 
#     (envi | species__) + (trait | site) +
#     (1 | species__@site),
#   data = dat, bayes = T,
#   cov_ranef = list(species = phy, site = Vspace),
#   s2.init = sqrt(z$ss)
# )
# 
pdf("designPlot.pdf", width = 10, height = 10)
pglmm.plot.re(
  x = z,
  sp.var = "species",
  site.var = "site",
  show.image = TRUE,
  show.sim.image = FALSE
)
dev.off()

# # finding the optimal `range`
# z.gaussian <- pglmm(abund ~ 1 + (1 | site__),
#                     data = dat,
#                     cov_ranef = list(site = Vspace))
# summary(z.gaussian)
# # if the ranef `site__` is large, it makes to investigate the value of range
# for (range in .1 * (1:9)) {
#   Vspace <- exp(-Dist / range)
#   rownames(Vspace) <- 1:nsite
#   colnames(Vspace) <- 1:nsite
#   
#   z.gaussian <- pglmm(abund ~ 1 + (1 | site__),
#                       data = dat,
#                       cov_ranef = list(site = Vspace))
#   show(c(range, z.gaussian$logLik))
# }
