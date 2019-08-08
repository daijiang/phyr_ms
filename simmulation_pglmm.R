# devtools::install_github("daijiang/phyr")

library(ape)
library(LatticeKrig)

#######################################################
############## Set up simulation dataset ##############
#######################################################

###### Set up species and site variables
nspp <- 10
nsite <- 20

# construct phylogeny from 'ape'
phy <- rcoal(nspp)
phy <- compute.brlen(phy, method = "Grafen", power = 0.5)

trait <- rTraitCont(phy, sigma = .1)

# generate environmental site variable
envi <- matrix(1:nsite, nrow = 1, ncol = nsite)
envi <- (envi - mean(envi))/sd(envi)

# create spatial coordinates for sites
x.coord <- runif(n=nsite)
y.coord <- runif(n=nsite)

dat = data.frame(species = rep(phy$tip.label, nsite),
                 site = rep(1:nsite, each = nspp),
                 x.coord = rep(x.coord, each = nspp),
                 y.coord = rep(y.coord, each = nspp),
                 envi = rep(envi, each = nspp),
                 trait = rep(trait, times=nsite))

###### Set up model parameters for simulations

# parameters
beta0 <- beta1 <- beta2 <- beta3 <- 1 # fixed terms
sd.b0 <- sd.b1 <- sd.b2 <- sd.nested <- sd.resid <- 1 # sd of random terms
# set sd.resid to zero for binary data
binary_data = FALSE
if(binary_data) sd.resid <- 0

# whether or not to include phylogenetic signal in B0 and B1
signal.b0 <- TRUE
signal.b1 <- TRUE

# set up species-specific regression coefficients as random effects
if (signal.b0) {
  b0 <- beta0 + rTraitCont(phy, sigma = sd.b0/10)
} else {
  b0 <- beta0 + rnorm(nspp, sd = sd.b0)
}
if (signal.b1) {
  b1 <- beta1 + rTraitCont(phy, sigma = sd.b1/10)
} else {
  b1 <- beta1 + rnorm(nspp, sd = sd.b1)
}
b2 <- beta2 + rnorm(nsite, sd = sd.b2)

# set up data.frame for parameters
B <- data.frame(b0=rep(b0, nsite),
                b1=rep(b1, nsite),
                b2=rep(b2, each = nspp),
                beta3=beta3)

###### Set up covaraince matrices for simulations of spatial and nested covariances

# Set up spatial correlation covariance matrix
sd.space <- 1
Dist <- matrix(0, nrow=nsite, ncol=nsite)
for(i in 1:nsite) for(j in 1:nsite) Dist[i,j] <- ((x.coord[i] - x.coord[j])^2 + (y.coord[i] - y.coord[j])^2)^.5

range <- 3
Vspace <- exp(-Dist/range)
rownames(Vspace) <- 1:nsite
colnames(Vspace) <- 1:nsite

cholVspace = t(chol(kronecker(Vspace, diag(1, nrow = nspp)))) # nested spatial correlation component

# set up phylogenetic attraction covariance matrix which requires a nested covariance matrix for phylogenetic attraction
sd.nested <- 1
Vphy <- vcv(phy)
Vphy <- Vphy/max(Vphy)
Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/nspp)

# Perform a Cholesky decomposition of Vphy. This is used to generate phylogenetic
# signal: a vector of independent normal random variables, when multiplied by the
# transpose of the Cholesky deposition of Vphy will have covariance matrix
# equal to Vphy. Another way is to use package mvnorm.
cholVnested = t(chol(kronecker(diag(1, nrow = nsite), Vphy))) # nested component

###### Simulate data

y <- B$b0 + 
  B$b1 * dat$env + 
  B$b2 * dat$trait + 
  beta3 * dat$env * dat$trait + 
  cholVnested %*% rnorm(nsite * nspp, sd = sd.nested) + 
  cholVspace %*% rnorm(nsite * nspp, sd = sd.space)

# Create both presence-absence and continuous abundance data
dat$presabs <- rbinom(n=length(y), size=1, prob=exp(y)/(1+exp(y)))
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

# without nested term
z.binary <- pglmm(presabs ~ 1 + envi + trait + envi:trait + 
                    (1|species__) + (1|site__), data = dat, 
                  cov_ranef = list(species = phy, site = Vspace), 
                  family = "binomial", add.obs.re = F)

z.gaussian <- pglmm(abund ~ 1 + envi + trait + envi:trait + 
                      (1|species__) +  (1|site__), data = dat, 
                    cov_ranef = list(species = phy, site = Vspace))

# with nested term
z.binary.nested <- pglmm(presabs ~ 1 + envi + trait + envi:trait + 
                           (1|species__) + (1|site__) + (1|species__@site), 
                         data = dat, cov_ranef = list(species = phy, site = Vspace), 
                         family = "binomial", add.obs.re = F)

z.gaussian.nested <- pglmm(abund ~ 1 + envi + trait + envi:trait + 
                             (1|species__) +  (1|site__) + (1|species__@site), 
                           data = dat, cov_ranef = list(species = phy, site = Vspace))

summary(z.binary)   
summary(z.gaussian)   
summary(z.binary.nested)   
summary(z.gaussian.nested)   

phyr::communityPGLMM.plot.re(x=z.binary.nested, sp.var = "species", site.var = "site",
                             show.image = TRUE, show.sim.image = FALSE)

# An example with phylogenetic random effects for species-specific slopes for the effect of the environment     
z.gaussian.re.env <- pglmm(abund ~ 1 + envi + trait + envi:trait + 
                             (1|species__) + (0 + envi|species__) + (0 + trait|site) +
                             (1|species__@site), data = dat, cov_ranef = list(species = phy))

p1 = communityPGLMM.plot.re(x = z.gaussian.re.env, sp.var = "species", site.var = "site", 
                            show.image = TRUE, show.sim.image = FALSE, colorkey = FALSE)
p1$vcv
