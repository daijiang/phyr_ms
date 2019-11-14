library(ape)

set.seed(1000)
# Set up parameter values for simulating data
n <- 50
phy <- rcoal(n, tip.label = 1:n)
trt_names <- paste0("par", 1:2)

R <- matrix(c(1, 0.7, 0.7, 1), nrow = 2, ncol = 2) # true var-corv matrix
d <- c(0.3, 0.95) # true phylo signal
B2 <- 1

Se <- c(0.2, 1) # measurement error
M <- matrix(Se, nrow = n, ncol = 2, byrow = TRUE)
colnames(M) <- trt_names

# Set up needed matrices for the simulations
p <- length(d)

star <- stree(n)
star$edge.length <- array(1, dim = c(n, 1))
star$tip.label <- phy$tip.label

Vphy <- vcv(phy)
Vphy <- Vphy/max(Vphy)
Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)

tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy
C <- matrix(0, nrow = p * n, ncol = p * n)
for (i in 1:p) for (j in 1:p) {
  Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
  C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
}
MM <- matrix(M^2, ncol = 1)
V <- C + diag(as.numeric(MM))
iD <- t(chol(V))

XX <- iD %*% rnorm(2 * n)
X <- matrix(XX, n, p)
colnames(X) <- trt_names

U <- list(cbind(rnorm(n, mean = 2, sd = 10)))
names(U) <- trt_names[2]

X[,2] <- X[,2] + B2[1] * U[[1]][,1] - B2[1] * mean(U[[1]][,1])

traits = as.data.frame(X)
names(traits) = c("trait_1", "trait_2")
traits$sp = phy$tip.label
traits$cov_trait_2 = U$par2[,1]
traits$se_trait_1 = M[,1]
traits$se_trait_2 = M[,2]
