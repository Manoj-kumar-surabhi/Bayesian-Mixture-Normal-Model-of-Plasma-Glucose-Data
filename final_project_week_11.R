##1
library(coda)
library(MASS)  

alpha <- 1
beta <- 1
mu0 <- 120
tau0_sq <- 200
sigma0_sq <- 1000
v0 <- 10
S <- 5000 

y <- rnorm(100, mean = 100, sd = 15)  
n <- length(y)
median_y <- median(y)

theta1 <- numeric(S)
theta2 <- numeric(S)
sigma1_sq <- numeric(S)
sigma2_sq <- numeric(S)
pi <- numeric(S)
X <- matrix(0, n, S)

X[,1] <- ifelse(y < median_y, 1, 2)
sigma1_sq[1] <- var(y[X[,1] == 1])
sigma2_sq[1] <- var(y[X[,1] == 2])
theta1[1] <- mean(y[X[,1] == 1])
theta2[1] <- mean(y[X[,1] == 2])
pi[1] <- mean(X[,1] == 1)

# Gibbs sampler
for (s in 2:S) {
  n1 <- sum(X[,s-1] == 1)
  n2 <- n - n1
  
  
  if (n1 > 0) {
    mean1 <- (n1 * mean(y[X[,s-1] == 1]) + mu0 / tau0_sq) / (n1 + 1 / tau0_sq)
  } else {
    mean1 <- mu0
  }
  
  if (n2 > 0) {
    mean2 <- (n2 * mean(y[X[,s-1] == 2]) + mu0 / tau0_sq) / (n2 + 1 / tau0_sq)
  } else {
    mean2 <- mu0
  }
  

  var1 <- sigma1_sq[s-1] / (n1 + 1 / tau0_sq) + 1e-10
  var2 <- sigma2_sq[s-1] / (n2 + 1 / tau0_sq) + 1e-10
  

  theta1[s] <- rnorm(1, mean = mean1, sd = sqrt(var1))
  theta2[s] <- rnorm(1, mean = mean2, sd = sqrt(var2))
  

  ss1 <- sum((y[X[,s-1] == 1] - theta1[s])^2)
  ss2 <- sum((y[X[,s-1] == 2] - theta2[s])^2)
  
  sigma1_sq[s] <- 1 / rgamma(1, alpha + n1 / 2, beta + ss1 / 2)
  sigma2_sq[s] <- 1 / rgamma(1, alpha + n2 / 2, beta + ss2 / 2)
  

  pi[s] <- rbeta(1, alpha + n1, beta + n2)
  
  for (i in 1:n) {
    p1 <- dnorm(y[i], mean = theta1[s], sd = sqrt(sigma1_sq[s])) * pi[s]
    p2 <- dnorm(y[i], mean = theta2[s], sd = sqrt(sigma2_sq[s])) * (1 - pi[s])
    
    if (is.na(p1) || is.na(p2)) {
      X[i, s] <- sample(1:2, 1)  # Random choice if probabilities are invalid
    } else {
      X[i, s] <- sample(1:2, 1, prob = c(p1, p2))
    }
  }
}

theta1_post <- pmin(theta1, theta2)  # theta(1)
theta2_post <- pmax(theta1, theta2)  # theta(2)

split_point <- S / 2
theta1_first <- theta1_post[1:split_point]
theta1_second <- theta1_post[(split_point + 1):S]

# Plot histogram and kernel density
par(mfrow = c(1, 1))
hist(theta1_first, probability = TRUE, main = "Theta1 - First Half")
lines(density(theta1_first), col = "blue")
hist(theta1_second, probability = TRUE, main = "Theta1 - Second Half")
lines(density(theta1_second), col = "blue")

# Posterior mean and 95% credible interval
theta1_mean_first <- mean(theta1_first)
theta1_ci_first <- quantile(theta1_first, probs = c(0.025, 0.975))
theta1_mean_second <- mean(theta1_second)
theta1_ci_second <- quantile(theta1_second, probs = c(0.025, 0.975))

cat("Posterior Mean and 95% CI (First Half):", theta1_mean_first, theta1_ci_first, "\n")
cat("Posterior Mean and 95% CI (Second Half):", theta1_mean_second, theta1_ci_second, "\n")

# Trace plot
plot(theta1_post, type = 'l', main = "Trace plot for Theta(1)", ylab = "Theta(1)")
plot(theta2, type = 'l', main = "Trace plot for Theta(2)", ylab = "Theta(2)")
# Autocorrelation and effective sample size
acf(theta1_post, main = "Autocorrelation for Theta(1)")
acf(theta2_post, main = "Autocorrelation for Theta(2)")
ess <- effectiveSize(theta1_post)
cat("Effective Sample Size:", ess, "\n")

##2

alpha <- 1
beta <- 1
mu0 <- 120
tau0_sq <- 200
sigma0_sq <- 1000
v0 <- 10
S <- 5000  
n_chains <- 3  
n <- 100  
y <- rnorm(n, mean = 100, sd = 15)  


initialize_chain <- function(y) {
  X <- ifelse(runif(length(y)) < 0.5, 1, 2)
  sigma1_sq <- var(y[X == 1])
  sigma2_sq <- var(y[X == 2])
  theta1 <- mean(y[X == 1])
  theta2 <- mean(y[X == 2])
  pi <- mean(X == 1)
  list(X = X, sigma1_sq = sigma1_sq, sigma2_sq = sigma2_sq, theta1 = theta1, theta2 = theta2, pi = pi)
}

gibbs_sampler <- function(y, S, init) {
  n <- length(y)
  theta1 <- numeric(S)
  theta2 <- numeric(S)
  sigma1_sq <- numeric(S)
  sigma2_sq <- numeric(S)
  pi <- numeric(S)
  X <- matrix(0, n, S)
  
  X[,1] <- init$X
  sigma1_sq[1] <- init$sigma1_sq
  sigma2_sq[1] <- init$sigma2_sq
  theta1[1] <- init$theta1
  theta2[1] <- init$theta2
  pi[1] <- init$pi
  
  for (s in 2:S) {
    # Count number of observations in each group
    n1 <- sum(X[,s-1] == 1)
    n2 <- n - n1
    
    # Sample theta1 and theta2
    if (n1 > 0) {
      mean1 <- (n1 * mean(y[X[,s-1] == 1]) + mu0 / tau0_sq) / (n1 + 1 / tau0_sq)
      var1 <- sigma1_sq[s-1] / (n1 + 1 / tau0_sq)
    } else {
      mean1 <- mu0
      var1 <- sigma1_sq[s-1] / (1 / tau0_sq)
    }
    theta1[s] <- rnorm(1, mean = mean1, sd = sqrt(var1))
    
    if (n2 > 0) {
      mean2 <- (n2 * mean(y[X[,s-1] == 2]) + mu0 / tau0_sq) / (n2 + 1 / tau0_sq)
      var2 <- sigma2_sq[s-1] / (n2 + 1 / tau0_sq)
    } else {
      mean2 <- mu0
      var2 <- sigma2_sq[s-1] / (1 / tau0_sq)
    }
    theta2[s] <- rnorm(1, mean = mean2, sd = sqrt(var2))
    

    ss1 <- sum((y[X[,s-1] == 1] - theta1[s])^2)
    ss2 <- sum((y[X[,s-1] == 2] - theta2[s])^2)
    sigma1_sq[s] <- 1 / rgamma(1, alpha + n1 / 2, beta + ss1 / 2)
    sigma2_sq[s] <- 1 / rgamma(1, alpha + n2 / 2, beta + ss2 / 2)
    

    pi[s] <- rbeta(1, alpha + n1, beta + n2)
    

    for (i in 1:n) {
      p1 <- dnorm(y[i], mean = theta1[s], sd = sqrt(sigma1_sq[s])) * pi[s]
      p2 <- dnorm(y[i], mean = theta2[s], sd = sqrt(sigma2_sq[s])) * (1 - pi[s])
      if (is.na(p1) || is.na(p2)) {
        X[i, s] <- sample(1:2, 1)
      } else {
        X[i, s] <- sample(1:2, 1, prob = c(p1, p2))
      }
    }
  }
  
  list(theta1 = theta1, theta2 = theta2, sigma1_sq = sigma1_sq, sigma2_sq = sigma2_sq, pi = pi, X = X)
}

chains <- vector("list", n_chains)
for (chain in 1:n_chains) {
  init <- initialize_chain(y)
  chains[[chain]] <- gibbs_sampler(y, S, init)
}

# Diagnostics for Theta(2) in each chain
par(mfrow = c(1, 1))
for (chain in 1:n_chains) {
  chain_theta2 <- chains[[chain]]$theta2
  
  # (a) Histogram and Kernel Density
  hist(chain_theta2, probability = TRUE, main = paste("Theta(2) Histogram - Chain", chain))
  lines(density(chain_theta2), col = "blue")
  
  # (b) Posterior Mean and 95% CI
  theta2_mean <- mean(chain_theta2)
  theta2_ci <- quantile(chain_theta2, probs = c(0.025, 0.975))
  cat("Chain", chain, "- Theta(2) Mean:", theta2_mean, "\n")
  cat("Chain", chain, "- Theta(2) 95% CI:", theta2_ci, "\n")
  
  # (c) Trace Plot
  plot(chain_theta2, type = 'l', main = paste("Trace plot for Theta(2) - Chain", chain), ylab = "Theta(2)")
  
  # (d) Autocorrelation and ESS
  acf(chain_theta2, main = paste("Autocorrelation for Theta(2) - Chain", chain))
  ess <- effectiveSize(chain_theta2)
  cat("Chain", chain, "- Effective Sample Size:", ess, "\n")
}


