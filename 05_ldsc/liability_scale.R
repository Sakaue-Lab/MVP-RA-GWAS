# convert into liability scale

# Inputs
h2_obs <- 0.0066
se_obs <- 0.0019
K <- 0.005 # population prevalence
P <- (5240) / 449398 # cases in sample

# Threshold on liability scale
t <- qnorm(1 - K)

# Standard normal density at threshold
z <- dnorm(t)

# Liability-scale heritability
h2_liab <- h2_obs * (K^2 * (1 - K)^2) / (z^2 * P * (1 - P))
se_liab<- se_obs * (K^2 * (1 - K)^2) / (z^2 * P * (1 - P))
