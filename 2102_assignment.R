### 2102 Assignment ###
install.packages("moments")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyr")
library(moments)
library(ggplot2)
library(dplyr)
library(tidyr)

# Given information
y_rate <- 0.18
o_rate <- 0.12

# Read in policyholder data
data <- read.csv(file = "AgeLevelData.csv", header = T)


# Discount schemes
discount_0 <- c(1.3, 1.15, 1, 0.9, 0.8)
discount_1 <- c(1.15, 1, 0.9)
discount_2 <- c(1.45, 1.3, 1.15, 1, 0.9, 0.8, 0.7)


## Task 1 ##

# Question 1 #

# Existing scheme:


# A function that returns the probability transition matrix for the option 0
# scheme based on the policy holders rate parameter

zero_P <- function(rate) {
  no_claim <- dpois(0, rate)
  one_claim <- dpois(1, rate)
  two_claim <- 1 - no_claim - one_claim
  p_claim <- c(no_claim, one_claim, two_claim)
  
  P <- matrix(c(p_claim[2] + p_claim[3], p_claim[1], 0, 0, 0,
                p_claim[2] + p_claim[3], 0, p_claim[1], 0, 0,
                p_claim[3], p_claim[2], 0, p_claim[1], 0,
                0, p_claim[3], p_claim[2], 0, p_claim[1],
                0, 0, p_claim[3], p_claim[2], p_claim[1]), 
                byrow = T, nrow = 5)
  return(P)
}

zero_P_old <- zero_P(o_rate)
zero_P_young <- zero_P(y_rate)

# P matrices under option 0
show(zero_P_old)
show(zero_P_young)


# Option 1:

# A function that returns the probability transition matrix for the option 1 
#scheme based on the policy holders rate parameter

one_P <- function(rate) {
  no_claim <- dpois(0, rate)
  one_claim <- dpois(1, rate)
  two_claim <- 1 - no_claim - one_claim
  
  p_claim <- c(no_claim, one_claim, two_claim)
  
  P <- matrix(c(p_claim[2] + p_claim[3], p_claim[1], 0,
                p_claim[2] + p_claim[3], 0, p_claim[1],
                p_claim[3], p_claim[2], p_claim[1]),
                byrow = T, nrow = 3)
  return(P)
}

one_P_old <- one_P(o_rate)
one_P_young <- one_P(y_rate)

# P matrices under option 1
show(one_P_old)
show(one_P_young)

# Option 2:

# A function that returns the probability transition matrix for the option 2 
#scheme based on the policy holders rate parameter
two_P <- function(rate) {
  no_claim <- dpois(0, rate)
  one_claim <- dpois(1, rate)
  two_claim <- 1 - no_claim - one_claim
  
  p_claim <- c(no_claim, one_claim, two_claim)
  
  P <- matrix(c(p_claim[2] + p_claim[3], p_claim[1], 0, 0, 0, 0, 0,
                p_claim[2] + p_claim[3], 0, p_claim[1], 0, 0, 0, 0,
                p_claim[3], p_claim[2], 0, p_claim[1], 0, 0, 0, 
                0, p_claim[3], p_claim[2], 0, p_claim[1], 0, 0,
                0, 0, p_claim[3], p_claim[2], 0, p_claim[1], 0,
                0, 0, 0, p_claim[3], p_claim[2], 0, p_claim[1],
                0, 0, 0, 0, p_claim[3], p_claim[2], p_claim[1]),
                byrow = T, nrow = 7)
  
  return(P)
}

two_P_young <- two_P(y_rate)

# P matrices under option 2
show(zero_P_old)
show(two_P_young)

# Question 2 #

# For simulating 1500 times
simlength <- 1500

# Initialise a count for young policy holders in each group
young_counts <- matrix(0, nrow = simlength, ncol = 5)
colnames(young_counts) <- c("-2", "-1", "0", "1", "2")

# Initialise a count for old policy holders in each group
old_counts <- matrix(0, nrow = simlength, ncol = 5)
colnames(old_counts) <- c("-2", "-1", "0", "1", "2")

# Initialise the total claims vector
total_claims <- numeric(simlength)

for (i in 1 : simlength) {
  data$claims <- ifelse(
    data$age >= 20 & data$age <= 25,
    rpois(n = nrow(data), lambda = 0.18),  # Young group with rate 0.18
    rpois(n = nrow(data), lambda = 0.12)   # Old group with rate 0.12
  )
  
  # Add a column to data for the new NCD level
  data$ncdlevel25 <- pmin(2, pmax(-2, data$ncdlevel24 + ifelse(data$claims == 0, 1, ifelse(data$claims == 1, -1, -2))))
  
  # Note the total number of claims in each sample path
  total_claims[i] <- sum(data$claims)
  
  # Separate young and old dataset
  young_data <- data[data$age >= 20 & data$age <= 25, ]
  old_data <- data[data$age >= 50 & data$age <= 55, ]
  
  # Count policyholders at each level for both groups and store in respective matrices
  young_counts[i, ] <- as.numeric(table(factor(young_data$ncdlevel25, levels = -2:2)))
  old_counts[i, ] <- as.numeric(table(factor(old_data$ncdlevel25, levels = -2:2)))

}

# Find the average distribution of policy holders in the levels
level_dist <- colMeans(young_counts + old_counts)
show(level_dist)


# Calculate the premiums paid by each group for each discount level
young_premium <- (500 * young_counts) %*% discount_0
old_premium <- (450 * old_counts) %*% discount_0
total_premium <- young_premium + old_premium

# Summary statistics for the total premium
total_premium_summary <- c(
  Mean = mean(total_premium),
  SD = sd(total_premium),
  Skew = skewness(total_premium),
  Kurtosis = kurtosis(total_premium),
  Min = min(total_premium),
  Max = max(total_premium)
)

show(total_premium_summary)

expected_profit <- mean(total_premium) - mean(total_claims) * 2200
show(expected_profit)

## Task 2 ##

# Question 1 #

# Functiions used as a part of calculating Loimaranta Efficiency
stationary_p <- function(scheme) {
  I <- diag(nrow(scheme))
  A <- t(scheme) - I
  A <- rbind(A, rep(1, nrow(scheme)))
  A <- A[-1, ]
  b <- c(rep(0, nrow(scheme) - 1), 1)
  
  stationary_p <- solve(A, b)
  return(stationary_p)
}
p_lambda <- function(schemeP, schemeD) {
  I <- diag(nrow(schemeP))
  A <- t(schemeP) - I
  A <- rbind(A, rep(1, nrow(schemeP)))
  A <- A[-1, ]
  b <- c(rep(0, nrow(schemeP) - 1), 1)
  
  stationary_p <- solve(A, b)
  
  return(as.numeric(stationary_p %*% schemeD))
}
dp_lambda <- function(p_lambda, lambda, delta, schemeP, schemeD) {
  return(as.numeric((p_lambda(schemeP(lambda + delta), schemeD) - p_lambda(schemeP(lambda), schemeD)) / delta))
}
loimaranta <- function(p_lambda, schemeP, schemeD, delta) {
  # Sequence of lambda values from 0.05 to 0.2, incremented by 0.005
  lambda_values <- seq(0, 1.5, by = 0.01)
  
  # Initialize a vector to store results
  loimaranta_values <- numeric(length(lambda_values))
  
  # Loop over each lambda value
  for (i in seq_along(lambda_values)) {
    lambda <- lambda_values[i]
    
    # Compute dp_lambda for the current lambda
    dp_lambda_val <- dp_lambda(p_lambda, lambda, delta, schemeP, schemeD)
    
    # Compute p_lambda for the current lambda
    p_lambda_val <- p_lambda(schemeP(lambda), schemeD)
    
    # Calculate lambda * dp_lambda / p_lambda
    loimaranta_values[i] <- (lambda / p_lambda_val) * dp_lambda_val
  }
  # Return a data frame with lambda values and their corresponding Loimaranta values
  result <- data.frame(lambda = lambda_values, loimaranta = loimaranta_values)
  return(result)
}

# Loim eff. for each scheme
loim0 <- loimaranta(p_lambda, zero_P, discount_0, 0.01)
loim1 <- loimaranta(p_lambda, one_P, discount_1, 0.01)
loim2 <- loimaranta(p_lambda, two_P, discount_2, 0.01)

loims <- data.frame(
  lambda = loim0$lambda,
  scheme1 = loim0$loimaranta,
  scheme2 = loim1$loimaranta,
  scheme3 = loim2$loimaranta
)

#Loim eff. plot
ggplot(loims, aes(x = lambda)) +
  geom_line(aes(y = scheme1, color = "Existing Scheme"), size = 1) +
  geom_line(aes(y = scheme2, color = "Simplified Scheme"), size = 1) +
  geom_line(aes(y = scheme3, color = "Complex Scheme"), size = 1) +
  labs(
    title = "Loimaranta Efficiency by Rate Parameter for Each Scheme",
    x = expression(paste("Rate Parameter ", lambda)),
    y = expression(paste("Loimaranta Efficiency ", eta(lambda)))
  ) +
  # Adding vertical lines with distinct color and linetype for a separate legend
  geom_vline(aes(xintercept = 0.18, linetype = "Young"), color = "black", size = 0.5) +
  geom_vline(aes(xintercept = 0.12, linetype = "Old"), color = "grey", size = 0.5) +
  
  # Customizing color for scheme lines only
  scale_color_manual(
    name = "Scheme",
    values = c(
      "Existing Scheme" = 'red',
      "Simplified Scheme" = 'blue',
      "Complex Scheme" = 'green'
    )
  ) +
  # Customizing linetype for Lambda lines only
  scale_linetype_manual(
    name = "Policy Holder Group",
    values = c("Young" = "dashed", "Old" = "dashed")
  ) +
  
  theme_minimal()


ggsave("actl2102_assignment_2.1.png", height = 10, width = 15, units = 'cm', dpi = 1080)

# Question 2 #

# Stationary probs
pi_0_old <- stationary_p(zero_P_old)
pi_0_young <- stationary_p(zero_P_young)
pi_1_old <- stationary_p(one_P_old)
pi_1_young <- stationary_p(one_P_young)
pi_2_young <- stationary_p(two_P_young)

# Expected long run ann. profit 
profit0y <- (5000 * pi_0_young) %*% (500 * discount_0) - 2200 * 900
profit0o <- (5000 * pi_0_old) %*% (450 * discount_0) - 2200 * 600
profit1y <- (5000 * pi_1_young) %*% (500 * discount_1) - 2200 * 900
profit1o <- (5000 * pi_1_old) %*% (450 * discount_1) - 2200 * 600
profit2y <- (5000 * pi_2_young) %*% (500 * discount_2) - 2200 * 900

long_run_profit0 <- profit0y + profit0o
long_run_profit1 <- profit1y + profit1o
long_run_profit2 <- profit2y + profit0o

show(long_run_profit0)
show(long_run_profit1)
show(long_run_profit2)
