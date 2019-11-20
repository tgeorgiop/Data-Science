library(mvnfast); library(mnormt); library(invgamma)
CalculatePosterior <- function(y, x, prior_parameters, g){
  # Args:
  #   y: a vector with length n containing the observed values of the response variable y
  #   x: an n × q matrix x containing the observed design matrix X
  #   prior_parameters: 
  #     a: prior parameter (scalar), a>0
  #     b: prior parameter (scalar), b>0
  #     mu: prior parameter (vector), mu
  #   g: positive scalar g > 0 that corresponds to the Zellner’s g hyper-parameter 
  #   
  #Returns:
  #   posterior_parameters:
  #     a: parameter (scalar), a_tilde > 0 of the posterior distribution
  #     b: parameter (scalar), b_tilde > 0 of the posterior distribution
  #     beta: parameter (vector), beta_tilde > 0 of the posterior distribution
  #     Sigma: parameter (q×q pos def sym matrix) posterior Sigma_tilde
  #   posterior_mean: posterior mean of betas and sigma^2
  #   posterior_var: posterior variance of betas and sigma^2
  #   credible_intervals: 95% and 99% credible_intervals of betas and sigma^2
  #   marginal_y: tha natural logarithm of the marginal likelihood m(y)
  V <- g * solve(t(x) %*% x)
  
  # Calculate posterior parameters
  posterior_Sigma <- solve(t(x) %*% x + solve(V))
  posterior_beta <- posterior_Sigma %*% 
                                 (t(x) %*% y + solve(V) %*% prior_parameters$mu)
  posterior_a <- prior_parameters$a + n/2
  posterior_b <- prior_parameters$b + 1/2*t(y - x %*% prior_parameters$mu) %*% 
                              solve(diag(n) + x %*% V %*% t(x)) %*% (y - x %*% prior_parameters$mu)
  posterior_parameters <- list(a = posterior_a, b = posterior_b,
                               beta = posterior_beta, Sigma = posterior_Sigma)
  
  # Calculate posterior mean and variance
  posterior_mean_beta <- posterior_beta
  posterior_mean_sigma_sq <- posterior_b / (posterior_a - 1)
  posterior_mean <- list(beta = posterior_mean_beta, sigma_sq = posterior_mean_sigma_sq)
  
  posterior_var_beta <- (posterior_b / (posterior_a - 1)) %*% diag(posterior_Sigma)
  posterior_var_sigma_sq <- posterior_b^2 / ((posterior_a - 2) * (posterior_a - 1)^2)
  posterior_var <- list(beta = posterior_var_beta, sigma_sq = posterior_var_sigma_sq)
  
  # Calulate beta and sigma square credible intervals
  cred_intervals_beta_975 <- t(posterior_beta) - qt(.025, 2 * posterior_a, lower.tail = FALSE) *
                            sqrt((posterior_b / posterior_a) %*% diag(posterior_Sigma))
  cred_intervals_beta_025 <- t(posterior_beta) + qt(.025, 2 * posterior_a, lower.tail = FALSE) *
                            sqrt((posterior_b / posterior_a) %*% diag(posterior_Sigma))
  cred_intervals_beta_995 <- t(posterior_beta) - qt(.005, 2 * posterior_a, lower.tail = FALSE) *
                            sqrt((posterior_b / posterior_a) %*% diag(posterior_Sigma))
  cred_intervals_beta_005 <- t(posterior_beta) + qt(.005, 2 * posterior_a, lower.tail = FALSE) *
                            sqrt((posterior_b / posterior_a) %*% diag(posterior_Sigma))
  cred_intervals_beta <- list("97.5" = cred_intervals_beta_975, "0.25" = cred_intervals_beta_025,
                         "99.5" = cred_intervals_beta_995, "0.05" = cred_intervals_beta_005)
  
  cred_intervals_sigma_sq_975 <- 1 / qgamma(.025, shape = posterior_a, rate = posterior_b,
                                            lower.tail = FALSE)
  cred_intervals_sigma_sq_025 <- 1 / qgamma(.975, shape = posterior_a, rate = posterior_b,
                                            lower.tail = FALSE)
  cred_intervals_sigma_sq_995 <- 1 / qgamma(.005, shape = posterior_a, rate = posterior_b,
                                            lower.tail = FALSE)
  cred_intervals_sigma_sq_005 <- 1 / qgamma(.995, shape = posterior_a, rate = posterior_b,
                                            lower.tail = FALSE)
  cred_intervals_sigma_sq <- list("97.5" = cred_intervals_sigma_sq_975,
                                  "0.25" = cred_intervals_sigma_sq_025,
                                  "99.5" = cred_intervals_sigma_sq_995,
                                  "0.05" = cred_intervals_sigma_sq_005)
  cred_intervals = list(beta = cred_intervals_beta, sigma_sq = cred_intervals_sigma_sq)
  
  # Calculate the logarithm of the marginal likelihood of the data
  m_y <- dmvt(X = y,
              mu = x %*% prior_parameters$mu,
              sigma = (prior_parameters$b / prior_parameters$a) * (diag(n) + x %*% V %*% t(x)),
              df = 2 * prior_parameters$a,
              log = TRUE)
  
  results <- list(posterior_parameters = posterior_parameters, posterior_mean = posterior_mean,
                  posterior_var = posterior_var, credible_intervals = cred_intervals,
                  marginal_likelihood_y = m_y)
  
  return(results)
}

# 2a - Read the data, compute prior parameters, and fit the full model
data <- read.csv("~/Dropbox/MSc Data Science/Year2_2/Bayesian Statistics/Assignments/assign1/assignment_data_conjugateBayesRegression.txt", sep="")
y <- data$y
x <- matrix(unlist(data), ncol = 10, byrow = FALSE)
x[,1] <- 1
a <- 0.01
b <- 0.01
q <- ncol(x)
n <- nrow(x)
g <- n
mu <- replicate(q,0)
prior_parameters <- list(a = a, b = b, mu = mu)
results <- CalculatePosterior(y, x, prior_parameters, g)

# 2b - Calculate the posterior probabilities
P1 <- 1 - pmt(0, mean = results$posterior_parameters$beta[2],
              S = (results$posterior_parameters$b / results$posterior_parameters$a) *
                      results$posterior_parameters$Sigma[2,2],
              df = 2 * results$posterior_parameters$a)
P2 <- 1 - pmt(0, mean = results$posterior_parameters$beta[4],
              S = (results$posterior_parameters$b/results$posterior_parameters$a) * 
                      results$posterior_parameters$Sigma[4,4],
              df = 2 * results$posterior_parameters$a)
P3 <- pinvgamma(3, shape = results$posterior_parameters$a,
                    rate = results$posterior_parameters$b,
                    lower.tail = FALSE,
                    log.p = FALSE)

# 2c - Hypothesis testing using Bayes factor
x_new <- x[,1:3]
V <- g * solve(t(x_new) %*% x_new)
m0_y <- dmvt(X = y,
            mu = x_new %*% prior_parameters$mu[1:3],
            sigma = (prior_parameters$b / prior_parameters$a) * (diag(n) + x_new %*% V %*% t(x_new)),
            df = 2 * prior_parameters$a,
            log = TRUE)
bayes_factor <- 2 * (m0_y - results$marginal_likelihood_y)

