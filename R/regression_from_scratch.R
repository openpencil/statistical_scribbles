#### Objective: Walkthrough a linear regression estimation ####

#### Load library ####
library(Hmisc)
library(MASS)

#### Examine data of interest  ####
describe(Boston)

#### Define the variables ####
y <- Boston$medv
x <- as.matrix(Boston[-(which(names(Boston) == "medv"))])

#### Add an intercept term ####
intercept <- rep(1, nrow(x))
x <- cbind(intercept, x)

#### Solve for the regression equation ####
## betas represent the linear relationship between the dependent variable y and each of the independent variables
## Goal of linear regression: minimize the residual sum of squares (RSS) (Least squares estimate)
## i.e. find \beta_hat that minimizes the least squares estimate
## RSS = \sigma_{i = 1}^{n} (y_i - \beta*X_i)^2
## Minimize RSS w.r.t \beta_hat (i.e. find the differential w.r.t to \beta_hat)
## d(RSS)/d(\beta_hat) = 0
## Apply chain rule and simplify
## 0 = 2*\sigma_{i = 1}^{n} (y_i - \beta_hat*X_i) * (-X_i)
## 0 = \sigma_{i = 1}^{n} (y_i*X_i - \beta_hat*(X_i)^2)
## \sigma_{i = 1}^{n} (y_i*X_i) = \beta_hat * \sigma_{i = 1}^{n} (X_i)^2
## \beta_hat = \sigma_{i = 1}^{n} (y_i*X_i) /  \sigma_{i = 1}^{n} (X_i)^2
## Convert to matrix format:
## \sigma_{i = 1}^{n} (y_i*X_i) = y * X^T = y^T * X
## \sigma_{i = 1}^{n} (X_i)^2 = X^T*X
## Therefore, solution of interest: \beta_hat = (X^T*X)^-1 * X^T * y
## Relevant base functions:
## 1) matrix multiplication: %*%
## 2) matrix transpose: t()
## 3) inverse of an invertible matrix: solve()
betas <- solve(t(x) %*% x) %*% t(x) %*% y
betas_rounded <- round(x = betas, digits = 2)

#### Verify results with the lm() functions ####
lm_model <- lm(medv ~ ., data = Boston)

## Round for easier viewing
lm_betas <- round(lm_model$coefficients, 2)

## Create data.frame of results
results <- data.frame(calculated_betas = betas_rounded, lm_betas = lm_betas)

#### Calculate the confidence intervals around \beta_hat ####
## Residual degree of freedom
residual_df <- nrow(x) - ncol(x)
## RSS = \sigma_{i = 1}^{n} (y_i - \beta_hat*X_i)^2
rss <- sum((y - x %*% betas)^2)
## variance of RSS
res_var <- rss/residual_df

#### Verify results with summary.lm function ####
lm_summary <- summary.lm(object = lm_model)
lm_resvar <- (lm_summary$sigma)^2

## Deriving the expression for variance of beta_hat
## \beta_hat = (X^T*X)^-1 * X^T * y
## = (X^T*X)^-1 * X^T * (X*\beta + \epsilon)
## = (X^T*X)^-1 * (X^T*X*\beta + X^T*\epsilon)
## = (X^T*X)^-1) * X^T*X*\beta + (X^T*X)^-1) * X^T*\epsilon
## = \beta + (X^T*X)^-1) * X^T*\epsilon
## variance_covariances_beta_hat | X = E[(\beta_hat - \beta) * (\beta_hat - \beta)^T | X] (Because: A %*% A^T = A^2))
## = E[(\beta + (X^T*X)^-1) * X^T*\epsilon - \beta) * (\beta + (X^T*X)^-1) * X^T*\epsilon - \beta)^T | X]
## = E[((X^T*X)^-1 * X^T*\epsilon) * ((X^T*X)^-1 * X^T*\epsilon)^T | X]
## = E[((X^T*X)^-1 * X^T*\epsilon) * (X^T*\epsilon)^T * ((X^T*X)^-1)^T | X] (Because: (A %*% B)^T = B^T %*% A^T)
## = E[((X^T*X)^-1 * X^T*\epsilon) * \epsilon^T * X * ((X^T*X)^-1)^T | X]
## = E[((X^T*X)^-1 * X^T*\epsilon) * \epsilon^T * X * ((X^T*X)^T)^-1 | X] (Because: (A^-1)^T = (A^T)^-1)
## = E[((X^T*X)^-1 * X^T*\epsilon) * \epsilon^T * X * (X^T*X)^-1 | X] (Because: (A %*% B)^T = B^T %*% A^T)
## = (X^T*X)^-1 * X^T*E[\epsilon * \epsilon^T | X] * X * (X^T*X)^-1 (Moving expectation around random bits)
## = (X^T*X)^-1 * X^T*\sigma^2 * I * X * (X^T*X)^-1 (Because E[\epsilon) * \epsilon^T | X] is error variance)
## = \sigma^2 * I * (X^T*X)^-1
## (X^T*X)^-1
second_term <- solve(t(x) %*% x)
variance_covariances_beta_hat_given_x <- res_var * second_term
variance_beta_hat <- diag(variance_covariances_beta_hat_given_x)
cibeta_hat_lower <- betas - 1.96*sqrt(variance_beta_hat)
cibeta_hat_upper <- betas + 1.96*sqrt(variance_beta_hat)
confint_beta <- cbind(cibeta_hat_lower, cibeta_hat_upper)

#### Verify against confidence intervals from confint() function ####
lm_confint <- confint(object = lm_model)
