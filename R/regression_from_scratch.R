#### Objective: Walkthrough regression estimation underlying lm() ####

#### Load library ####
library(Hmisc)

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
## i.e. find \beta that minimizes the least squares estimate
## RSS = \sigma_{i = 1}^{n} (y_i - \beta*X_i)^2
## Minimize RSS w.r.t \beta (i.e. find the differential w.r.t to \beta)
## d(RSS)/d(\beta) = 0
## 0 = 2*\sigma_{i = 1}^{n} (y_i - \beta*X_i) * (-X_i) (Applying chain rule)
## 0 = \sigma_{i = 1}^{n} (y_i*X_i - \beta*(X_i)^2)
## \sigma_{i = 1}^{n} (y_i*X_i) = \beta * \sigma_{i = 1}^{n} (X_i)^2
## \beta = \sigma_{i = 1}^{n} (y_i*X_i) /  \sigma_{i = 1}^{n} (X_i)^2
## Convert to matrix format:
## \sigma_{i = 1}^{n} (y_i*X_i) = y * X^T = y^T * X
## \sigma_{i = 1}^{n} (X_i)^2 = X^T*X
## Therefore, solution of interest: \beta = (X^T*X)^-1 * X^T * y
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
