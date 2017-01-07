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
## Solution of interest: beta = (X^T X)^-1 * X^T * y
## Relevant base functions:
## 1) matrix multiplication: %*%
## 2) matrix transpose: t()
## 3) inverse of an invertible matrix: solve()
##
## betas embody the linear relationship between the dependent variable y and each of the independent variables
betas <- solve(t(x) %*% x) %*% t(x) %*% y
betas_rounded <- round(x = betas, digits = 2)

#### Verify results with the lm() functions ####
lm_model <- lm(medv ~ ., data = Boston)

## Round for easier viewing
lm_betas <- round(lm_model$coefficients, 2)

## Create data.frame of results
results <- data.frame(calculated_betas = betas_rounded, lm_betas = lm_betas)


