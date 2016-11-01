#' ### Understanding a single hypothesis test ###
#' ### An R-code-based walkthrough of concepts and intuition. ###

#' ##### Load libraries #####
library(ggplot2)
library(doMC)
library(Hmisc)
# source("https://bioconductor.org/biocLite.R")
# biocLite("limma")
library(limma)

#' #### SECTION I: The nuts and bolts a hypothesis test ####
#' ##### 1. Generate some expression data for a single gene #####
#' Seeds represent sampling variability
#' set.seed(149754)
#' set.seed(188549)
#' set.seed(2838)
#' generated data comes from normal distributions
#' here group1 and group2 are labels for the sample
#' simulated values represent gene-expression values
group1 <- log2(round(rnorm(n = 1000, mean = 500, sd = 50), 0))
group2 <- log2(round(rnorm(n = 1000, mean = 525, sd = 40), 0))

#' ##### 2. Visualize the distribution of data ######
#' Construct the dataframe for plotting
ggdata <- data.frame(Counts = c(group1, group2), stringsAsFactors = F)
#' Add a column with the group label
ggdata$group <- c(rep(x = "group1", times = length(group1)),
                  rep(x = "group2", times = length(group2)))


#' Call ggplot, define the dataframe for plotting, define the x and y axis
p <- ggplot(data = ggdata, mapping = aes(x = group, y = Counts))
#' Draw a boxplot
p <- p + geom_boxplot(outlier.size = 2.0, fatten = 0.5)
#' Or a violinplot
#' p <- p + geom_violin()
#' Remove the label of the x-axis
p <- p + xlab("")
#' Display the graph
p
#' Let's look a the density distribution of the generated data
p <- ggplot(data = ggdata, mapping = aes_string(x = "Counts", y = "..density.."))
# Draw density plots
p <- p + geom_density(aes(fill = group), alpha = 0.4)
# Draw the means
p <- p + geom_vline(xintercept = c(mean(group1), mean(group2)), linetype = "dashed",colour = "grey30")
p <- p + xlab("Value") + ylab("Density")
p

#' ##### 2. Describing the data #####
#' The five number summary
fivenumsummary <- fivenum(x = group1)
names(fivenumsummary) <- c("minimum", "lower-hinge", "median", "upper-hinge", "maximum")
fivenumsummary
#' Summary from Hmisc
describe(x = group1)
describe(x = group2)
#' Summary from base
summary(object = group1)
summary(object = group2)
#' Variance
var(x = group1)
var(x = group2)
#' Ratio fold change
fold_change_ratio <- mean(group1)/mean(group2)
#' Difference fold change
fold_change_difference <- mean(group1) - mean(group2)

#' ##### 3. Compute the test-statistics #####
#' Ordinary two sample t-statistic
#' Difference between mean values | Numerical change in average gene-expression
tnumerator <- mean(group1) - mean(group2)
#' Pooled standard deviation for a single gene | Standardize for noisy genes in which variance is more
#' so large absolute differences mean less.
#' Z-score: compute the number of standard deviations of the shift between treatment and control
tdenominator <- sqrt((sd(group1))^2/length(group1) + (sd(group2))^2/length(group2))
tstatistic <- tnumerator / tdenominator

#' Moderated t-statistic according to LIMMA
#' moderated_tstatistic <- tnumerator / tdenominator + some_constant
#' some_constant is chosen to shrink the tstatistic or minimize the coefficient of variation (sd/mean) of the tstatistic
#' As this constant is increased the resulting gene-ordering approaches that obtained with fold_change_difference
#' Shape the data into the weird format that LIMMA needs it to be
#' list of three elements:
#' 1) a dataframe with expression data of cases and controls | columns = cases/controls, rows = genes
#' expression_data <- transposed_data <- data.table(t(ggdata), keep.rownames = T)
#' setnames(x = expression_data, as.character(transposed_data[rn == "group"]))
#' # remove last row
#' expression_data <- expression_data[-nrow(expression_data)]
#' # remove first column
#' expression_data <- expression_data[,-1, with = F]
#' #' 2) a factor vector with group1/group2 designation of the samples
#' sample_class_factor <- factor(x = as.numeric(as.factor(colnames(expression_data)[1:ncol(expression_data)])),
#'                               levels = c(1, 2))
#' #' 3) a limma design matrix: a two column matrix with as many rows as samples
#' limma_design_matrix <- model.matrix(formula(~0 + sample_class_factor))
#' colnames(limma_design_matrix) <- c("group1", "group2")
#' # The contrast matrix generated below specifies the direction of the testing,
#' # e.g. group1 - group2. (these are columns of the design matrix)
#' limmacontrast <- makeContrasts(contrast = group1 - group2, levels = limma_design_matrix)
#' # fit a linear model: estimate the average "M-value" for each gene
#' expression_data <- rbind(expression_data, expression_data, expression_data)
#' limmafit <- lmFit(object = data.frame(expression_data), design = limma_design_matrix)

#' ##### 4. Thinking about relevance #####
#' What kind of expression differences between control and treatment have biological relevance?
#' - are large absolute changes relevant?
#' - are changes relative to underlying noise relevant?
#' What kind of fold-changes will noisy genes have?
#' What is more important: reproducibility or accuracy?
#' Are we more interested in reproducible measures of noise or differential expression?
#' Some simulation-based findings from Witten and Tibshirani
#' - fold_change_difference and tstatistic_moderated are modified tstatistic with different constants in the denominator
#' - unmodified tstatistic is never superior to moderated tstatistic in terms of accuracy or reproducibility
#' - reproducibility does not equate to accuracy | accuracy evaluates performance
#' - choice between fold-change and moderated t-statistic depends upon biological relevance


