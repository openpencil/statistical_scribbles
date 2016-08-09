#' ## Multiple hypothesis testing and correcting for the resulting errors. ##
#' ### A R-code-based walkthrough of concepts and intuition. ###

#' ##### Load libraries #####
library(ggplot2)
library(doMC)

#' ##### Functions #####
visualize_distribution <- function(datasub, name_of_column_to_plot, fillwhat = "blue"){
  #'
  #' @param datasub dataset which contains the values to be plotted
  #' @param name_of_column_to_plot
  #'
  #' @output will give you a graph object that you can either display or save.
  #'
  #' call ggplot, define the dataframe for plotting, define the x and y axis
  p <- ggplot(data = datasub, mapping = aes_string(x = name_of_column_to_plot, y = "..density.."))
  #' draw a density plot
  # p <- p + geom_density(alpha = 0.4, fill = "#80d827", colour = "#d8f3bd")
  if (fillwhat == "blue"){
    p <- p + geom_histogram(bins =  50, fill = "#35aeba", colour = "#d8f3bd")
  } else {
    p <- p + geom_histogram(bins =  50, aes(fill = group))
  }
  p <- p + geom_vline(xintercept = 0.05, linetype = "dashed",colour = "grey30")
  return(p)
  }

#' #### SECTION I: What is hypothesis testing? ####

#' ##### 1. Generate some data #####
#' Seeds introduce variability ##

#' set.seed(149754)

#' set.seed(188549)

#' set.seed(2838)

group1 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
group2 <- round(rnorm(n = 100, mean = 525, sd = 40), 0)

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
#' Remove the label of the x-axis
p <- p + xlab("")
#' Display the graph
p
#' Also display the two count distributions as histograms and colour by group
visualize_distribution(datasub = ggdata, name_of_column_to_plot = "Counts", fillwhat = "group")

#' ##### Stop and reflect #####
#' What do you observe?
#' What happens when you regenerate the data with different random number seed?
#' THINK: Do you need formal statistics to establish differences? Why or Why not?
#' What is the null hypothesis?
#' What is the alternative hypothesis?

#' #### SECTION II: What is multiple hypothesis testing? ####

#' ##### 1. Generate entire dataset and conduct 100 hypotheses tests #####
#' 90% measurements with no differences in the underlying distribution.
no_real_difference <- unlist(mclapply(1:9000, function(x){
  #' Run once
  group1 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
  group2 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
  tstat <- t.test(x = group1, y = group2)
  #' tstat
  return(tstat$p.value)
}, mc.cores = 4))

#' 10% measurements with differences in the mean of the distribution they are drawn from.
real_differences <- unlist(mclapply(1:1000, function(x){
  #' Run once
  group1 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
  group2 <- round(rnorm(n = 100, mean = 525, sd = 50), 0)
  tstat <- t.test(x = group1, y = group2)
  #' tstat
  return(tstat$p.value)
}, mc.cores = 4))

#' Put all the p-values together
pdata <- data.frame(pvalues = c(no_real_difference, real_differences), stringsAsFactors = F)
pdata$annot <- c(rep(x = "nodiff", times = length(no_real_difference)),
                 rep(x = "diff", times = length(real_differences)))
#' ##### Stop and check #####
dim(pdata)
head(pdata)
table(pdata$annot)

#' ##### 2. Examine the errors after multiple hypothesis testing #####
#' We want only 0.05 errors per test. (i.e. <1 errors per test).
#' BUT when we do 10K tests, these errors compound to 0.05 * 10000 = 500 errors.
#' Let's look at the data to see how many errors we netted by doing 10000 tests.

pdata$nocorrection_h0reject <- ifelse(pdata$pvalues <= 0.05, "rejecth0", "accepth0")

pdata$nocorrection_error <- ifelse(pdata$annot == "nodiff" & pdata$nocorrection_h0reject == "rejecth0",
                                   "FalsePositives",
                                   #' If there were no differences between the groups (i.e. null hypothesis was correct)
                                   #' and you rejected the null
                                   #' False positives
                                   ifelse(pdata$annot == "diff" & pdata$nocorrection_h0reject == "accepth0",
                                          "FalseNegatives",
                                          # If there was a difference between the groups and you didn't detect this
                                          # (i.e. you accepted the null that there was no difference - the null)
                                          # False negatives
                                          "Correct")
                                          # If it is neither a false positive or a false negative, it is a correct call.
                                          )

table(pdata$nocorrection_error)

#' ##### Stop and reflect #####
#' How many errors do you observe?
#' Given that you know the truth, is this result acceptable to you?
#' How would your answer change if you had not known the truth?
#' THINK ALOUD: Can you think of examples where you would accept this result and where you would not?

#'### SECTION III: What is multiple hypothesis testing correction? ####
#' We will explore 3 most popular ways of correcting some of the errors
#' and achieving control of Type I errors.

#'### 1. The Bonferroni correction ####
#' Described by Olive Jean Dunn in 1959 and 1961.
#' Bonferroni formulated the Bonferroni inequality which helped Dunn formulate the correction.
#' We want only 0.05 errors per test. (i.e. <1 errors per test).
#' BUT when we do 10K tests, these errors compound to 0.05 * 10000 = 500 errors.
#' How do we drive back the error rate to 0.05 per 10000 tests?
#' Divide the 0.05 significance level by the number of tests.
#' Controls the "Family-wise" error rate.

pdata$bonferroni_h0reject <- ifelse(pdata$pvalues <= 0.05/nrow(pdata), "rejecth0", "accepth0")
pdata$bonferroni_error <- ifelse(pdata$annot == "nodiff" & pdata$bonferroni_h0reject == "rejecth0",
                                 "FalsePositives",
                                 ifelse(pdata$annot == "diff" & pdata$bonferroni_h0reject == "accepth0",
                                        "FalseNegatives",
                                        "Correct"))
table(pdata$bonferroni_error)
#' Look at the error table without corrections:
table(pdata$nocorrection_error)
#' STOP: How do results with the bonferroni correction compare to the results with no corrections?
#' THINK: Is this result acceptable to you?
#' THINK ALOUD: What has happened to power - i.e. True positive rate.



#'### 2a. The Benjamini Hochberg threshold ####
#' Was presented in a paper in 1995, inspired by substantial body of works from 1979 onwards.
#' It took more than 20 years for theory behind the threshold to be formalized and gain acceptance.
#' For context see: https://en.wikipedia.org/wiki/False_discovery_rate#Literature

#' Can we improve power of our tests? We can tolerate some errors in return.
#' So instead of comparing all significant p-values to the equally strict
#' 0.05/(number of tests), can we be most strict for lowest p-values and
#' least strict for the higher p-values, with varying degrees of strictness
#' in between? In other words, can we apply a more nuanced Bonferroni correction?

#' order p-values in increasing order (i.e. lowest p-value is at the top)
pdata <- pdata[order(pdata$pvalues),]
#' assign ranks to the p-values (i.e. lowest p-value gets the lowest rank)
pdata$BHrank <- seq_along(pdata$pvalues)
#' calculate the Benjamini-Hochberg threshold: rank/length(p-values) * 0.05
#' lowest p-value will be compared to the lowest threshold, i.e. the harshest threshold
#' highest p-value with be compared to the most lenient (highest) threshold (rank/length(p-values) * 0.050 = 1 * 0.05 = 0.05)
pdata$BH_threshold <- (0.05/nrow(pdata))*(pdata$BHrank)
#' test at which rank, p-values are less than or equal to rank/length(p-values) * 0.05
BH_turningpoint <- max(which(pdata$pvalues <= pdata$BH_threshold))
pdata$BH_h0reject <- ifelse(pdata$BHrank <= BH_turningpoint, "rejecth0", "accepth0")
pdata$BH_error <- ifelse(pdata$annot == "nodiff" & pdata$BH_h0reject == "rejecth0",
                         "FalsePositive",
                         ifelse(pdata$annot == "diff" & pdata$BH_h0reject == "accepth0",
                                "FalseNegative",
                                "Correct"))
table(pdata$BH_error)
#' Look at the error table without corrections and with the Bonferroni corrections:
table(pdata$bonferroni_error)
table(pdata$nocorrection_error)

#' STOP: How do results employing the BH threshold compare to the Bonferroni and the uncorrected 0.05 threshold?
#' THINK: Is this result acceptable to you?
#' THINK ALOUD: Can you think about examples where you would pick this method rather than Bonferroni or no corrections at all?


#'### 2b. Smyth's interpretation of BH threshold + StoreyTibs monotonic correction: "The BH adjusted p-value" ####
#' Smyth (of LIMMA fame) took the BH idea and blended it with Storey and Tibshirani's method of
#' achieving monotonicity to contribute a function called "p.adjust()" to R in 2002.

#' In the BH threshold, the p-value is compared to (0.05/n) * rank
#' i.e. p-value ~ (0.05/n) * rank
#' Twist this around:
#' p-value * n/rank ~ 0.05 --- (Equation:Smythhack)
#' So by applying Smythhack you could change the p-value and
#' always compare it to 0.05. In other words, Smyth's hack made
#' it possible for people forget about the problems of the
#' 0.05 significance level in multiple testing.

#' Order the original p-values in decreasing order, i.e. highest p-value at the top.
pdata <- pdata[order(pdata$pvalues, decreasing = T),]
#' assign ranks to the p-values (i.e. highest p-value gets highest rank)
pdata$Smythrank <- rank(pdata$pvalues, ties.method = "max")
#' On these ordered p-values, we will now apply the Smythhack equation
pdata$Smythhack <- pdata$pvalues * (nrow(pdata)/pdata$Smythrank)
#' LOOK at the values. Specifically, look at the order of the values
head(pdata$Smythhack, 50)
head(rank(pdata$Smythhack, ties.method = "max"), 50)
#' Compare these to the rank of the original p-values:
head(pdata$Smythrank, 50)
#' What is wrong here?
#' Even through the original p-values were ordered, the SmythHack results are not ordered (i.e. not monotonic)
#' Why is this bad?
#' You cannot set a threshold and call every value below it significant.
#' This is an important property of the original BH threshold - and this
#' non-monotonicity is why BH formulated the threshold and not an
#' adjustment. So Smyth looked around and found a method proposed
#' by Storey and Tibshirani the year before to bring monotonicity
#' to the non-monotonic values. He coded this monotonic conversion
#' into the p.adjust() function and contributed it to R base.
#'
#' Apply Storey method of monotonicity to SmythHack
pdata$Smythhack_monotonic <- rep(x = 0, times = nrow(pdata))
for (rank in pdata$Smythrank) {
  if (rank == length(pdata$pvalues)) {
    # for the highest ranked Smythhack value, keep the Smythhack value
    pdata$Smythhack_monotonic[which(pdata$Smythrank == rank)] <- pdata$Smythhack[which(pdata$Smythrank == rank)]
  } else {
    # from the second highest ranked Smythhack value onwards
    # compare the Smythhack value to the Smythhack value that is ranked one rank above
    # assign the Smythhack value to the lesser of the two Smythhack values
    pdata$Smythhack_monotonic[which(pdata$Smythrank == rank)] <-
      min(pdata$Smythhack[which(pdata$Smythrank == rank)], pdata$Smythhack_monotonic[which(pdata$Smythrank == (rank + 1))])
  }
}
#' SANITY CHECK: Compare what we implemented with the p.adjust() function.
padjust_check <- p.adjust(p = pdata$pvalues, method = "BH")
pdata[which(abs(pdata$Smythhack_monotonic - pdata$padjust_check) > 0.01), ]

#' Compute the FP, FN and correct calls.
pdata$smythhack_h0reject <- ifelse(pdata$Smythhack_monotonic <= 0.05, "rejecth0", "accepth0")
pdata$smythhack_error <- ifelse(pdata$annot == "nodiff" & pdata$smythhack_h0reject == "rejecth0",
                                "FalsePositive",
                                 ifelse(pdata$annot == "diff" & pdata$smythhack_h0reject == "accepth0",
                                        "FalseNegative",
                                        "Correct"))
table(pdata$smythhack_error)
#' Compare with the earlier error tables:
table(pdata$BH_error)
table(pdata$bonferroni_error)
table(pdata$nocorrection_error)
#' OBSERVE: Is there a difference between the results from the SmythHack and from applying the BH threshold?
#' THINK ALOUD: Why use the Smythhack adjustment instead of the original BH threshold?


#'### 3. Storey's Q-value ####
#' In 2001, Storey and Tibshirani (Storey's gradschool mentor) took a more empirical viewpoint
#' and extract some more information from the data. Instead of explicitly minimizing the
#' errors and maximizing the power, they determined the false discovery rate, at each p-value,
#' if that p-value were the significance level of the test. The minimum FDR at each p-value
#' is called its q-value.

#' Before we explore the q-value, we need to think about some concepts.
#'
#' CONCEPT I: What does the distribution of original p-values look like?
#' Visualize the distribution of values
visualize_distribution(datasub = pdata, name_of_column_to_plot = "pvalues")
#' DISCUSS: What kind of distribution is this? Which are null p-values? Which are the significant ones?
#' Try increasing or decreasing the sample sizes in the mclapply loop and re-run the visualization above.
#' Does the shape of the distribution change?
#' We arrive at an empirical truth: Truly null p-values are _______ly distributed.
#' So if I ask you how many null p-values occur below say, 0.75, how will you determine this?

#' CONCEPT II: What is the False Discovery Rate?
#' False discovery rate (FDR) is the proportion of discoveries that are false.
#' i.e. FDR = (number of false discoveries) / (total number of discoveries)
#' For any threshold (such as 0.05):
#' FDR = (number of null p-values < threshold) / (total number of p-values < threshold)
#' Total number of null p-values = Proportion of null p-values (P0) in data * total number of p-values
#' Since the distribution of null p-values is ________ (fill in the answer here from the observation exercise above),
#' the number of null p-values at a threshold is Total number of null p- values * threshold.
#' i.e. FDR = P0 * length(pdata$pvalues) * threshold / length(which(pdata$pvalues < threshold))

#' CONCEPT III. What is a q-value?
#' The Q-value of any given p-value is the FDR if the given p-value were used as the threshold for calling significance.
#' For example: The Q-value of a p-value = 0.6 is the FDR if we called all p-values below 0.6 as significant.

#' So, to compute the Q-value for every p-value (or the FDR at every p-value), we need 4 quantities:
#' 1) P0: The proportion of null p-values in the data.
#' 2) length(pdata$pvalues) or the number of tests
#' 3) threshold: the p-value itself
#' 4) length(which(pdata$pvalues < threshold)): number of significant finding at that threshold
#' and plug into the formula:
#' Qvalue = FDR_at_a_given_pvalue = P0 * length(pdata$pvalues) * given_pvalue / length(which(pdata$pvalues < given_pvalue))

#' So how do we estimate P0: of the proportion of null p-values in the data?
#'###  STEP 1: Estimate P0, the proportion of truly null p-values ####
#' Let's first think about it and then compute it.
#' recall the distribution of p-values
p <- visualize_distribution(datasub = pdata, name_of_column_to_plot = "pvalues")
#' let's draw a bunch of thresholds on the distribution
sequence_of_thresholds <- seq(from = 0, to = 0.95, by = 0.05)
p <- p + geom_vline(xintercept = sequence_of_thresholds, linetype = "dashed",colour = "grey30")
p
#' What is the proportion of true null pvalues mixed in with the truly non-null values at each threshold?
#' As the threshold approaches 1, the proportion of true nulls will become more stable.

proportion_true_nulls_at_each_threshold <- sapply(sequence_of_thresholds, function(threshold_value){
  #' threshold_value <- 0.85
  #' Number of observed p-values greater than the threshold / Total number of p-values greater than the threshold
  number_true_nulls <- sum(pdata$pvalues >= threshold_value)
  total_pvalues_greater_than_threshold <- (length(pdata$pvalues) * (1 - threshold_value))
  prop_true_nulls <- number_true_nulls/total_pvalues_greater_than_threshold
  return(prop_true_nulls)
})

#' Let's plot the distribution of true nulls
qplot(x = sequence_of_thresholds, y = proportion_true_nulls_at_each_threshold)
#' What does it look like? Non-smooth. Non-linear. Non-monotonic
#' So the question is : At what threshold do you get the least variability or the most stability in the true
#' proportion of nulls? In other words, what threshold will lead to a good estimator of the true proportion
#' of nulls? To determine this, we build a spline object or in other words, fit a curve through these points.
cubic_spline_object <- smooth.spline(x = sequence_of_thresholds,
                                     y = proportion_true_nulls_at_each_threshold,
                                     df = 3)
cubic_spline_predicted_pvalue <- predict(cubic_spline_object, x = sequence_of_thresholds)$y
#' Look at the p-value predicted by the cubic spline, i.e. the smoothed p-value
qplot(x = sequence_of_thresholds, y = cubic_spline_predicted_pvalue)
# As the threshold approaches 1, the proportion of true nulls becomes more stable so pick the minimum prediction
# Since the proportion cannot exceed 1, we take a minimum of the prediction of the largest threshold and 1.
estimate_of_true_proportion_of_nulls <- min(cubic_spline_predicted_pvalue[length(sequence_of_thresholds)], 1)

#'#### STEP 2: Compute the Q-value or the FDR at every p-value ####
#' order p-values in decreasing order (i.e. highest p-value is at the top)
pdata <- pdata[order(pdata$pvalues, decreasing = T),]
#' FDR = P0 * length(pdata$pvalues) * threshold / length(which(pdata$pvalues < threshold))
pdata$Storeyrank <- rank(pdata$pvalues, ties.method = "max")
#' THINK: pdata$Storeyrank of a given p-value is basically the number of pvalues less than the given pvalue
pdata$qvalues <- (estimate_of_true_proportion_of_nulls * length(pdata$pvalues) * pdata$pvalues)/ pdata$Storeyrank

#' Look at the order of the computed qvalues
head(pdata$qvalues, 50)
head(order(pdata$qvalues), 50)
#' Not monotonic, even though we ordered the p-values before computing the qvalue!
#' So next we make them monotonic.| i.e. the method that Storey and Tibshirani used to make qvalues monotonic
#' and the one that Smyth copied in p.adjust()
#
#'#### STEP 3: Make Q-values monotonic ####
#' Compare each qvalue to the qvalue above.
#' Assign the qvalue to the minimum of the two.

pdata$qvalues_monotonic <- rep(x = 0, times = nrow(pdata))
for (rank in pdata$Storeyrank) {
  if (rank == length(pdata$pvalues)) {
    #' for the highest ranked q-value, keep the qvalue
    pdata$qvalues_monotonic[which(pdata$Storeyrank == rank)] <- pdata$qvalues[which(pdata$Storeyrank == rank)]
  } else {
    #' from the second highest ranked q-value onwards
    #' compare the q-value to the q-value that is ranked one rank above
    #' assign the q-value to the lesser of the two q-values
    pdata$qvalues_monotonic[which(pdata$Storeyrank == rank)] <-
      min(pdata$qvalues[which(pdata$Storeyrank == rank)], pdata$qvalues_monotonic[which(pdata$Storeyrank == (rank + 1))])
  }
}

#' SANITY CHECK: Compare what we implemented with the qvalue function
#' source("https://bioconductor.org/biocLite.R")
#' biocLite("qvalue")
#' library(qvalue)
pdata$qvalues_from_package <- qvalue::qvalue(p = pdata$pvalues)$qvalues
pdata[which(abs(pdata$qvalues_monotonic - pdata$qvalues_from_package) > 0.01), ]

#' Compute the FP, FN and correct calls.
pdata$qvalue_h0reject <- ifelse(pdata$qvalues_monotonic <= 0.05, "rejecth0", "accepth0")
pdata$qvalue_error <- ifelse(pdata$annot == "nodiff" & pdata$qvalue_h0reject == "rejecth0",
                                "FalsePositive",
                                ifelse(pdata$annot == "diff" & pdata$qvalue_h0reject == "accepth0",
                                       "FalseNegative",
                                       "Correct"))
table(pdata$qvalue_error)
#' Compare with the earlier error tables:
table(pdata$smythhack_error)
table(pdata$BH_error)
table(pdata$bonferroni_error)
table(pdata$nocorrection_error)
#' THINK ALOUD: When would you use qvalue rather than any other form of correction?
#' THINK ALOUD: What can happen when you don't have enough tests and data?
