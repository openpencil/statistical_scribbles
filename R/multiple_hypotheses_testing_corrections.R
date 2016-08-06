#### Load libraries ####
library(ggplot2)
library(doMC)

#### 1. Generate some data ####
group1 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
group2 <- round(rnorm(n = 100, mean = 525, sd = 40), 0)

#### 2. Visualize the distribution of data ####
# construct the dataframe for plotting
ggdata <- data.frame(Counts = c(group1, group2), stringsAsFactors = F)
# add a column with the group label
ggdata$group <- c(rep(x = "group1", times = length(group1)),
                  rep(x = "group2", times = length(group2)))
# call ggplot, define the dataframe for plotting, define the x and y axis
p <- ggplot(data = ggdata, mapping = aes(x = group, y = Counts))
# draw a boxplot
p <- p + geom_boxplot(outlier.size = 2.0, fatten = 0.5)
# remove the label of the x-axis
p <- p + xlab("")
# display the graph
p
# STOP: What do you observe?
# THINK: Do you need formal statistics to establish differences?
# THINK ALOUD: Why or Why not?

#### 3. Generate entire dataset and conduct 100 hypotheses tests ####
# 90 measurements with no differences in the underlying distribution.
no_real_difference <- unlist(mclapply(1:90, function(x){
  group1 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
  group2 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
  tstat <- t.test(x = group1, y = group2)
  return(tstat$p.value)
}, mc.cores = 4))

# 10 measurements with differences in the mean of the distribution they are drawn from.
real_differences <- unlist(mclapply(1:10, function(x){
  group1 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
  group2 <- round(rnorm(n = 100, mean = 525, sd = 50), 0)
  tstat <- t.test(x = group1, y = group2)
  return(tstat$p.value)
}, mc.cores = 4))

## Put all the p-values together
pdata <- data.frame(pvalues = c(no_real_difference, real_differences), stringsAsFactors = F)
pdata$annot <- c(rep(x = "nodiff", times = length(no_real_difference)),
                 rep(x = "diff", times = length(real_differences)))

## Visualize the distribution of values
#  display_density_graph <- function(datasub, name_of_column_to_plot){
#   #'
#   #' @param datasub dataset which contains the values to be plotted
#   #' @param name_of_column_to_plot
#   #'
#   #' @output will give you a graph object that you can either display or save.
#   #'
#   # call ggplot, define the dataframe for plotting, define the x and y axis
#   p <- ggplot(data = datasub, mapping = aes_string(x = name_of_column_to_plot, y = "..density.."))
#   # draw a density plot
#   p <- p + geom_density(alpha = 0.4, fill = "#80d827", colour = "#d8f3bd")
#   p <- p + geom_vline(xintercept = 0.05, linetype = "dashed",colour = "grey30")
#   return(p)
# }


#### 4. Examine errors with no multiple hypothesis testing correction ####
pdata$nocorrection_h0reject <- ifelse(pdata$pvalues <= 0.05, "rejecth0", "accepth0")
# if there was no difference between the groups (i.e. null hypothesis was correct)
# and you rejected the null = False positive
pdata$nocorrection_error <- ifelse(pdata$annot == "nodiff" & pdata$nocorrection_h0reject == "rejecth0",  "FP",
                                   # if there was a difference between the groups and you didn't detect this
                                   # (i.e you accepted the null that there was no difference - the null): False negative
                                   ifelse(pdata$annot == "diff" & pdata$nocorrection_h0reject == "accepth0",  "FN", "Correct"))
table(pdata$nocorrection_error)
# STOP: What do you observe?
# THINK: Given that you know the truth, is this result acceptable to you?
# THINK AGAIN: How would your answer change if you had not known the truth?
# THINK ALOUD: Can you think of examples where you would accept this result and where you would not?

#### 4. The Bonferroni correction ####
# apply Bonferroni correction (control "Family-wise" error rate: probability of at least 1 type 1 error)
pdata$bonferroni_h0reject <- ifelse(pdata$pvalues <= 0.05/nrow(pdata), "rejecth0", "accepth0")
pdata$bonferroni_error <- ifelse(pdata$annot == "nodiff" & pdata$bonferroni_h0reject == "rejecth0",  "FP",
                                 ifelse(pdata$annot == "diff" & pdata$bonferroni_h0reject == "accepth0",  "FN",
                                        "Correct"))
table(pdata$bonferroni_error)
# STOP: How do results with the bonferroni correction compare to the results with no corrections?
# THINK: Is this result acceptable to you?
# THINK ALOUD: Why or Why not?
# THINK AGAIN: How does knowing and not-knowing the truth about the data change your perception of these results?

#### 4. The Benjamini Hochberg threshold ####
# apply Benjamini Hochberg threshold (control false discovery rate: expected proportion of type 1 errors)
# order p-values in increasing order (i.e. lowest p-value is at the top)
pdata <- pdata[order(pdata$pvalues),]
# assign ranks to the p-values (i.e. lowest p-value gets the lowest rank)
pdata$BHrank <- seq_along(pdata$pvalues)
# calculate the Benjamini-Hochberg threshold: rank/length(p-values) * 0.05
# lowest p-value will be compared to the lowest threshold.
# highest p-value with be compared to the most lenient (highest) threshold (rank/length(p-values) * 0.050 = 1 * 0.05 = 0.05)
pdata$BH_threshold <- (pdata$BHrank/nrow(pdata))*0.05
# test at which rank, p-values are less than or equal to rank/length(p-values) * 0.05
BH_turningpoint <- max(which(pdata$pvalues <= pdata$BH_threshold))
pdata$BH_h0reject <- ifelse(pdata$BHrank <= BH_turningpoint, "rejecth0", "accepth0")
pdata$BH_error <- ifelse(pdata$annot == "nodiff" & pdata$BH_h0reject == "rejecth0",  "FP",
                         ifelse(pdata$annot == "diff" & pdata$BH_h0reject == "accepth0",  "FN", "Correct"))
table(pdata$BH_error)
# STOP: How do results employing the BH threshold compare to the uncorrected 0.05 threshold?
# THINK: Is this result acceptable to you?
# THINK ALOUD: Why or Why not?
# THINK AGAIN: Can you talk about examples where you would pick this method rather than Bonferroni or no corrections at all?


# Note: p-value is compared to rank/n * 0.05
# p-value ~ rank/n * 0.05
# Twist this around:
# p-value * n/rank ~ 0.05
# you could change the p-value and always compare it to 0.05 (i.e. Smyth's hack to avoid thinking about the adjustment.)
pdata$smyth_hackyBH_adjustment <- pdata$pvalues * (nrow(pdata)/pdata$BHrank)
# LOOK at the values.
# They are not monotonic.
# SO:
# for the highest p-value, use formula from L117, i.e. compare it to 0.05
pvalue_decreasing_rank <- order(pdata$pvalues, decreasing = T)
# create an empty array
container_for_padjustment <- rep(0, nrow(pdata))
for (rank in pvalue_decreasing_rank) {
  #' rank <- 100
  pval <- pdata$pvalues[rank]
  if (rank == length(pvalue_decreasing_rank)) {
    # highest p-value is the original p-value, unchanged
    container_for_padjustment[rank] <- pval
  } else {
    # for the next highest, we will use formula which is min(L117, and the one bigger than L117)
    container_for_padjustment[rank] <- min(pval * (nrow(pdata))/rank, container_for_padjustment[rank + 1])
  }
}
pdata$smyth_hackyBH_adjustment <- container_for_padjustment
pdata$padjust_check <- p.adjust(p = pdata$pvalues, method = "BH")
pdata[which(abs(pdata$smyth_hackyBH_adjustment[order(pdata$smyth_hackyBH_adjustment)] - pdata$padjust_check) > 0.01), ]


#### 5. Smyth's hacky "BH" p-value adjustment ####
# o <- order(p, decreasing = TRUE)
pdata$smyth_rank <- order(pdata$pvalues, decreasing = TRUE)
# ro <- order(o) (order of the p-values)
ro <- order(pdata$smyth_rank)
# i <- lp:1L
i <- (nrow(pdata):1L)
po <- pdata$pvalues[pdata$smyth_rank]
n <- nrow(pdata)
# (parallel minima) | compare first value to 1?
# pmin(1, cummin(n/i * p[o]))[ro]
pdata$smyth_bhadjustment <- pmin(1, cummin(n/i * po))[ro]
pdata$padjust_check <- p.adjust(p = pdata$pvalues, method = "BH")
which(pdata$smyth_bhadjustment != pdata$padjust_check, arr.ind = T)

pdata$BH_h0reject <- ifelse(pdata$BHrank <= BH_turningpoint, "rejecth0", "accepth0")
pdata$BH_error <- ifelse(pdata$annot == "nodiff" & pdata$BH_h0reject == "rejecth0",  "FP",
                         ifelse(pdata$annot == "diff" & pdata$BH_h0reject == "accepth0",  "FN", "Correct"))
table(pdata$BH_error)

# apply Storey's correction | qvalue (control positive false discovery rate: rate that the discoveries are false)



