#### Set work directory ####
workdir <- "~/mhst_correction/"
dir.create(workdir)
setwd(workdir)

#### Load libraries ####
library(ggplot2)
library(doMC)
library(plyr)

#### Generate some data ####
group1 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
group2 <- round(rnorm(n = 100, mean = 525, sd = 40), 0)

#### Visualize the distributions ####
ggdata <- data.frame(counts = c(group1, group2), stringsAsFactors = F)
ggdata$group <- c(rep(x = "group1", times = length(group1)), rep(x = "group2", times = length(group2)))
p <- ggplot(data = ggdata, mapping = aes(x = group, y = counts))
p <- p + geom_boxplot(outlier.size = 2.0, fatten = 0.5)

#### Generate entire dataset and conduct 100 hypotheses tests ####
no_real_difference <- unlist(mclapply(1:90, function(x){
  group1 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
  group2 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
  tstat <- t.test(x = group1, y = group2)
  return(tstat$p.value)
}, mc.cores = 4))
real_differences <- unlist(mclapply(1:10, function(x){
  group1 <- round(rnorm(n = 100, mean = 500, sd = 50), 0)
  group2 <- round(rnorm(n = 100, mean = 525, sd = 50), 0)
  tstat <- t.test(x = group1, y = group2)
  return(tstat$p.value)
}, mc.cores = 4))

## Put data together
pdata <- data.frame(pvalues = c(no_real_difference, real_differences), stringsAsFactors = F)
pdata$annot <- c(rep(x = "nodiff", times = length(no_real_difference)),
                 rep(x = "diff", times = length(real_differences)))
# apply bonferroni correction (control "Family-wise" error rate: probability of at least 1 type 1 error)
pdata$bonferroni_h0reject <- ifelse(pdata$pvalues <= 0.05/nrow(pdata), "rejecth0", "accepth0")
pdata$bonferroni_error <- ifelse(pdata$annot == "nodiff" & pdata$bonferroni_h0reject == "rejecth0",  "FP",
                                 ifelse(pdata$annot == "diff" & pdata$bonferroni_h0reject == "accepth0",  "FN",
                                        "Correct"))
table(pdata$bonferroni_error)

# apply Benjamini Hochberg threshold (control false discovery rate: expected proportion of type 1 errors)
# order p-values in increasing order
pdata <- pdata[order(pdata$pvalues),]
# test at which p-value is less than or equal to rank/length(p-values) * 0.05
pdata$BHrank <- seq_along(pdata$pvalues)
pdata$BH_threshold <- (pdata$BHrank/nrow(pdata))*0.05
BH_turningpoint <- max(which(pdata$pvalues <= pdata$BH_threshold))
pdata$BH_h0reject <- ifelse(pdata$BHrank <= BH_turningpoint, "rejecth0", "accepth0")

# apply Limma bloke's hacky BH correction (what the hell.)

# apply Storey's correction | qvalue (control positive false discovery rate: rate that the discoveries are false)



