#Figure 4 code
library(OlinkAnalyze)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(sva)
library(tibble)
library(tidyr)

setwd("") #set your file path
data <- readRDS("") #read in data 

data$time_specific_VAS_global <- ifelse(data$`Disease status`=="TN", data$month_0_physcore, data$month_6_physcore)
head(data$time_specific_VAS_global)
table(data$time_specific_VAS_global, data$`Disease status`)
table(data$time_specific_VAS_global, data$Site.x)
table(data$time_specific_VAS_global, data$common_ID) #common_ID equivalent to "SampleID"


da_lmer_alldata <- olink_lmer(df = data,
                              variable = "time_specific_VAS_global",
                              random = c("Site.x", "common_ID"), verbose = TRUE)  
View(da_lmer_alldata)

da_all_sig <- da_lmer_alldata %>%
  filter(Adjusted_pval <0.05) %>%
  pull(OlinkID)

da_lmer_alldata_posthoc <- olink_lmer_posthoc(df = data,
                                              olinkid_list = da_all_sig,
                                              variable = "time_specific_VAS_global",
                                              random = c("Site.x", "common_ID"),
                                              effect = "time_specific_VAS_global")
View(da_lmer_alldata_posthoc)
write.csv(da_lmer_alldata_posthoc, file = "da_lmer_all_data.csv")


#make point range plot of top 50 by greatest estimate
da_proteins <- da_lmer_alldata_posthoc
top_50_data <- da_proteins %>%
  mutate(abs_estimate = abs(estimate)) %>%  # Calculate absolute value of estimates
  group_by(Assay) %>%  # Group by Assay
  slice_max(order_by = abs_estimate, n = 1, with_ties = FALSE) %>%  # Take top 1 per Assay
  ungroup() %>%  # Ungroup for overall sorting
  arrange(desc(abs_estimate)) %>%  # Sort by absolute estimate
  slice(1:50) %>%  # Select the top 50 rows overall
  arrange(estimate) %>%  # Sort by actual estimate value
  mutate(Assay = factor(Assay, levels = rev(Assay)))  # Reverse factor levels
View(top_50_data)

ggplot(top_50_data, aes(x = Assay, y = estimate)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), color = "steelblue", size = 0.8) +
  coord_flip() +  # Flip coordinates for better readability
  labs(
    title = "",
    x = "Assay",
    y = "Estimate"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )
