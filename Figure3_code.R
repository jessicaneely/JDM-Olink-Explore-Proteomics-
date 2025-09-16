#Figure 3 code

library(dplyr)
library(ggplot2)

#read in BL v CTL results and FU v CTL results
setwd("") #replace with your directory
bl_proteins <- read.csv(file ="lmer_tn_hc_posthoc.csv")
fu_proteins <- read.csv(file = "FU_ctl_posthoc.csv")
View(bl_proteins)

bl_proteins_log <- bl_proteins[ , c("Assay", "estimate..difference.in.NPX.value.between.contrast.terms.", "Adjusted_pval", "Threshold", "significant_lmer_no_posthoc")]
head(bl_proteins_log)
bl_proteins_log$sig_overall <- ifelse(bl_proteins_log$Threshold=="Significant" & bl_proteins_log$significant_lmer_no_posthoc=="Significant", "sig", "ns")
table(bl_proteins_log$sig_overall)

bl_proteins_log <- bl_proteins_log %>%
  rename(NPX_diff_BL = estimate..difference.in.NPX.value.between.contrast.terms.)
bl_proteins_log$sig_direction <- ifelse(bl_proteins_log$NPX_diff_BL>0 & bl_proteins_log$sig_overall=="sig", "sigUp", ifelse(
  bl_proteins_log$NPX_diff_BL<0 & bl_proteins_log$sig_overall=="sig", "sigDown", "ns"
))
table(bl_proteins_log$sig_direction)

#look at overlapping proteins for this
fu_proteins_log <- fu_proteins[ , c("Assay", "estimate..difference.in.NPX.value.between.contrast.terms.", "Adjusted_pval", "Threshold", "significant_lmer_no_posthoc")]
head(fu_proteins_log)
fu_proteins_log$sig_overall <- ifelse(fu_proteins_log$Threshold=="Significant" & fu_proteins_log$significant_lmer_no_posthoc=="Significant", "sig", "ns")
table(fu_proteins_log$sig_overall)

fu_proteins_log <- fu_proteins_log %>%
  rename(NPX_diff_fu = estimate..difference.in.NPX.value.between.contrast.terms.)
fu_proteins_log$sig_direction <- ifelse(fu_proteins_log$NPX_diff_fu>0 & fu_proteins_log$sig_overall=="sig", "sigUp", ifelse(
  fu_proteins_log$NPX_diff_fu<0 & fu_proteins_log$sig_overall=="sig", "sigDown", "ns"
))
table(fu_proteins_log$sig_direction)

#take only intersection of proteins significant in both to see if direction is same or different
bl_proteins_sig <- bl_proteins_log[bl_proteins_log$sig_overall=="sig", ]$Assay
fu_proteins_sig <- fu_proteins_log[fu_proteins_log$sig_overall=="sig", ]$Assay
common_sig <- intersect(bl_proteins_sig, fu_proteins_sig)

#new dataframe to merge NPX bl and NPX fu
bl_common <- bl_proteins_log[bl_proteins_log$Assay %in% common_sig, ]
bl_filtered <- bl_common %>%
  group_by(Assay) %>%
  slice_max(order_by = NPX_diff_BL, n = 1) %>%
  ungroup()

fu_common <- fu_proteins_log[fu_proteins_log$Assay %in% common_sig, ]
fu_filtered <- fu_common %>%
  group_by(Assay) %>%
  slice_max(order_by = NPX_diff_fu, n = 1) %>%
  ungroup()

#need to remove duplicates before merging

common_df <- merge(bl_filtered, fu_filtered, by = "Assay")
head(common_df)

#define significance direction concordance
common_dir <- common_df %>%
  mutate(
    Significance = case_when(
      NPX_diff_BL>0 & NPX_diff_fu>0 ~ "Up in Both",
      NPX_diff_BL<0 & NPX_diff_fu<0 ~ "Down in Both",
      NPX_diff_BL>0 & NPX_diff_fu<0 ~ "Opposite",
      NPX_diff_BL<0 & NPX_diff_fu>0 ~ "Opposite",
    )
  )
head(common_dir)

# Define colors for significance categories
sig_colors <- c("Up in Both" = "red", "Down in Both" = "blue", "Opposite" = "purple")

# Create the log-log scatter plot
ggplot(common_dir, aes(x = NPX_diff_BL, y = NPX_diff_fu, color = Significance)) +
  geom_point(alpha = 0.7) +  # Scatter points
  scale_color_manual(values = sig_colors) +  # Use defined colors
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal reference
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical reference
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "gray40") +  # y = x line
  geom_text_repel(data = common_dir %>% filter(abs(NPX_diff_BL) > 1 | abs(NPX_diff_fu) > 1), 
                  aes(label = Assay), size = 5, box.padding = 0.5) +  # Highlight extreme proteins
  theme_minimal() +
  labs(x = "NPX Diff BL v CTL", y = "NPX Diff FU v CTL", color = "Legend") +
  ggtitle("")


#Paired analysis 
#paired t-test between 0 and 6 month
#reduce data to only TN and FU samples significant TN samples
data_tn_fu <- data_meta[data_meta$`Disease status`=="TN" | data_meta$`Disease status`=="FU", ]
table(data_tn_fu$`Disease status`)
data_tn_fu_sig <- data_tn_fu[data_tn_fu$Assay %in% tn_hc_sig_names, ]

#get common_ID's of all patients with follow up samples ###Replace with SampleID for GEO users 
common_id_withFU <- unique(data_meta[data_meta$`Disease status`=="FU", ]$common_ID)

#reduce data to only these common IDs
data_tn_fu_sig <- data_tn_fu_sig[data_tn_fu_sig$common_ID %in% common_id_withFU, ]
View(data_tn_fu_sig)

paired_t <- olink_ttest(data_tn_fu_sig, variable = "Disease status", pair_id = "common_ID")
View(paired_t)
write.csv(paired_t, file = "paired_t_test.csv")
#FU is first level
paired_t <- read.csv(file = "paired_t_test.csv")
View(paired_t)

#subset increased proteins 
increased_paired <- paired_t[paired_t$Threshold=="Significant" & paired_t$estimate>0, ]
View(increased_paired)
length(increased_paired$Assay)
length(unique(increased_paired$Assay))


#subset decreased proteins
decreased_paired <- paired_t[paired_t$Threshold=="Significant" & paired_t$estimate<=0, ]
View(decreased_paired)
length(decreased_paired$Assay)
length(unique(decreased_paired$Assay))
sig_decr <- decreased_paired$Assay
overlap_sig_decr_BL <- intersect(sig_decr, bl_proteins_up)
overlap_sig_decr_BL_down <- intersect(sig_decr, bl_proteins_down)

#subset unchanged
unchanged <- paired_t[paired_t$Threshold=="Non-significant", ]
View(unchanged)
length(unique(unchanged$Assay))
length(unchanged$Assay)