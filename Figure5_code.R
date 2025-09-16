#Figure 5 code

#read in da-associated proteins
da_proteins <- read.csv(file = "") #supplementary Table 4

#read in skeletal muscle HPA enriched and enhanced specific list
muscle_hpa_all <- read.csv(file = "") #supplementary Table 5

muscle_overlap_all <- intersect(da_proteins$Assay, muscle_hpa_all$Gene)
print(muscle_overlap_all)

#read in skin tissue HPA enriched enhanced specific list
skin_all <- read.csv(file = "") #supplementary Table 5

skin_overlap_all <- intersect(da_proteins$Assay, skin_all$Gene)
print(skin_overlap_all)

#read in endothelial cell list 
endo_all <- read.csv(file = "") #supplementary Table 5

endo_all_overlap <- intersect(da_proteins$Assay, endo_all$Gene) 

#now create a dataframe of tissue specificity for disease-activity associated proteins
muscle_da_proteins <- data.frame(Protein = muscle_overlap_all, Tissue = "muscle")
skin_da_proteins <- data.frame(Protein = skin_overlap_all, Tissue = "skin")
endothelial_da_proteins <- data.frame(Protein = endo_all_overlap, Tissue = "endothelial")

tissue_sp_da_proteins <- rbind(muscle_da_proteins, skin_da_proteins, endothelial_da_proteins)
tissue_da_protein_list <- tissue_sp_da_proteins$Protein
View(tissue_sp_da_proteins)
write.csv(tissue_sp_da_proteins, file = "")


#subset lmer results to get estimates and p values then merge with tissue specific list to get tissue
da_proteins_subset <- da_proteins[da_proteins$Assay %in% tissue_da_protein_list, ]
da_protein_subset <- merge(da_proteins_subset, tissue_sp_da_proteins, by.x = "Assay", by.y = "Protein")
View(da_protein_subset)
write.csv(da_protein_subset, file = "tissue_da_proteins_wlmer.csv")

#there seems to be a significant enrichment of endothelial proteins; perform a proportion test to see if these are significantly different
tissue_protein_numbers <- c(46, 30, 20)
tissue_hpa_numbers <- c(290, 652, 503)
prop.test(tissue_protein_numbers, tissue_hpa_numbers)


###create a venn diagram ov the overlaps between skin, muscle, endo, da, all proteins
library(eulerr)

set_sizes <- c(
  skin = 503,
  muscle = 652,
  endo = 290,
  da = 1062,
  all = 2879,
  "skin&da"=20,
  "muscle&da"=30,
  "endo&da" = 46,
  "skin&all" = 61,
  "muscle&all" = 82,
  "endo&all" = 101,
  "da&all" = 1062
)

#generate Euler diagram
fit <- euler(set_sizes)
plot(fit) ##does not reflect overlaps
da_proteins <- da_proteins_all$Assay
all_proteins <- data_merged$Assay

devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
list_data <- list(
  HPA_Skin = skin_hpa_proteins,
  HPA_Muscle = muscle_hpa_proteins,
  HPA_Endothelial = endo_hpa_proteins,
  DA_proteins = da_proteins,
  Olink_proteins = all_proteins
)
list_data_unique <- lapply(list_data, unique)
ggvenn(list_data_unique)

list_data_skin <- list_data_unique[c(1, 4, 5)]
ggvenn(list_data_skin)

plot(euler(list_data_skin),
     quantities = TRUE,
     fills = c("#FFFF00", "#66CC99", "#6699CC"),
     edges = list(col = "black", lwd = 2),
     labels = list(font =2, cex = 1.2),
     main = "Skin Overlap")


list_data_muscle <- list_data_unique[c(2,4,5)]
plot(euler(list_data_muscle),
     quantities = TRUE,
     fills = c("#FF6666", "#66CC99", "#6699CC"),
     edges = list(col = "black", lwd = 2),
     labels = list(font =2, cex = 1.2),
     main = "Muscle Overlap")

list_data_endo <- list_data_unique[c(3,4,5)]
plot(euler(list_data_endo),
     quantities = TRUE,
     fills = c("#FFA500", "#66CC99", "#6699CC"),
     edges = list(col = "black", lwd = 2),
     labels = list(font =2, cex = 1.2),
     main = "Endothelial Overlap")


#Panel D, associating disease activity proteins with CHAQ, Physician Muscle VAS and Physician Skin VAS
skin <- read.csv(file = "skin_enriched_enhanced.csv") #Supplementary Table 5
skin_hpa <- skin$Gene
muscle_h <- read.csv(file = "muscle_enriched_enhanced.csv") #Supplementary Table 5
muscle_hpa <- muscle_h$Gene
endo <- read.csv(file = "endo_enriched_enchanced.csv") #Supplementary Table 5
endo_hpa <- endo$Gene
all_hpa <- c(skin_hpa, muscle_hpa, endo_hpa)
all_hpa_df <- data.frame(Proteins = all_hpa)
all_hpa_df$Tissue <- ifelse(all_hpa_df$Proteins %in% skin_hpa, "skin", ifelse(
  all_hpa_df$Proteins %in% muscle_hpa, "muscle", "endothelial"
))

View(all_hpa_df)
dim(all_hpa_df)

#read in protein data
data <- readRDS("") #replace with file name

#subset data to tissue specific proteins and remove CTL
data_hpa <- data[data$Assay %in% all_hpa & data$`Disease status`!="CTL", ]


#need to convert 97 to NA
data_hpa$phygskin_baseline_arm_1 <- ifelse(data_hpa$phygskin_baseline_arm_1>10, NA, data_hpa$phygskin_baseline_arm_1)
data_hpa$phygskin_6__2_month_fu_arm_1 <- ifelse(data_hpa$phygskin_6__2_month_fu_arm_1>10, NA, data_hpa$phygskin_6__2_month_fu_arm_1)

#now create a single time specific variable
data_hpa$time_sp_physkin <- ifelse(data_hpa$`Disease status`=="TN", data_hpa$phygskin_baseline_arm_1, data_hpa$phygskin_6__2_month_fu_arm_1)

physkin_lmer_hpa <- olink_lmer(df = data_hpa,
                               variable = "time_sp_physkin",
                               random = c("Site.x", "common_ID"), verbose = TRUE) 
View(physkin_lmer_hpa)
physkin_lmer_hpa_sig <- physkin_lmer_hpa[physkin_lmer_hpa$Threshold=="Significant", ]

physkin_posthoc_hpa <- olink_lmer_posthoc(df = data_hpa,
                                          olinkid_list = physkin_lmer_hpa_sig$OlinkID,
                                          variable = "time_sp_physkin",
                                          random = c("Site.x", "common_ID"),
                                          effect = "time_sp_physkin")
View(physkin_posthoc_hpa)
physkin_hpa_proteins <- merge(physkin_posthoc_hpa, all_hpa_df, by.x = "Assay", by.y = "Proteins")
View(physkin_hpa_proteins)
write.csv(physkin_hpa_proteins, file = "physkin_hpa_proteins.csv")


#order by tissue
physkin_top_20 <- physkin_hpa_proteins %>%
  group_by(Assay) %>%  # Group by Assay to handle duplicates
  slice_min(order_by = estimate, n = 1, with_ties = FALSE) %>%  # Take top 1 per Assay (most negative estimate)
  ungroup() %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order)) %>%  # Set tissue order
  arrange(estimate) %>%  # Sort by the most negative estimates
  mutate(Assay = factor(Assay, levels = rev(unique(Assay))))

tissue_colors <- c(
  "muscle" = "#F8766D",    # Medium Red
  "skin" = "#FFD700",      # Medium Yellow (Gold)
  "endothelial" = "#FFA500" # Medium Orange
)
physkin_plot <- ggplot(physkin_top_20, aes(x = Assay, y = estimate, color = Tissue)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.8) +
  scale_color_manual(values = tissue_colors) +
  coord_flip() +  # Flip coordinates for easier readability
  labs(
    title = "", #Proteins associated with Physician Skin VAS
    x = "Protein",
    y = "Estimate",
    color = "Tissue"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )
physkin_plot


##Muscle VAS
data_hpa$phygbmus_baseline_arm_1 <- ifelse(data_hpa$phygbmus_baseline_arm_1>10, NA, data_hpa$phygbmus_baseline_arm_1)
data_hpa$phygbmus_6__2_month_fu_arm_1 <- ifelse(data_hpa$phygbmus_6__2_month_fu_arm_1>10, NA, data_hpa$phygbmus_6__2_month_fu_arm_1)

#now create a single time specific variable
data_hpa$time_sp_phymusc <- ifelse(data_hpa$`Disease status`=="TN", data_hpa$phygbmus_baseline_arm_1, data_hpa$phygbmus_6__2_month_fu_arm_1)

phymusc_lmer_hpa <- olink_lmer(df = data_hpa,
                               variable = "time_sp_phymusc",
                               random = c("Site.x", "common_ID"), verbose = TRUE)  
View(phymusc_lmer_hpa)
phymusc_lmer_hpa_sig <- phymusc_lmer_hpa[phymusc_lmer_hpa$Threshold=="Significant", ]

phymusc_posthoc_hpa <- olink_lmer_posthoc(df = data_hpa,
                                          olinkid_list = phymusc_lmer_hpa_sig$OlinkID,
                                          variable = "time_sp_phymusc",
                                          random = c("Site.x", "common_ID"),
                                          effect = "time_sp_phymusc")
View(phymusc_posthoc_hpa)
phymusc_proteins_hpa <- merge(phymusc_posthoc_hpa, all_hpa_df, by.x = "Assay", by.y = "Proteins")
View(phymusc_proteins_hpa)
phymuscle_hpa_sig = phymusc_posthoc_hpa[phymusc_posthoc_hpa$Threshold=="Significant", ]
dim(phymusc_lmer_hpa_sig)
length(unique(phymusc_lmer_hpa_sig$Assay))
write.csv(phymusc_proteins_hpa, file = "physmuscle_hpa_sig.csv") #46 endothelial, 31 muscle, 16 skin

#take top 50 for visualization
phymuscle_top_50 <- phymusc_proteins_hpa %>%
  group_by(Assay) %>%  # Group by Assay to handle duplicates
  slice_min(order_by = estimate, n = 1, with_ties = FALSE) %>%  # Take top 1 per Assay (most negative estimate)
  ungroup() %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order)) %>%  # Set tissue order
  arrange(estimate) %>%  # Sort by the most negative estimates
  mutate(Assay = factor(Assay, levels = rev(unique(Assay)))) %>%
  slice(1:50)

tissue_colors <- c(
  "muscle" = "#F8766D",    # Medium Red
  "skin" = "#FFD700",      # Medium Yellow (Gold)
  "endothelial" = "#FFA500" # Medium Orange
)
phymuscle_plot <- ggplot(phymuscle_top_50, aes(x = Assay, y = estimate, color = Tissue)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.8, fatten = 0.6) +
  scale_color_manual(values = tissue_colors) +
  coord_flip() +  # Flip coordinates for easier readability
  labs(
    title = "", #Top 20 Proteins associated with Physician Muscle VAS
    x = "Protein",
    y = "Estimate",
    color = "Tissue"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )
phymuscle_plot



#CHAQ
data_hpa$chaqindx_baseline_arm_1 <- as.numeric(data_hpa$chaqindx_baseline_arm_1)
summary(data_hpa$chaqindx_baseline_arm_1)
data_hpa$chaqindx_6__2_month_fu_arm_1 <- as.numeric(data_hpa$chaqindx_6__2_month_fu_arm_1)
summary(data_hpa$chaqindx_6__2_month_fu_arm_1)

#now create a single time specific variable
data_hpa$time_sp_chaq <- ifelse(data_hpa$`Disease status`=="TN", data_hpa$chaqindx_baseline_arm_1, data_hpa$chaqindx_6__2_month_fu_arm_1)

chaq_lmer_hpa <- olink_lmer(df = data_hpa,
                            variable = "time_sp_chaq",
                            random = c("Site.x", "common_ID"), verbose = TRUE)  
View(chaq_lmer_hpa) #11 missing data
chaqlmer_sig <- chaq_lmer_hpa[chaq_lmer_hpa$Threshold=="Significant", ]

chaq_posthoc_hpa <- olink_lmer_posthoc(df = data_hpa,
                                       olinkid_list = chaqlmer_sig$OlinkID,
                                       variable = "time_sp_chaq",
                                       random = c("Site.x", "common_ID"),
                                       effect = "time_sp_chaq")
View(chaq_posthoc_hpa)

chaq_proteins_hpa <- merge(chaq_posthoc_hpa, all_hpa_df, by.x = "Assay", by.y = "Proteins")
View(chaq_proteins_hpa)
write.csv(chaq_proteins_hpa, file = "chaq_proteins_hpa.csv")

table(chaq_proteins_hpa$Threshold) #71
table(chaq_proteins_hpa$Threshold, chaq_proteins_hpa$Tissue) #33 muscle, 25 endothelial, 13 skin
chaq_proteins_hpa_sig <- chaq_proteins_hpa[chaq_proteins_hpa$Threshold=="Significant", ]
dim(chaq_proteins_hpa_sig)
length(unique(chaq_proteins_hpa_sig$Assay))

#take top 50 for visualization
chaq_top50 <- chaq_proteins_hpa %>%
  group_by(Assay) %>%  # Group by Assay to handle duplicates
  slice_min(order_by = estimate, n = 1, with_ties = FALSE) %>%  # Take top 1 per Assay (most negative estimate)
  ungroup() %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order)) %>%  # Set tissue order
  arrange(estimate) %>%  # Sort by the most negative estimates
  mutate(Assay = factor(Assay, levels = rev(unique(Assay)))) %>%
  slice(1:50)

tissue_colors <- c(
  "muscle" = "#F8766D",    # Medium Red
  "skin" = "#FFD700",      # Medium Yellow (Gold)
  "endothelial" = "#FFA500" # Medium Orange
)
chaq_plot <- ggplot(chaq_top50, aes(x = Assay, y = estimate, color = Tissue)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.8, fatten = 0.6) +
  scale_color_manual(values = tissue_colors) + 
  coord_flip() +  # Flip coordinates for easier readability
  labs(
    title = "", #Top 20 Proteins associated with CHAQ
    x = "Protein",
    y = "Estimate",
    color = "Tissue"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )
chaq_plot



