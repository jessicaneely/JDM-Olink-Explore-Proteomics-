#Figure 6 code 
library(OlinkAnalyze)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)


#read in protein data
data <- readRDS("merged_data_protein_meta_5_24_24.rds")

#subset to only TN
table(data$`Disease status`)
data_TN <- data_msa[data_msa$`Disease status`=="TN", ]

#remove duplicate sample: CLG13
data_TN <- data_TN[data_TN$common_ID!="CLG13", ] 

#run lmer with NXP2, TIF1y and MDA5
nxp2_lmer <- olink_lmer(df = data_TN,
                        variable = "antimmj_omrf",
                        random = "Site.x")
#take "Assay" and "Threshold columns"
nxp2_lmer_cols <- nxp2_lmer[ , c("Assay", "Adjusted_pval")]
View(nxp2_lmer_cols)

#renanxp2_lmer_cols#rename "threshold" column to facilitate merging
names(nxp2_lmer_cols)[2] <- "Adjusted_pval_lmer"

p155_lmer <- olink_lmer(df = data_TN,
                        variable = "antip155_omrf",
                        random = "Site.x") #no significant proteins

mda5_lmer <- olink_lmer(df = data_TN,
                        variable = "antimda5_omrf",
                        random = "Site.x")
#take "Assay" and "Threshold columns"
mda5_lmer_cols <- mda5_lmer[ , c("Assay", "Adjusted_pval")]

#rename "threshold" column to facilitate merging
names(mda5_lmer_cols)[2] <- "Adjusted_pval_lmer"

#run post-hoc on nxp2 and mda5 lmer
nxp2_lmer_posthoc <- olink_lmer_posthoc(df = data_TN,
                                        olinkid_list = data_TN$OlinkID,
                                        variable = "antimmj_omrf",
                                        random = c("Site.x"),
                                        effect = "antimmj_omrf")
#add results together to identify those proteins significant in both 
nxp2_lmer_posthoc <- right_join(nxp2_lmer_posthoc, nxp2_lmer_cols, by = "Assay")

#Create a column that indicates significance
nxp2_lmer_posthoc$sig_proteins <- ifelse(nxp2_lmer_posthoc$Threshold=="Significant" & nxp2_lmer_posthoc$Adjusted_pval_lmer<0.1, "Significant", "Not_significant")
write.csv(nxp2_lmer_posthoc, file = "")

blood_vessel <- c("CLEC14A","ACVRL1","PGF", "LAMA4", "COL15A1","COL4A1")
nxp2_lmer_posthoc$enrichment <- ifelse(nxp2_lmer_posthoc$Assay %in% blood_vessel & nxp2_lmer_posthoc$sig_proteins=="Significant", "GOBP: blood vessel development", ifelse(
  nxp2_lmer_posthoc$sig_proteins=="Significant", "Significant", "Not_Significant"
)
)

# Define custom color palette
custom_colors <- c("Significant" = "blue", "Not_Significant" = "gray", "GOBP: blood vessel development" = "red")

ggplot(nxp2_lmer_posthoc, aes(x = estimate, y = -log10(Adjusted_pval), color = enrichment)) +
  geom_point() +
  
  # Highlight significant proteins with text labels
  geom_text_repel(data = nxp2_lmer_posthoc[nxp2_lmer_posthoc$sig_proteins == 'Significant', ],
                  aes(label = Assay),
                  box.padding = 0.25,
                  point.padding = 0.5,
                  min.segment.length = 0.5,
                  size = 4) +
  
  # Define custom colors for the combined categories
  scale_color_manual(values = custom_colors) +
  
  # Additional plot formatting
  theme_minimal() +
  labs(title = "NXP2-associated proteins in BL JDM samples",
       x = "Estimate (NPX Difference)",
       y = "-log10(Adjusted p-value)",
       color = "Adjusted p-val lmer<0.1") +
  
  # Adjust the position of the legend and title
  theme(legend.position = "right",
        plot.title = element_text(size = 11, hjust = 0.5))


#mda5 posthoc
mda5_lmer_posthoc <- olink_lmer_posthoc(df = data_TN,
                                        olinkid_list = data_TN$OlinkID,
                                        variable = "antimda5_omrf",
                                        random = c("Site.x"),
                                        effect = "antimda5_omrf")
View(mda5_lmer_posthoc)

#add original lmer results
mda5_lmer_posthoc <- right_join(mda5_lmer_posthoc, mda5_lmer_cols, by = "Assay")

#Create a column that indicates significance
mda5_lmer_posthoc$sig_proteins <- ifelse(mda5_lmer_posthoc$Threshold=="Significant" & mda5_lmer_posthoc$Adjusted_pval_lmer<0.1, "Significant", "Not_significant")

write.csv(mda5_lmer_posthoc, file = "")

surfactant_proteins <- c("SFTPA1", "SFTPA2")
typeIII_ifn_proteins <- c("IFNL1")
mda5_lmer_posthoc$enrichment <- ifelse(mda5_lmer_posthoc$Assay %in% surfactant_proteins & mda5_lmer_posthoc$sig_proteins=="Significant", "Surfactant Proteins", ifelse(
  mda5_lmer_posthoc$Assay %in% typeIII_ifn_proteins & mda5_lmer_posthoc$sig_proteins=="Significant", "Type III IFN", ifelse(
    mda5_lmer_posthoc$sig_proteins=="Significant", "Significant", "Not Significant"
  )
))

custom_colors_mda5 <- c("Not Significant" = "grey",
                        "Significant" = "blue",
                        "Surfactant Proteins" = "#FFA500",
                        "Type III IFN" = "red")

ggplot(mda5_lmer_posthoc, aes(x = estimate, y = -log10(Adjusted_pval), color = enrichment)) +
  geom_point() + # Now uses the 'Significant' column for coloring
  geom_text_repel(data = mda5_lmer_posthoc[mda5_lmer_posthoc$sig_proteins == 'Significant', ], 
                  aes(label = Assay),  # Ensure this column exists
                  box.padding = 0.25, 
                  point.padding = 0.5, 
                  min.segment.length = 0.5, 
                  size = 3) +
  scale_color_manual(values = custom_colors_mda5) +
  theme_minimal() +
  labs(title = "MDA5-associated proteins in BL JDM samples",
       x = "Estimate (NPX Difference)",
       y = "-log10(Adjusted p-value)",
       color = "Adjusted p-val lmer<0.1") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))




#Point range plots 
#function edited from OlinkAnalyze package for generating point range plots:
# function for olink lmer plot edited to adjust x-axis order and plot titles
olink_lmer_plot_edited_msa <- function (df, variable, outcome = "NPX", random, olinkid_list = NULL,
                                        covariates = NULL, x_axis_variable, col_variable = NULL,
                                        number_of_proteins_per_plot = 8, verbose = FALSE, ..., x_axis_order = NULL, plot_titles = NULL) {
  if (missing(df) | missing(variable) | missing(x_axis_variable) | missing(random)) {
    stop("The df, variable, random and x_axis_variable arguments need to be specified.")
  }
  if (!all(x_axis_variable %in% unique(unlist(strsplit(variable, "[\\*:]"))))) {
    stop("The x axis variable must be included in the variable argument.")
  }
  if (!is.null(col_variable)) {
    if (!all(col_variable %in% unique(unlist(strsplit(variable, "[\\*:]"))))) {
      stop("The color variable must be included in the variable argument.")
    }
  }
  if (length(list(...)) > 0) {
    ellipsis_variables <- names(list(...))
    if (length(ellipsis_variables) == 1) {
      if (!(ellipsis_variables == "coloroption")) {
        stop(paste0("The ... option only takes the coloroption argument. ... currently contains the variable ", ellipsis_variables, "."))
      }
    } else {
      stop(paste0("The ... option only takes one argument. ... currently contains the variables ", paste(ellipsis_variables, collapse = ", "), "."))
    }
  }
  df <- df %>% dplyr::filter(stringr::str_detect(OlinkID, "OID[0-9]{5}"))
  if (is.null(olinkid_list)) {
    olinkid_list <- df %>% dplyr::select(OlinkID) %>% dplyr::distinct() %>% dplyr::pull()
  }
  if (is.null(col_variable)) {
    current_fixed_effect <- x_axis_variable
    color_for_plot <- x_axis_variable
  } else {
    current_fixed_effect <- paste0(x_axis_variable, ":", col_variable)
    color_for_plot <- col_variable
  }
  lm.means <- olink_lmer_posthoc(df = df, variable = variable, random = random, outcome = outcome, olinkid_list = olinkid_list,
                                 covariates = covariates, effect = current_fixed_effect, mean_return = TRUE, verbose = verbose) %>%
    dplyr::mutate(Name_Assay = paste0(Assay, "_", OlinkID))
  assay_name_list <- lm.means %>%
    dplyr::mutate(OlinkID = factor(OlinkID, levels = olinkid_list)) %>%
    dplyr::arrange(OlinkID) %>%
    dplyr::pull(Name_Assay) %>%
    unique()
  lm.means <- lm.means %>%
    dplyr::mutate(Name_Assay = factor(Name_Assay, levels = assay_name_list))
  if (!is.null(x_axis_order)) {
    lm.means <- lm.means %>%
      dplyr::mutate(!!rlang::ensym(x_axis_variable) := factor(!!rlang::ensym(x_axis_variable), levels = x_axis_order))
  }
  topX <- length(assay_name_list)
  protein_index <- seq(from = 1, to = topX, by = number_of_proteins_per_plot)
  list_of_plots <- list()
  COUNTER <- 1
  for (i in seq_along(protein_index)) {
    from_protein <- protein_index[i]
    to_protein <- if ((protein_index[i] + number_of_proteins_per_plot) > topX) topX + 1 else protein_index[i + 1]
    assays_for_plotting <- assay_name_list[c(from_protein:(to_protein - 1))]
    
    # Map plot_titles to the corresponding facets
    lm.means <- lm.means %>%
      dplyr::mutate(Facet_Title = factor(Name_Assay, levels = assays_for_plotting,
                                         labels = plot_titles[from_protein:(to_protein - 1)]))
    
    lmerplot <- lm.means %>%
      dplyr::filter(Name_Assay %in% assays_for_plotting) %>%
      ggplot2::ggplot() +
      ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
      ggplot2::ylab("NPX") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10)) +
      ggplot2::geom_pointrange(ggplot2::aes(x = as.factor(!!rlang::ensym(x_axis_variable)),
                                            y = emmean, ymin = conf.low, ymax = conf.high,
                                            color = as.factor(!!rlang::ensym(color_for_plot))),
                               position = ggplot2::position_dodge(width = 0.4),
                               size = 0.8) +
      ggplot2::facet_wrap(~Facet_Title, scales = "free_y", nrow = 2, ncol = 4) +  # Use Facet_Title for individual titles
      OlinkAnalyze::olink_color_discrete(...) +
      OlinkAnalyze::set_plot_theme() +
      ggplot2::labs(x = "anti-MDA5 status", color = "anti-MDA5 status") + 
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),  # Remove box around facet labels
        strip.text = ggplot2::element_text(hjust = 0.5, size = 14) # Center and style facet titles
      )
    
    list_of_plots[[COUNTER]] <- lmerplot
    COUNTER <- COUNTER + 1
  }
  return(invisible(list_of_plots))
}


ifnl1_mda5_plot <- olink_lmer_plot_edited_msa(df=data_TN,
                                              variable = "antimda5_omrf",
                                              outcome = "NPX",
                                              random = "Site.x",
                                              x_axis_variable = "antimda5_omrf",
                                              olinkid_list = c("OID20795", "OID30129"),
                                              plot_titles = c("IFNL1", "IFNW1"))
ifnl1_mda5_plot

surfactant_plot <- olink_lmer_plot_edited_msa(df=data_TN,
                                              variable = "antimda5_omrf",
                                              outcome = "NPX",
                                              random = "Site.x",
                                              x_axis_variable = "antimda5_omrf",
                                              olinkid_list = c("OID21165","OID21218"),
                                              plot_titles = c("SFTPA1", "SFTPA2"))
surfactant_plot
