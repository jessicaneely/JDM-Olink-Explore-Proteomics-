#Figure 2
library(OlinkAnalyze)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggrepel)

#set working directory
data <- readRDS("")

#subset data to TN and CTL
data_TN_CTL <- data[(data$`Disease status`=="TN" | data$`Disease status`=="CTL"), ]

# Run Linear mixed model using OlinkAnalyze package
lmer.TN_CTL.site <- olink_lmer(df = data_TN_CTL,
                                variable = 'Disease_status',
                                random = 'Site.x')


lmer.TN_CTL.site_posthoc <- olink_lmer_posthoc(df = data_TN_CTL,
                                        olinkid_list = data_TN_CTL$OlinkID,
                                        variable = 'Disease_status',
                                        random = c("Site.x"),
                                        effect = 'Disease_status')

# specify if Assay significant in post-hoc or initial LMER and add data to post-hoc results for a single dataframe 

#save results
write.csv(lmer.TN_CTL.site_posthoc, file = "lmer_tn_hc_posthoc.csv")

#load in list of differentially expressed proteins between TN and HC 
tn_hc_results <- read.csv(file = "lmer_tn_hc_posthoc.csv")

#rename the estimated difference between contrast terms 
tn_hc_results  <- tn_hc_results  %>%
  rename(NPX_difference = estimate..difference.in.NPX.value.between.contrast.terms.)
View(tn_hc_results)

# deal with duplicated proteins
tn_hc_filt <- tn_hc_results %>%
  group_by(Assay) %>%
  slice(which.max(abs(NPX_difference))) %>%
  ungroup()
View(tn_hc_filt)

#take only proteins with both lmer and post hoc results significant
tn_hc_sig <- tn_hc_filt[tn_hc_filt$Threshold=="Significant" & tn_hc_filt$significant_lmer_no_posthoc=="Significant", ]
View(tn_hc_sig)


#set proteins to label
top_proteins <- tn_hc_sig[order(abs(tn_hc_sig$NPX_difference), decreasing = TRUE), ][1:50, ]
head(top_proteins)


ggplot(tn_hc_results, aes(x = NPX_difference, y = -log10(Adjusted_pval))) +
  geom_point(aes(color =  (Threshold=="Significant" & significant_lmer_no_posthoc=="Significant"))) +
  geom_text_repel(data = top_proteins, aes(label = Assay), 
                  box.padding = 0.25, point.padding = 0.5, 
                  min.segment.length = 0.5, size = 3)  +
  scale_color_manual(values = c("grey", "blue")) +
  #coord_cartesian(xlim = c(-6, 6)) + 
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "NPX Difference",
       y = "-log10(adj pvalue)",
       color = "Significant") +
  theme(legend.position = "top")

#Repeat analysis for FU v CTL
data_FU_CTL <- data[(data$`Disease status`=="FU" | data$`Disease status`=="CTL"), ]

# Run Linear mixed model using OlinkAnalyze package
lmer.FU_CTL.site <- olink_lmer(df = data_FU_CTL,
                               variable = 'Disease_status',
                               random = 'Site.x')


lmer.FU_CTL.site_posthoc <- olink_lmer_posthoc(df = data_FU_CTL,
                                               olinkid_list = data_FU_CTL$OlinkID,
                                               variable = 'Disease_status',
                                               random = c("Site.x"),
                                               effect = 'Disease_status')

# specify if Assay significant in post-hoc or initial LMER and add data to post-hoc results for a single dataframe 

#take those significant in both lmer and post hoc and save results
write.csv(lmer.FU_CTL.site_posthoc , file = "FU_hc_posthoc.csv")

#now make volcano plot of FU v CTL results
fu_hc_results <- read.csv(file = "FU_ctl_posthoc.csv"))

fu_hc_results  <- fu_hc_results  %>%
  rename(NPX_difference = estimate..difference.in.NPX.value.between.contrast.terms.)

#deal with duplicated proteins
fu_hc_filt <- fu_hc_results %>%
  group_by(Assay) %>%
  slice(which.max(abs(NPX_difference))) %>%
  ungroup()

#take only proteins with both lmer and post hoc results significant
fu_hc_sig <- fu_hc_filt[fu_hc_filt$Threshold=="Significant" & fu_hc_filt$significant_lmer_no_posthoc=="Significant", ]
dim(fu_hc_sig)
View(fu_hc_sig)

top_proteins_fu <- fu_hc_sig[order(abs(fu_hc_sig$NPX_difference), decreasing = TRUE), ][1:50, ]
head(top_proteins_fu)

ggplot(fu_hc_results, aes(x = NPX_difference, y = -log10(Adjusted_pval))) +
  geom_point(aes(color =  (Threshold=="Significant" & significant_lmer_no_posthoc=="Significant"))) +
  geom_text_repel(data = top_proteins_fu, aes(label = Assay), 
                  box.padding = 0.25, point.padding = 0.5, 
                  min.segment.length = 0.5, size = 3)  +
  scale_color_manual(values = c("grey", "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "NPX Difference",
       y = "-log10(adj pvalue)",
       color = "Significant") +
  theme(legend.position = "top")

#generate point-range plots for top differentially-expressed proteins in each
#first need to generate function edited from OlinkAnalyze package:
# function for olink lmer plot edited to adjust x-axis order and plot titles
olink_lmer_plot_edited <- function (df, variable, outcome = "NPX", random, olinkid_list = NULL,
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
      ggplot2::labs(x = "Disease status", color = "Disease status") + 
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),  # Remove box around facet labels
        strip.text = ggplot2::element_text(hjust = 0.5, size = 14) # Center and style facet titles
      )
    
    list_of_plots[[COUNTER]] <- lmerplot
    COUNTER <- COUNTER + 1
  }
  return(invisible(list_of_plots))
}

#read in protein data
data <- readRDS("merged_data_protein_meta_5_24_24.rds") #replace with name of Olink dataset 

table(data$`Disease status`)
data <- data %>%
  rename(disease_status = "Disease status")

data <- data %>% dplyr::mutate(disease_status = dplyr::recode(disease_status, "TN" = "BL"))
table(data$disease_status)

data_bl_ctl <- data[data$disease_status!="FU", ]
data_fu_ctl <- data[data$disease_status!="BL", ]


#the dots represents the fixed effect estimate and the line the 95% confidence interval summarizing the NPX difference 
#between groups accounting for the random effect


# top 8 proteins with greatest NPX change in BL v CTL
bl_ctl_top8 <- olink_lmer_plot_edited(data_bl_ctl,
                                      variable = "disease_status",
                                      random = "Site.x",
                                      x_axis_variable = "disease_status",
                                      olinkid_list = c("OID20697", "OID20050", "OID31128", "OID31290", "OID31254", "OID20523", "OID30131", "OID30251"),
                                      x_axis_order = c("BL", "CTL"),
                                      plot_titles = c("CXCL10", "TNNI3", "GBP1", "TLR2", "RAD51", "CCL7", "MYOM2", "MYL3"))
bl_ctl_top8

# top 8 proteins with greatest NPX change in FU v CTL
fu_ctl_top8 <- olink_lmer_plot_edited(data_fu_ctl,
                                      variable = "disease_status",
                                      random = "Site.x",
                                      x_axis_variable = "disease_status",
                                      olinkid_list = c("OID30836", "OID31290", "OID31454", "OID30862", "OID21277", "OID20794", "OID20427", "OID31353"),
                                      x_axis_order = c("FU", "CTL"),
                                      plot_titles = c("DOCK9", "TLR2", "BCL2", "IGDCC3", "STAT5B", "PTEN", "IL1B", "SMAD2"))
fu_ctl_top8                                            

