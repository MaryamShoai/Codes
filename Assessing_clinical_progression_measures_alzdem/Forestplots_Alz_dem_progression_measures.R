library(readr)
library(metafor)
library(dplyr)
library(ggforestplot)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(emmeans)
library(cowplot)
library(ggstance)
library(ggplot2)


allmeasures<- read_delim("~/MShoaee/Students_Consultations/MSc students/JonathanMc/allmeasures_reworked_ADstaged.txt", 
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Create a named vector with the divisors for each measure
divisors <- c("ADAS-Cog-11" = 70, "ADAS-Cog-13" = 85, "ADAS-Cog-14" = 90, 
              "ADCS-ADL-MCI" = 53, "ADCS-iADL" = 56, "CDR-SB" = 18, 
              "FAQ" = 30, "iADRS" = 146, "MMSE" = 30, "NPI" = 144, 
              "NTB" = 6, "QOL-AD" = 39, "ADCS-ADL" = 78)

# Search for each measure, standardise based on range of scale. 
allmeasures <- allmeasures %>% 
  arrange(measure) %>%
  mutate(Change = abs(Change),
         Change = case_when(
           measure %in% names(divisors) ~ Change / divisors[measure],
           TRUE ~ Change
         ),
         SE = case_when(
           measure %in% names(divisors) ~ SE / divisors[measure],
           TRUE ~ SE
         ))

#subset the measures individually.
for (meas in names(divisors)) {
  assign(meas, allmeasures %>% filter(measure == meas))
}


#Plot individual forest plots Suppl. Figure S4-S14
##### 

create_meta_analysis_plot <- function(dat, formula, measure_name, metagroup) {
  
colnames(dat) <- c("Study", "type_of_change", "Change", "SE", 
                     "NumberofPlaceboPt.s", "measure", "PercentageF", 
                     "Age", "AD_stage")
  
  dat$sd <- as.numeric(dat$SE) * sqrt(as.numeric(dat$NumberofPlaceboPt.s))
  dat$m2i <- 0 
  dat$sd2i <- 0
  dat$ri <- 0
  
  # Transform Leaste Square Mean to Mean Change
  dat <- escalc(measure = "MC", m1i = as.numeric(dat$Change), 
                m2i = as.numeric(dat$m2i), sd1i = as.numeric(dat$sd), 
                sd2i = as.numeric(dat$sd2i), ri = as.numeric(dat$ri), 
                ni = as.numeric(dat$NumberofPlaceboPt.s), data = dat)
  
  res_mod <- rma(yi, vi, mods = formula, data = dat, 
                 slab = Study, control = list(maxiter = 1000))
  res <- rma(yi, vi, data = dat, slab = Study)
  sav <- emmprep(res_mod)
  TMP <- emmeans(sav, specs = "1", type = "response", weights = "proportional")
  
  dat <- dat %>% 
    add_row(Study = "Meta-regressed", Change = summary(TMP)$emmean, 
            SE = summary(TMP)$SE, measure = measure_name, 
            metagroup = metagroup)
  
  # ... (Calculate z-score and p-values for meta-regressed results)
  z <- summary(TMP)$emmean / summary(TMP)$SE
  a <- -0.717 * z
  b <- 0.416 * z^2
  p_value <- exp(a - b) 
  
  fplot <- forestplot(
    df = dat, labeltext = c("Study", "AD_stage"), 
    estimate = Change,
    title = measure_name,
    colour = metagroup, 
    name = Study, se = SE, ci = 0.95, clip = c(-.01, 0.04), xlim = c(-.01, .04)
  ) +
    theme(axis.title.y = element_text(family = 'Arial', size = 10)) + 
    theme(axis.title.x = element_text(family = 'Arial', size = 8)) + 
    theme(legend.text = element_text(family = 'Arial', size = 8)) + 
    theme(axis.text = element_text(family = 'Arial', size = 8), 
          axis.title.x = element_text(colour = "white")) +
    theme(plot.title = element_text(family = 'Arial', size = 12)) +
    theme(legend.title = element_text(family = 'Arial', size = 8)) +
    scale_colour_manual(values = c(
      Functional = "blue", Cognitive = "red", 
      Composite = "green", Neuropsychiatric = "purple"
    )) +
    scale_fill_manual(values = c(
      Functional = "lightblue", Cognitive = "pink", 
      Composite = "lightgreen", Neuropsychiatric = "lavender"
    )) +
    ggforestplot::geom_effect(
      aes(x = Change, xmin = .data$.xmin, xmax = .data$.xmax, 
          colour = metagroup, fill = metagroup),
      position = ggstance::position_dodgev(height = 0.5),
      fatten = 1
    )
  
  return(list(plot = fplot, p_value = p_value)) 
}


meta_regressedres <- as.data.frame(matrix(ncol = 11))
names(meta_regressedres) <- c("Test", "Beta_original", "SE_original", "L95_original", "U95_original", "Pval_original", "Beta_mod", "SE_mod", "L95_mod", "U95_mod", "Pval_mod")

measureslist <- c("ADCS-ADL-MCI", "ADCS-iADL", "ADCS-ADL", "FAQ", "QOL-AD", "ADAS-Cog-11", "ADAS-Cog-13", "ADAS-Cog-14", "MMSE", "CDR-SB", "iADRS", "NPI")
formulaslist <- list(
  "ADCS-ADL-MCI" = "~ Age + PercentageF",
  "ADCS-iADL" = "~ Age + PercentageF + AD_stage",
  "ADCS-ADL" = "~ Age + PercentageF + AD_stage",
  "FAQ" = "~ Age + PercentageF + AD_stage",
  "QOL-AD" = "~ AD_stage",
  "ADAS-Cog-11" = "~ Age + PercentageF + AD_stage",
  "ADAS-Cog-13" = "~ Age + PercentageF + AD_stage",
  "ADAS-Cog-14" = "~ AD_stage",
  "MMSE" = "~ Age + PercentageF + AD_stage",
  "CDR-SB" = "~ Age + PercentageF + AD_stage",
  "iADRS" = "~ Age + PercentageF + AD_stage",
  "NPI" = "~ Age + PercentageF + AD_stage"
)

# Define metagroups
metagroups <- list(
  Functional = c("ADCS-ADL-MCI", "ADCS-iADL", "ADCS-ADL", "FAQ", "QOL-AD"),
  Cognitive = c("ADAS-Cog-11", "ADAS-Cog-13", "ADAS-Cog-14", "MMSE"),
  Composite = c("CDR-SB", "iADRS"),
  Neuropsychiatric = c("NPI")
)

# Create a reverse mapping from measure to metagroup
measure_to_metagroup <- unlist(lapply(names(metagroups), function(group) {
  setNames(rep(group, length(metagroups[[group]])), metagroups[[group]])
}))




for (i in 1:length(measureslist)) {
  # Use the constructed name with backticks to get the dataframe
  dat <- get(measureslist[i]) 
  colnames(dat) <- c("Study", "type_of_change", "Change", "SE", "NumberofPlaceboPt.s", "measure", "PercentageF", "Age", "AD_stage")
  dat$sd <- as.numeric(dat$SE) * sqrt(as.numeric(dat$NumberofPlaceboPt.s))
  
  dat$m2i <- 1 * 0
  dat$sd2i <- 1 * 0
  dat$ri <- 1 * 0
  # Transform LSM to MC
  dat <- escalc(measure = "MC", m1i = as.numeric(dat$Change), m2i = as.numeric(dat$m2i), sd1i = as.numeric(dat$sd), sd2i = as.numeric(dat$sd2i), ri = as.numeric(dat$ri), ni = as.numeric(dat$NumberofPlaceboPt.s), data = dat)
  
  # Construct the formula for the current measure
  THEFORMULA <- as.formula(formulaslist[[measureslist[i]]])
  
  res_mod <- rma(yi, vi, mods = THEFORMULA, data = dat, slab = Study, control = list(maxiter = 1000))
  res <- rma(yi, vi, data = dat, slab = Study)
  sav <- emmprep(res_mod)
  TMP <- emmeans(sav, specs = "1", type = "response", weights = "proportional")
  
  dat[(nrow(dat) + 1), 1] <- "Meta-regressed"
  dat[nrow(dat), 3] <- summary(TMP)$emmean
  dat[nrow(dat), 4] <- summary(TMP)$SE
  dat[nrow(dat), 6] <- measureslist[i]
  meta_regressedres[i, 1] <- measureslist[i]
  meta_regressedres[i, 2] <- res$beta[1]
  meta_regressedres[i, 3] <- res$se
  meta_regressedres[i, 4] <- res$ci.lb[1]
  meta_regressedres[i, 5] <- res$ci.ub[1]
  meta_regressedres[i, 6] <- res$pval[1]
  meta_regressedres[i, 7] <- summary(TMP)$emmean
  meta_regressedres[i, 8] <- summary(TMP)$SE
  meta_regressedres[i, 9] <- summary(TMP)$asymp.LCL
  meta_regressedres[i, 10] <- summary(TMP)$asymp.UCL
  
  # Calculate z-score and p-values for metaregressed results
  z <- summary(TMP)$emmean / summary(TMP)$SE
  a <- -0.717 * z
  b <- 0.416 * z^2
  meta_regressedres[i, 11] <- exp(a - b)
  
  # Add metagroup to dat
  dat$metagroup <- measure_to_metagroup[[measureslist[i]]]
  
  fplot <- forestplot(
    df = dat, labeltext = c("Study", "AD_stage"), 
    estimate = Change,
    title = measureslist[i],
    colour = metagroup,  # Color by metagroup
    name = Study, se = SE, ci = 0.95, clip = c(-.01, 0.04), xlim = c(-.01, .04)
  ) +
    theme(axis.title.y = element_text(family = 'Arial', size = 10)) + 
    theme(axis.title.x = element_text(family = 'Arial', size = 8)) + 
    theme(legend.text = element_text(family = 'Arial', size = 8)) + 
    theme(axis.text = element_text(family = 'Arial', size = 8), axis.title.x = element_text(colour = "white")) +
    theme(plot.title = element_text(family = 'Arial', size = 12)) +
    theme(legend.title = element_text(family = 'Arial', size = 8)) +
    scale_colour_manual(values = c(Functional = "blue", Cognitive = "red", Composite = "green", Neuropsychiatric = "purple")) +
    scale_fill_manual(values = c(Functional = "lightblue", Cognitive = "pink", Composite = "lightgreen", Neuropsychiatric = "lavender")) +  # Add scale_fill_manual
    ggforestplot::geom_effect(aes(x = Change, xmin = .data$.xmin, xmax = .data$.xmax, colour = metagroup, fill = metagroup),  # Define colour and fill inside aes
                              position = ggstance::position_dodgev(height = 0.5),
                              fatten = 1)
  ggsave(filename = paste0("forest_", measureslist[i], ".tiff"), plot = fplot, 
         width = 20, height = 20, units = "cm", dpi = 300)
  assign(paste0("forest_", measureslist[i]), fplot)
}



#Plot based on meta group Figure 2. 
#####

library(readr)
library(metafor)
library(dplyr)
library(ggforestplot)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(emmeans)
library(cowplot)
library(ggstance) 



# Create a function to extract and bind data from forest_X lists where X is the name of the progression measure. group_name refers to the metagroup eg Functional etc...

bind_forest_data <- function(group_name) {
  object_names <- metagroups[[group_name]]
  
  data_list <- list()
    for (i in seq_along(object_names)) {
    forest_list_name <- paste0("forest_", object_names[i])
    # Extract the "data" element from the forest plot list
    data_list[[i]] <- get(forest_list_name)[["data"]]
  }
  do.call(rbind, data_list)
}

Composite_data <- bind_forest_data("Composite")
Functional_data <- bind_forest_data("Functional")
Cognitive_data <- bind_forest_data("Cognitive")
Neuropsychiatric_data <- bind_forest_data("Neuropsychiatric")


Functional_metaregress<- subset(Functional_data, Functional_data$Study=="Meta-regressed")
Functional_data<- Functional_data %>%  subset(.,Functional_data$Study!="Meta-regressed")


Functional_forest<- forestplot(
  df = Functional_data, labeltext = c(Study), 
  estimate = Change,
  title= "Functional",
  colour= measure,
  name = Study, se = SE, ci = 0.95,clip = c(-.01, 0.04),xlim=c(-.01,.04)) +
  theme(axis.title.y = element_text(family = 'Arial', size = 10))+ 
  theme(axis.title.x = element_text(family = 'Arial',size = 8))+ 
  theme(legend.text = element_text(family = 'Arial',size = 8)) + 
  theme(axis.text = element_text(family = 'Arial',size = 8), axis.title.x=element_text(colour="white")) +
  theme(plot.title = element_text(family = 'Arial', size = 12)) +
  theme(legend.title = element_text(family = 'Arial', size = 8)) +
  scale_colour_discrete(na.translate = F)+ 
  ggforestplot::geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax, 
                                         colour = .data$measure, filled = .data$.filled), 
                            position = ggstance::position_dodgev(height = 0.5),
                            fatten=1)

Functional_metaregress_forest<- forestplot(
  df = Functional_metaregress, 
  estimate = Change,
  colour= measure,
  xlab = "Weighted Mean Change in 12 weeks [95% CI]",
  name = Study, se = SE, ci = 0.95,  clip = c(-.01, 0.04),xlim=c(0,.025)) +
  theme(axis.title.y = element_text(family = 'Arial', size = 10))+ 
  theme(axis.title.x = element_text(family = 'Arial',size = 8))+ 
  theme(axis.text = element_text(family = 'Arial',size = 8)) +
  theme(plot.title = element_text(family = 'Arial', size = 12)) +
  theme(legend.position = "none") +
  scale_colour_discrete(na.translate = F)+ 
  ggforestplot::geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax, 
                                         colour = .data$measure, filled = .data$.filled), 
                            position = ggstance::position_dodgev(height = 0.5),
                            fatten=1)


Composite_metaregress<- subset(Composite_data, Composite_data$Study=="Meta-regressed")
Composite_data<- Composite_data %>%  subset(.,Composite_data$Study!="Meta-regressed")


Composite_forest<- forestplot(
  df = Composite_data, labeltext = c(Study), 
  estimate = Change,
  title= "Composite",
  colour= measure,
  name = Study, se = SE, ci = 0.95,clip = c(-.01, 0.04),xlim=c(-.01,.04)) +
  theme(axis.title.y = element_text(family = 'Arial', size = 10))+ 
  theme(axis.title.x = element_text(family = 'Arial',size = 8))+ 
  theme(legend.text = element_text(family = 'Arial',size = 8)) + 
  theme(axis.text = element_text(family = 'Arial',size = 8), axis.title.x=element_text(colour="white")) +
  theme(plot.title = element_text(family = 'Arial', size = 12)) +
  theme(legend.title = element_text(family = 'Arial', size = 8)) +
  scale_colour_discrete(na.translate = F)+ 
  ggforestplot::geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax, 
                                         colour = .data$measure, filled = .data$.filled), 
                            position = ggstance::position_dodgev(height = 0.5),
                            fatten=1)

Composite_metaregress_forest<- forestplot(
  df = Composite_metaregress, 
  estimate = Change,
  colour= measure,
  xlab = "Weighted Mean Change in 12 weeks [95% CI]",
  name = Study, se = SE, ci = 0.95,  clip = c(-.01, 0.04),xlim=c(0,.025)) +
  theme(axis.title.y = element_text(family = 'Arial', size = 10))+ 
  theme(axis.title.x = element_text(family = 'Arial',size = 8))+ 
  theme(axis.text = element_text(family = 'Arial',size = 8)) +
  theme(plot.title = element_text(family = 'Arial', size = 12)) +
  theme(legend.position = "none") +
  scale_colour_discrete(na.translate = F)+ 
  ggforestplot::geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax, 
                                         colour = .data$measure, filled = .data$.filled), 
                            position = ggstance::position_dodgev(height = 0.5),
                            fatten=1)


Cognitive_metaregress<- subset(Cognitive_data, Cognitive_data$Study=="Meta-regressed")
Cognitive_data<- Cognitive_data %>%  subset(.,Cognitive_data$Study!="Meta-regressed")


Cognitive_forest<- forestplot(
  df = Cognitive_data, labeltext = c(Study), 
  estimate = Change,
  title= "Cognitive",
  colour= measure,
  name = Study, se = SE, ci = 0.95,clip = c(-.01, 0.04),xlim=c(-.01,.04)) +
  theme(axis.title.y = element_text(family = 'Arial', size = 10))+ 
  theme(axis.title.x = element_text(family = 'Arial',size = 8))+ 
  theme(legend.text = element_text(family = 'Arial',size = 8)) + 
  theme(axis.text = element_text(family = 'Arial',size = 8), axis.title.x=element_text(colour="white")) +
  theme(plot.title = element_text(family = 'Arial', size = 12)) +
  theme(legend.title = element_text(family = 'Arial', size = 8)) +
  scale_colour_discrete(na.translate = F)+ 
  ggforestplot::geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax, 
                                         colour = .data$measure, filled = .data$.filled), 
                            position = ggstance::position_dodgev(height = 0.5),
                            fatten=1)

Cognitive_metaregress_forest<- forestplot(
  df = Cognitive_metaregress, 
  estimate = Change,
  colour= measure,
  xlab = "Weighted Mean Change in 12 weeks [95% CI]",
  name = Study, se = SE, ci = 0.95,  clip = c(-.01, 0.04),xlim=c(0,.025)) +
  theme(axis.title.y = element_text(family = 'Arial', size = 10))+ 
  theme(axis.title.x = element_text(family = 'Arial',size = 8))+ 
  theme(axis.text = element_text(family = 'Arial',size = 8)) +
  theme(plot.title = element_text(family = 'Arial', size = 12)) +
  theme(legend.position = "none") +
  scale_colour_discrete(na.translate = F)+ 
  ggforestplot::geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax, 
                                         colour = .data$measure, filled = .data$.filled), 
                            position = ggstance::position_dodgev(height = 0.5),
                            fatten=1)


Neuropsychiatric_metaregress<- subset(Neuropsychiatric_data, Neuropsychiatric_data$Study=="Meta-regressed")
Neuropsychiatric_data<- Neuropsychiatric_data %>%  subset(.,Neuropsychiatric_data$Study!="Meta-regressed")


Neuropsychiatric_forest<- forestplot(
  df = Neuropsychiatric_data, labeltext = c(Study), 
  estimate = Change,
  title= "Neuropsychiatric",
  colour= measure,
  name = Study, se = SE, ci = 0.95,clip = c(-.01, 0.04),xlim=c(-.01,.04)) +
  theme(axis.title.y = element_text(family = 'Arial', size = 10))+ 
  theme(axis.title.x = element_text(family = 'Arial',size = 8))+ 
  theme(legend.text = element_text(family = 'Arial',size = 8)) + 
  theme(axis.text = element_text(family = 'Arial',size = 8), axis.title.x=element_text(colour="white")) +
  theme(plot.title = element_text(family = 'Arial', size = 12)) +
  theme(legend.title = element_text(family = 'Arial', size = 8)) +
  scale_colour_discrete(na.translate = F)+ 
  ggforestplot::geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax, 
                                         colour = .data$measure, filled = .data$.filled), 
                            position = ggstance::position_dodgev(height = 0.5),
                            fatten=1)

Neuropsychiatric_metaregress_forest<- forestplot(
  df = Neuropsychiatric_metaregress, 
  estimate = Change,
  colour= measure,
  xlab = "Weighted Mean Change in 12 weeks [95% CI]",
  name = Study, se = SE, ci = 0.95,  clip = c(-.01, 0.04),xlim=c(0,.025)) +
  theme(axis.title.y = element_text(family = 'Arial', size = 10))+ 
  theme(axis.title.x = element_text(family = 'Arial',size = 8))+ 
  theme(axis.text = element_text(family = 'Arial',size = 8)) +
  theme(plot.title = element_text(family = 'Arial', size = 12)) +
  theme(legend.position = "none") +
  scale_colour_discrete(na.translate = F)+ 
  ggforestplot::geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax, 
                                         colour = .data$measure, filled = .data$.filled), 
                            position = ggstance::position_dodgev(height = 0.5),
                            fatten=1)

ggsave(file="multiforestplot.eps", device = cairo_ps,  dpi = "retina", width=25, height=30, units = "cm", 
       arrangeGrob(Functional_forest,Cognitive_forest,Functional_metaregress_forest, Cognitive_metaregress_forest, Composite_forest,Neuropsychiatric_forest, Composite_metaregress_forest, Neuropsychiatric_metaregress_forest, ncol=2, nrow=4, heights=c(30,6,30,4)))

