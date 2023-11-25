###
# In RE more negative orientation means the disc is increasingly 
# tilted towards fovea, while in LE this means the disc is increasingly
# tilted away from fovea
###
library(vctrs)
library(sjPlot)
library(segmented)
library(ggplot2)
library(hrbrthemes)
library(car)
library(lme4)
library(quantreg)
library(RColorBrewer)
library(dplyr)
library(colorspace)
library(CGPfunctions)


rm(list=ls())
setwd("/Users/fabianyii/OneDrive - University of Edinburgh/Projects/UKB_full")
# setwd('../') # set working directory to project parent directory

# N=34455 (17828 RE, 16627 LE)
d <- read.csv("data/cleaned_data_long_SER_cohort.csv")

# Derived (imaging) data (89216 eyes)
fullDerivedData <- read.csv("outputs/csv/fullDerivedData.csv")

# N=34455 (17828 RE, 16627 LE)
d <- merge(d, fullDerivedData, by="fundus") 

# Convert FPI to numeric
d$scaled_fovea30r_intensity_gray <- as.numeric(d$scaled_fovea30r_intensity_gray)

# Calculate OD area (https://www.sciencedirect.com/science/article/pii/S0161642095309050)
d$adj_od_area_ellipse <- d$adj_minor_length * d$adj_major_length * pi/4

# Adjust CRAE & CRVE for ocular magnification
littmann <- function(measured_size, SER, CR){
  # Function: calculates object size by accounting for the effect of 
  # magnification due to ametropia. Takes image_size, spherical 
  # equivalent refraction (SER) and mean (b/w strong and week meridians) 
  # corneal radius (CR, in mm), and returns true (object) size. This is based
  # on Littmann's original magnification formula for fundus camera.
  # READ: https://link.springer.com/content/pdf/10.1007/BF00175988.pdf
  a = 0.01 + 0.00236 * (CR - 8)
  b = 0.6126 + 0.0968 * (CR - 8)
  c = 30.52 + 2.57 * (CR - 8)
  q = (a*SER^2 - b*SER + c) / 100
  true_size = 1.37 * q * measured_size
  return(true_size) }
d$adj_CRAE_Knudtson <- littmann(d$CRAE_Knudtson, d$SER, d$meanCornealRadius)
d$adj_CRVE_Knudtson <- littmann(d$CRVE_Knudtson, d$SER, d$meanCornealRadius)

# # Bengston
# # https://www.sciencedirect.com/science/article/pii/S0002939407006630?via%3Dihub
# d$adj_CRAE_Knudtson <- d$CRAE_Knudtson * (1 - 0.017 * d$SER)
# d$adj_CRVE_Knudtson <- d$CRVE_Knudtson * (1 - 0.017 * d$SER)
# d$adj_major_length <- d$major_length * (1 - 0.017 * d$SER)
# d$adj_minor_length <- d$minor_length * (1 - 0.017 * d$SER)
# d$adj_dist <- d$dist * (1 - 0.017 * d$SER)

# Factorise id, eye and sex
d$id <- factor(d$id)
d$eye <- factor(d$eye)
d$sex <- factor(d$sex)

####################################################################################
################################## Start analyses ##################################
####################################################################################

## RE: Linear quantile regression 
taus <- seq(0.005, 0.995, 0.03)
# taus <- seq(0.01, 0.99, 0.02)
n_taus <- length(taus)

RE <- rq(SER ~ sex +
           scale(age) + 
           scale(adj_dist) + 
           scale(vertical_angle) +
           scale(orientation) + 
           scale(scaled_fovea30r_intensity_gray) + 
           scale(adj_major_length/adj_minor_length) +
           scale(adj_od_area_ellipse) +
           scale(adj_CRAE_Knudtson) +
           scale(adj_CRVE_Knudtson) +
           scale(Tortuosity_density_combined) +
           scale(FD_combined) +
           scale(conc_rp_artery) +
           scale(conc_rp_vein),
         na.action = na.omit,
         tau = taus,
         data = subset(d, eye=="RE") )
REsummary <- summary(RE)

## RE: OLS linear regression (myopes)
RE_OLS_myopes <- lm(SER ~ sex +
                      scale(age) + 
                      scale(adj_dist) + 
                      scale(vertical_angle) +
                      scale(orientation) + 
                      scale(scaled_fovea30r_intensity_gray) + 
                      scale(adj_major_length/adj_minor_length) +
                      scale(adj_od_area_ellipse) +
                      scale(adj_CRAE_Knudtson) +
                      scale(adj_CRVE_Knudtson) +
                      scale(Tortuosity_density_combined) +
                      scale(FD_combined) +
                      scale(conc_rp_artery) +
                      scale(conc_rp_vein),
                    data = subset(d, eye=="RE" & SER < -0.25) )
RE_OLS_myopes_summary <- summary(RE_OLS_myopes)

## RE: OLS linear regression (non-myopes)
RE_OLS_nonmyopes <- lm(SER ~ sex +
                         scale(age) + 
                         scale(adj_dist) + 
                         scale(vertical_angle) +
                         scale(orientation) + 
                         scale(scaled_fovea30r_intensity_gray) + 
                         scale(adj_major_length/adj_minor_length) +
                         scale(adj_od_area_ellipse) +
                         scale(adj_CRAE_Knudtson) +
                         scale(adj_CRVE_Knudtson) +
                         scale(Tortuosity_density_combined) +
                         scale(FD_combined) +
                         scale(conc_rp_artery) +
                         scale(conc_rp_vein),
                       data = subset(d, eye=="RE" & SER >= -0.25) )
RE_OLS_nonmyopes_summary <- summary(RE_OLS_nonmyopes)

## LE: Linear quantile regression 
LE <- rq(SER ~ sex +
           scale(age) + 
           scale(adj_dist) + 
           scale(vertical_angle) +
           scale(orientation) + 
           scale(scaled_fovea30r_intensity_gray) + 
           scale(adj_major_length/adj_minor_length) +
           scale(adj_od_area_ellipse) +
           scale(adj_CRAE_Knudtson) +
           scale(adj_CRVE_Knudtson) +
           scale(Tortuosity_density_combined) +
           scale(FD_combined) +
           scale(conc_rp_artery) +
           scale(conc_rp_vein),
         na.action = na.omit,
         tau = taus,
         data = subset(d, eye=="LE"))
LEsummary <- summary(LE)

## LE: OLS linear regression (myopes)
LE_OLS_myopes <- lm(SER ~ sex +
                      scale(age) + 
                      scale(adj_dist) + 
                      scale(vertical_angle) +
                      scale(orientation) + 
                      scale(scaled_fovea30r_intensity_gray) + 
                      scale(adj_major_length/adj_minor_length) +
                      scale(adj_od_area_ellipse) +
                      scale(adj_CRAE_Knudtson) +
                      scale(adj_CRVE_Knudtson) +
                      scale(Tortuosity_density_combined) +
                      scale(FD_combined) +
                      scale(conc_rp_artery) +
                      scale(conc_rp_vein),
                    data = subset(d, eye=="LE" & SER < -0.25) )
LE_OLS_myopes_summary <- summary(LE_OLS_myopes)

## LE: OLS linear regression (non-myopes)
LE_OLS_nonmyopes <- lm(SER ~ sex +
                         scale(age) + 
                         scale(adj_dist) + 
                         scale(vertical_angle) +
                         scale(orientation) + 
                         scale(scaled_fovea30r_intensity_gray) + 
                         scale(adj_major_length/adj_minor_length) +
                         scale(adj_od_area_ellipse) +
                         scale(adj_CRAE_Knudtson) +
                         scale(adj_CRVE_Knudtson) +
                         scale(Tortuosity_density_combined) +
                         scale(FD_combined) +
                         scale(conc_rp_artery) +
                         scale(conc_rp_vein),
                       data = subset(d, eye=="LE" & SER >= -0.25) )
LE_OLS_nonmyopes_summary <- summary(LE_OLS_nonmyopes)


## Save results in a data frame
extract_results <- function(summaryObject){
  
  # Names of features to be extracted
  features <- c("Intercept", "Male", "Age", "OD-fovea distance", "OD-fovea angle", 
                "OD orientation", "FPI", "OD ovality", "OD area", "CRAE",
                "CRVE", "Vessel tortuosity", "Vessel FD", "Artery concavity", 
                "Vein concavity")
  
  n_features <- length(features)
  
  # Create empty dataframe
  df <- data.frame("features"=vec_rep_each(features, n_taus),
                   "taus"=vec_rep(taus, n_features),
                   "coef"=rep(NA, n_features*n_taus),
                   "SE"=rep(NA, n_features*n_taus),
                   "p_val"=rep(NA, n_features*n_taus))
  
  # Start extracting results
  for(feature in features){
    for(tau in taus){
      df_row <- (df$features == feature) & (df$taus == tau)
      tau_ind <- match(tau, taus)
      feature_ind <- match(feature, features)
      df[df_row,]$coef <- summaryObject[tau_ind][[1]]$coefficient[feature_ind, 1]
      df[df_row,]$SE <- summaryObject[tau_ind][[1]]$coefficient[feature_ind, 2]
      df[df_row,]$p_val <- summaryObject[tau_ind][[1]]$coefficient[feature_ind, 4]
    }
  }
  # Return dataframe
  return(df)
}

REdf <- extract_results(REsummary)
REdf$eye <- "RE"
LEdf <- extract_results(LEsummary)
LEdf$eye<- "LE"

df <- rbind(REdf, LEdf)
df$eye <- factor(df$eye)

## Specify background and foreground (plot panel) colours
bg_col <- rgb(0.6,0.1,0, alpha=0.01)
fg_col <- rgb(0.8,0.4,0, alpha=0.03)

## Intercept
df %>% filter(features == "Intercept") %>% ggplot(aes(x=taus, y=coef)) +
  geom_segment( aes(x=taus, xend=taus, y=0, yend=coef), colour="gray", alpha=0.6) +
  geom_point(size=2, colour=rep(hcl.colors(n_taus),2)) +
  facet_wrap(~eye, nrow=1, scales="fixed") +
  theme_light() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm") ) +
  geom_hline(yintercept=0, size=0.5, linetype="dashed", alpha=0.3, colour="blue") +
  geom_vline(xintercept=0.515, linetype="dashed", size=0.5, alpha=0.3, colour="blue") +
  ylab("SER (D)") +
  xlab("Refractive quantile") +
  scale_x_continuous(labels = seq(0, 1, 0.1), breaks = seq(0,1,0.1)) +
  scale_y_continuous(labels = seq(-7.5, 5, 1), breaks = seq(-7.5, 5, 1))
ggsave("intercept.png", width=7, height=6, units="in", bg="white")

## Coefficient plots (BE)
plot_df <- filter(df, !features %in% c("Age", "Male", "Intercept"))
point_cols <- ifelse(plot_df$p_val>0.05, "black", NA)
ggplot(data=plot_df, aes(x=taus, y=coef)) + 
  geom_point(size=1, alpha=1, aes(colour=eye) ) +
  geom_point(size=1, alpha=1, colour=point_cols) +
  geom_smooth(size=0.5, alpha=0.1, se=FALSE, aes(group=eye, colour=eye)) +
  geom_blank(aes(y = 0)) +
  geom_ribbon(aes(ymin=coef-SE*1.96, ymax=coef+SE*1.96, fill=eye), linetype=2, alpha=0.1) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5, alpha=0.25) +
  geom_vline(xintercept=0.5, size=0.5, alpha=0.3, colour="darkgreen") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  xlab("Refractive quantile") +
  ylab("Standardised beta") +
  facet_wrap(~features, nrow=3, scales="free_y") +
  theme_test() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.margin.y = unit(0.8, "lines"),
        plot.background = element_rect(fill=bg_col),
        strip.background = element_rect(colour=bg_col, fill=fg_col),
        panel.background = element_rect(fill=fg_col),
        strip.text = element_text(face="bold"),
        legend.position = "top",
        legend.title = element_blank()) +
  scale_color_manual(values=c("RE"="red", "LE"="blue")) +
  scale_fill_manual(values=c("RE"="red", "LE"="blue")) +
  labs(caption = "*Black points: p > 0.05")
ggsave("coef.png", width=9, height=9, units="in", bg="white")

## Rank features by magnitude of association
filter(df, taus %in% c(0.005, 0.095, 0.215, 0.305, 0.395, 0.515, 0.605, 0.725, 0.815, 0.905, 0.995) & !features %in% c("Age", "Male", "Intercept")) %>% 
  mutate(taus = factor(taus)) %>% 
  # mutate(taus = factor(taus, labels = c("High myopia","Myopia","Emmetropia","Hyperopia", "High hyperopia"))) %>% 
  mutate(rounded_coef = round(abs(coef),2)) %>% 
  newggslopegraph(taus, rounded_coef, features, 
                  Title=NULL, 
                  ThemeChoice="wsj",
                  Caption="Absolute beta vs refractive quantile",
                  SubTitle=NULL,
                  WiderLabels=TRUE,
                  YTextSize=3,
                  XTextSize = 8,
                  CaptionJustify=0,
                  CaptionTextSize=10,
                  LineThickness=0.8) +
  facet_wrap(~eye, nrow=2, scales="fixed") +
  theme(plot.margin=margin(0.5, 0.5, 0.5, 0.5, "cm"),
        strip.text = element_text(face="bold"))
ggsave("features_ranked.png", width=8, height=7, units="in", bg="white")











  
