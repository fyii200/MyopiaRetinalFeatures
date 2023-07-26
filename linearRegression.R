# Relative importance of different optic disc and foveal 
# parameters in the prediction of refractive error

# Author: Fabian SL Yii
# Email: fabian.yii@ed.ac.uk
################################################################################
## clear workspace
rm(list=ls())
## load useful packages
library(ggplot2)
library(car)
library(tidyverse)
library(dplyr)
library(segmented)
library(sjPlot)

####################################################################################
################################# Data preparation #################################
####################################################################################

## Read data
d <- read.csv('data/derivedData.csv')

## Factorise the following variables
d$sex <- as.factor(d$sex)
d$md <- as.factor(d$md)                     # macular degeneration (yes/no)
d$glaucoma <- as.factor(d$glaucoma)         # glaucoma (yes/no)
d$cataract <- as.factor(d$cataract)         # cataract (yes/no)
d$diabetic_eye <- as.factor(d$diabetic_eye) # diabetes-related eye condition (yes/no)

## Compute optic disc tilt (ratio of major axis length to minor axis length)
d$tilt <- d$major_length / d$minor_length 

## Remove eyes with self-reported macular degeneration (MD) and glaucoma
normal_d <- subset(d, md=='no' & glaucoma=='no')
n_removed <- nrow(d)-nrow(normal_d)
paste(n_removed, 'images were removed')
if(length(which(d$md=='yes' & d$glaucoma=='yes')) == 0){
  n_md <- length(which(d$md=='yes'))
  n_glaucoma <- length(which(d$glaucoma=='yes'))
  paste(n_md, 'removed due to MD, while', n_glaucoma, 'removed due to glaucoma')
}

## Remove images with reflection artifacts (unusually high foveal pixel intensity); 
## visual inspection confirmed the presence of reflection artifacts in these images
## but not the remaining ones
n_artifact <- sum(normal_d$scaled_macula_intensity > 5)
paste(n_artifact, 'removed due to reflection artifacts')
normal_d <- normal_d[normal_d$scaled_macula_intensity <= 5,]

## Remove images from participants with a distance VA poorer than 0.10 LogMAR
n_poorVA <- sum(normal_d$va > 0.1)
paste(n_poorVA, 'removed due to poor VA')
normal_d <- normal_d[normal_d$va <= 0.1,]

## Of the remaining eyes, 5 & 8 had self-reported diabetes-related condition
## and cataract, respectively. They were not excluded as these conditions were
## not judged to influence the retinal parameters considered herein.
table(normal_d$diabetic_eye)
table(normal_d$cataract)

## Final N = 273
paste('Left with' ,nrow(normal_d), 'images')


####################################################################################
################################## Data analysis ###################################
####################################################################################

## Compare included and excluded particpants
removed <- d[!(d$name %in% normal_d$name),]
for (column in c('age', 'sex', 'ser', 'adj_dist', 'adj_od_area', 'scaled_macula_intensity', 'orientation', 'tilt', 'vertical_angle')){
  print(paste('==================', column, '=================='))
  if (column == 'sex'){
    m <- cbind(table(normal_d$sex),table(removed$sex))
    print(chisq.test(m))
  }
  else{ print(t.test(normal_d[column], removed[column])) }
}

## Multiple (multivariable) linear regression (age and sex additionally as covairates)
## "scale": function that standardises the values of each predictor, i.e. mean=0 & SD=1

# All eyes included
fullModel <- lm(ser ~ scale(adj_dist) + 
                      scale(orientation) + 
                      scale(scaled_macula_intensity) + 
                      scale(age) + 
                      scale(adj_od_area) + 
                      sex + 
                      scale(vertical_angle) + 
                      scale(tilt), 
               normal_d)
tab_model(fullModel, pred.labels=c("Intercept", "OD-fovea Distance", "OD Orientation", "FPI", "Age", "OD Area", "Sex", "OD-macula angle", "OD Tilt"))

# Myopes only
myopes <- subset(normal_d, ser < -0.25)
myopiaModel <- lm(ser ~ scale(adj_dist) + 
                        scale(orientation) + 
                        scale(scaled_macula_intensity) + 
                        scale(age) + 
                        scale(adj_od_area) + 
                        sex + 
                        scale(vertical_angle) + 
                        scale(tilt), 
                  myopes)
tab_model(myopiaModel, pred.labels=c("Intercept", "OD-fovea Distance", "OD Orientation", "FPI", "Age", "OD Area", "Sex", "OD-macula angle", "OD Tilt"))

# Hyperopes only
hyperopes <- subset(normal_d, ser >= -0.25)
hyperopiaModel <- lm(ser ~ scale(adj_dist) + 
                           scale(orientation) + 
                           scale(scaled_macula_intensity) + 
                           scale(age) + 
                           scale(adj_od_area) + 
                           sex + 
                           scale(vertical_angle) + 
                           scale(tilt), 
                     hyperopes)
tab_model(hyperopiaModel, pred.labels=c("Intercept", "OD-fovea Distance", "OD Orientation", "FPI", "Age", "OD Area", "Sex", "OD-macula angle", "OD Tilt"))


## Univariable linear regression 
data <- normal_d     # all eyes included
# data <- myopes.    # myopes only
# data <- hyperopes. # hyperopes only

columns <- c(2, 8, 9, 11, 20, 28, 33, 41)
for (i in columns){
  print(names(data)[i])
  if(i == 33){                          # column 33 corresponds to 'sex'
    m <- lm(data$ser ~ data[,i])        # don't standardise the data if it's 'sex'
    print(summary(m))
    print(confint(m))
  } 
  else{ 
    m <- lm(data$ser ~ scale(data[,i])) # standardise other continuous predictors
    print(summary(m))
    print(confint(m))
  }
}

####################################################################################
######################## Additional & sensitivity analyses #########################
####################################################################################

############################# Part 1 ############################# 
## The relationship between SER and OD-fovea distance is nonlinear
# Standard linear regression
lin <- lm(adj_dist ~ ser, normal_d)
# Third-degree polynomial
nonlin <- lm(adj_dist ~ poly(ser,3), normal_d)
# Likelihood ratio test (third-degree polynomial has better fit)
anova(nonlin, lin, test='LRT')

## Segmented linear regression: which SER breakpoint is optimal?
segmented.fit <- segmented(lin, seg.Z = ~ser, psi=-2)
pred <- predict(segmented.fit)
ix <- sort(normal_d$ser, index.return=T)$ix
plotD <- data.frame('x'=normal_d$ser[ix], 'y'=pred[ix])
summary(segmented.fit) # optimal breakpoint: -2.50 D

## Plot OD-fovea distance vs SER using third-degree polynomial (black)
## and segmented linear regression (maroon) models
nonLinearPlot <- 
  ggplot(normal_d, aes(x=ser, y=adj_dist)) + 
  geom_point(alpha=0.15, size=3) +
  stat_smooth(method="lm", se=TRUE, fill=NA,
              formula=y ~ poly(x, 3, raw=TRUE),
              colour="black") +
  theme_ipsum(ticks=TRUE) +
  geom_line(data=plotD, aes(x=x, y=y, colour='maroon'), size=1.2) +
  ylim(c(150,280)) +
  scale_y_continuous( breaks = seq(150,280,20)) +
  labs(x = 'SER (D)', y='OD-fovea distance (pixels)', subtitle='n=273') +
  theme(legend.position = "none",
        panel.background = element_rect(fill="#fbf9f4", color="#fbf9f4"),
        plot.background = element_rect(fill="#fbf9f4", color="#fbf9f4"),
        panel.grid.major.x = element_line(color="grey", linetype="dashed"),
        panel.grid.minor.x = element_line(color="grey", linetype="dashed"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(0.5,0,0.5,0.5, "cm")) 
nonLinearPlot

## Sensitivity analysis:repeat the same (segmented linear regression) analysis
## including eyes that were excluded from the main analyses due to poor VA,
## as they were not judged to influence OD-fovea distance
poor_d <- subset(d, md=='no' & glaucoma=='no' & scaled_macula_intensity <= 5)
lin <- lm(adj_dist ~ ser, poor_d)
nonlin <- lm(adj_dist ~ poly(ser,3), poor_d)
# Likelihood ratio test (third-degree polynomial still has better fit)
anova(nonlin, lin, test='LRT')

# Segmented linear regression: -2.50 D is still the most optimal breakpoint
segmented.fit <- segmented(lin, seg.Z = ~ser, psi=-2)
pred <- predict(segmented.fit)
ix <- sort(normal_d$ser, index.return=T)$ix
plotD <- data.frame('x'=normal_d$ser[ix], 'y'=pred[ix])
summary(segmented.fit) 


############################# Part 2 ############################# 
## Sensitivity analysis: is foveal pixel intensity still associated with
## SER after including only eyes with good VA (VA < 0) in myopia?

# Myopia
myopeModelGoodVA <- lm(ser ~ scale(adj_dist) + 
                             scale(orientation) + 
                             scale(scaled_macula_intensity) + 
                             scale(age) + 
                             scale(adj_od_area) + 
                             sex + 
                             scale(vertical_angle) + 
                             scale(tilt), 
                             subset(myopes, va <= 0) )
tab_model(myopeModelGoodVA, pred.labels=c("Intercept", "OD-fovea Distance", "OD Orientation", "FPI", "Age", "OD Area", "Sex", "OD-macula angle", "OD Tilt"))

# Hyperopia
hyperopeModelGoodVA <- lm(ser ~ scale(adj_dist) + 
                                scale(orientation) + 
                                scale(scaled_macula_intensity) + 
                                scale(age) + 
                                scale(adj_od_area) + 
                                sex + 
                                scale(vertical_angle) + 
                                scale(tilt), 
                          subset(hyperopes, va <= 0) )
tab_model(hyperopeModelGoodVA, pred.labels=c("Intercept", "OD-fovea Distance", "OD Orientation", "FPI", "Age", "OD Area", "Sex", "OD-macula angle", "OD Tilt"))


############################# Part 3 ############################# 
## Association between foveal pixel intensity and SER in the red, 
## green & blue channels (multivariable linear regression) in myopes

# Red channel
myopeGoodVAModelR <- lm(ser ~ scale(adj_dist) + 
                              scale(orientation) + 
                              scale(scaled_macula_intensity_R) + 
                              scale(age) + 
                              scale(adj_od_area) + 
                              sex + 
                              scale(vertical_angle) + 
                              scale(tilt), 
                        subset(myopes, va <= 0) )
tab_model(myopeGoodVAModelR, pred.labels=c("Intercept", "OD-fovea Distance", "OD Orientation", "FPI_red", "Age", "OD Area", "Sex", "OD-macula angle", "OD Tilt"))

# Green channel
myopeGoodVAModelG <- lm(ser ~ scale(adj_dist) + 
                              scale(orientation) + 
                              scale(scaled_macula_intensity_G) + 
                              scale(age) + 
                              scale(adj_od_area) + 
                              sex + 
                              scale(vertical_angle) + 
                              scale(tilt), 
                        subset(myopes, va <= 0) )
tab_model(myopeGoodVAModelG, pred.labels=c("Intercept", "OD-fovea Distance", "OD Orientation", "FPI_green", "Age", "OD Area", "Sex", "OD-macula angle", "OD Tilt"))

# Blue channel
myopeGoodVAModelB <- lm(ser ~ scale(adj_dist) + 
                              scale(orientation) + 
                              scale(scaled_macula_intensity_B) + 
                              scale(age) + 
                              scale(adj_od_area) + 
                              sex + 
                              scale(vertical_angle) + 
                              scale(tilt), 
                        subset(myopes, va <= 0) )
tab_model(myopeGoodVAModelB, pred.labels=c("Intercept", "OD-fovea Distance", "OD Orientation", "FPI_blue", "Age", "OD Area", "Sex", "OD-macula angle", "OD Tilt"))














