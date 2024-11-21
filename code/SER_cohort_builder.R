################################################################################
#                                Fabian Yii @ 2024                             #
#                               fabian.yii@ed.ac.uk                            #
################################################################################

## Clear workspace
rm(list=ls())

## Set working directory to parent directory
setwd("")

## Read the cleaned tabular data in long format
d  <- read.csv(file.path("data", "cleaned_data_long_all.csv"))

## Read the csv file with fundus image quality grading and
## combine it with "d" (N=137016; 68508 RE & 68508 LE)                     
quality                           <- read.csv(file.path("data", "fundus_image_quality_result.csv"))
d$fundusQuality                   <- NA 
rowsWithFundus                    <- which(d$fundus != "") # 1229 eyes don't have fundus photos
d[rowsWithFundus, ]$fundusQuality <- quality$quality

#############################################################################################
############################# Quality control & building cohort #############################
#############################################################################################

## Include only eyes with fundus photo
# 135787 eyes left (68108 RE, 67679 LE)
d <- d[d$fundus != "", ]

## Include only eyes with good or usable fundus image quality
# 90191 eyes left (47272 RE, 42919 LE)
good_or_usable_rows      <- d$fundusQuality != "Reject"
d <- d[good_or_usable_rows, ]

## Include only eyes where spherical refractive error is available
# 89216 eyes left (46819 RE, 42397 LE)
d <- d[!is.na(d$SER), ] 

## Include only eyes where corneal radius is available
# 85456 eyes left (44715 RE, 40741 LE)
d <- d[!is.na(d$meanCornealRadius), ]

## Include only eyes with reasonably normal CR
quantile(d$meanCornealRadius, probs = c(0.005, 0.995), na.rm=TRUE)
lowerCutoff   <- 7.11
upperCutoff   <- 8.63
# 395 eyes (211 RE, 184 LE)
removeLowCR   <- which(d$meanCornealRadius < lowerCutoff) 
# 434 eyes (228 RE, 206 LE)
removeHighCR  <- which(d$meanCornealRadius > upperCutoff) 
# 84627 eyes left (44276 RE, 40351 LE)
d             <- d[-c(removeLowCR, removeHighCR), ]

## Include only eyes where VA is available
# 84409 eyes left (44173 RE, 40236 LE)
d <- d[!is.na(d$VA), ] 

## Include only eyes with VA better than 0.1 LogMAR
# 49267 eyes left (25582 RE, 23685 LE)
d <- d[d$VA <= 0, ] 

## Include only indiviudlas with good systemic health
# 35300 eyes left (18268 RE, 17032 LE)
# 13967 eyes excluded due to at least one form of systemic condition (note: even split b/w RE & LE)
# Of which, most (10619 eyes) have hypertension, 
# followed by diabetes (1784 eyes),
# followed by myocardial infarction (747 eyes)
# Together (at least one of the conditions present), they (11615 eyes) account for 83.2% of removals
d <- d[rowSums(d[,23:44]) == 0, ] 

## Include only indiviudlas with good ocular health
# 34455 eyes left (17828 RE, 16627 LE)
# 845 eyes excluded due to at least one form of ocular condition in that eye or the fellow eye, or both
# Of which, most have (346 eyes) have glaucoma, 
# followed by some form of chorioretinal disorder (311 eyes),
# followed by globe/scleral disorder including degenerative myopia and staphyloma (124 eyes);
# together (at least one of these conditions is present), they (761 eyes) account for 90.1% of all removals.
# Note 1: 58 eyes have strabismus
# Note 2: 29 eyes have some form of optic nerve disorder
# Note 3: 3 eyes have nystagmus
d <- d[rowSums(d[,45:51]) == 0, ] 

## Save final data frame as csv
write.csv(d, 
          file.path("data", "cleaned_data_long_SER_cohort.csv"),
          row.names = FALSE) 


