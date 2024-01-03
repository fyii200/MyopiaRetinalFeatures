## Intro

Left eye (simulated)|Right eye (simulated)
--|--
<img src="videos/simulated_LE.gif" width="450" />|<img src="videos/simulated_RE.gif" width="450" />

Left eye (average mask)|Right eye (average mask)
--|--
<img src="videos/average_LE.gif" width="450" />|<img src="videos/average_RE.gif" width="450" />

Left eye (top eigenvector)|Right eye (top eigenvector)
--|--
<img src="videos/eigen1_LE.gif" width="450" />|<img src="videos/eigen1_RE.gif" width="450" />


### 1) ***SER_cohort_builder.R***
R script detailing each step of the participant selection process.

### 2) ***ODfovea_analysis.m***
MATLAB script computing optic disc (OD) and foveal parameters, which include OD major axis length, OD minor axis length, OD orientation, OD-fovea distance and OD-fovea angle. Note that OD major and minor axis lengths are used to compute OD area (carried out in line 36 of *linear_and_quantile_regression.R* below).

### 3) ***linear_and_quantile_regression.R***
R script used to perform multiple linear regression and quantile regression, with retinal parameters as the independent variables and spherical equivalent refraction as the dependent variable.
