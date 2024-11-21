## This repository contains the source code for the following [work](https://iovs.arvojournals.org/article.aspx?articleid=2793667):
```
Yii F, Bernabeu MO, Dhillon D, Strang N, MacGillivray T. Retinal changes from hyperopia to myopia: Not all dioptres are created equal. Invest. Ophthalmol. Vis. Sci. 2024;65(5):25. https://doi.org/10.1167/iovs.65.5.25.
```

An overview of the workflow is provided below (source code int he directory named "code").

### Step 1: ***SER_cohort_builder.R***
R script detailing each step of the participant selection process.

### Step 2: ***ODfovea_analysis.m***
MATLAB script for computing optic disc (OD) and foveal parameters, which include OD major axis length, OD minor axis length, OD orientation, OD-fovea distance and OD-fovea angle. Note that OD major and minor axis lengths are used to compute OD area (in line 36 of *linear_and_quantile_regression.R* below).

### Step 3: ***linear_and_quantile_regression.R***
R script used to perform multiple linear regression and quantile regression, with retinal parameters as the independent variables and spherical equivalent refraction as the dependent variable, controlling for age, sex and corneal radius of curvature.

### Step 4: ***arcade_analysis***
Directory containing python scripts to compute vessel concavity:
##### *preprocess.py* is used to preprocess the vessel mask (artery and vein separately) and extract/detect the papillomacular vascular arcade
##### *arcade_model.py* is used to fit a second-degree polynomial function (parabola) to the preprocessed vascular arcade using either the least squares or RANSAC method. Note that RANSAC was used in this work due to its robustness to outliers.
##### *main.py* is the main script that calls *preprocess.py* and *arcade_model.py* to compute vessel concavity automatically.


#### Simulations using the mean value of each retinal parameter* stratified by refractive error (from high hyperopia to high myopia)
##### *Red vascular arcade: artery; green vascular arcade: vein; broken line: major axis of the optic disc*
##### *adjusted for ocular magnification where necessary
Left eye |Right eye 
--|--
<img src="videos/simulated_LE.gif" width="450" />|<img src="videos/simulated_RE.gif" width="450" />

#### Average segmentation mask* stratified by refractive error (from high hyperopia to high myopia)
##### *ocular magnification cannot be accounted for
Left eye |Right eye
--|--
<img src="videos/average_LE.gif" width="450" />|<img src="videos/average_RE.gif" width="450" />








