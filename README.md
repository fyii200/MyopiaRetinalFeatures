## This repository contains the source code for the following [work](https://iovs.arvojournals.org/article.aspx?articleid=2793667):
```
Yii F, Bernabeu MO, Dhillon D, Strang N, MacGillivray T. Retinal changes from hyperopia to myopia: Not all dioptres are created equal. Invest. Ophthalmol. Vis. Sci. 2024;65(5):25. https://doi.org/10.1167/iovs.65.5.25.
```

An overview of the workflow is provided below (source code in the directory named "code").

### Step 1: Participant selection
<pre>
SER_cohort_builder.R : R script for selecting eligible eyes/participants for the study.
</pre>

<p align="center">
  <img src="https://github.com/user-attachments/assets/198eed3c-606b-429f-bb60-964c1917c647" width="550" />
</p>

### Step 2: Segmentation of regions of interest
<pre>
OD_segmentation : Folder containing a custom-trained deep learning model (DeepLabV3 with a MobileNetV3 large backbone) used for optic disc segmentation. 

Optic disc segmentation model is freely and openly available to any interested reaserchers, provided that appropriate credit is given (see Note 1).
</pre>
Note 1: More details about the disc segmentation model are available at [Yii et al.](https://doi.org/10.1007/978-3-031-47425-5_30). Model weights are available to download from this repository. 

Note 2: Artery/vein segmentation was done using AutoMorph, which is freely and openly available [elsewhere](https://github.com/rmaphoh/AutoMorph/tree/main)

Note 3: Foveal segmentation/localisation was done using a DU-Net model available from [VAMPIRE](https://vampire.computing.dundee.ac.uk/index.html)



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








