## Relative importance of different optic disc and foveal parameters in the prediction of refractive error
Anonymised imaging data derived from colour fundus photographs and other relevant variables can be found in the 'data' directory (*derivedDataAnonymised.csv*). The dataset contans 374 rows and 40 columns. Each row corresponds to an eye from a unique individual. The codebook below summarises the content of the dataset:

| Column name | Data type | What it represents |
| :---:   | :---: | :---: |
| **name** | Integer | Each integer represents one of the 28 studies (papers) included in the meta-regression.   |
| **age** | Integer | Each integer represents one of the 26 unique datasets used by the included studies. Martinez et al. and Philip et al. worked with Sydney Myopia Study data, while Li et al. and Li et al. worked with Anyang Childhood Eye Study data. Data from different age groups were used. |
| **ser** | Character | Last name of first author and year of publication. |
| **cr** | Character | Full title.  |
| **sph** | Character | Country or city in which the participants were recruited. Note: "hk" refers to Hong Kong.  |
| **cyl** | Character | Study design where "cs" refers to cross-sectional, while "long" refers to longitudinal. |
| **dist** | Integer | Number of eyes. |
| **adj_dist** | Integer | Number of females (if provided). |
| **vertical_angle** | Integer | Number of males (if provided). |
| **od_area** | Numeric | Mean axial length (specific to emmetropes).  |
| **adj_od_area** | Numeric | Standard deviation of axial length (specific to emmetropes) |
| **major_length** | Numeric | Mean sample age (specific to emmetropes)   |
| **adj_major_length** | Numeric | Standard deviation of sample age (specific to emmetropes)   |
| **minor_length** | Character | Optical biometry: IOLMaster or Lenstar |
| **adj_minor_length** | Boolean (y/n) | If "y", cycloplegic refraction was performed; if "n" non-cycloplegic refraction was performed.   |
| **macula_intensity** | Numeric | Mean spherical equivalent refraction, SER (specific to emmetropes, if provided).  |
| **macula_intensity_R** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **macula_intensity_G** | Numeric (range) | Definition of emmetropia: ±0.50D SER; -0.25D to +0.50D SER; -0.25D to +0.75D SER; -0.25D to +1.25D SER; -0.50D to +0.75D SER; -0.50D to +1.00D SER; -0.50D to +1.25D; ±0.50D spherical power. |
| **macula_intensity_B** | Integer | 1 to 8 where each number corresponds to each definition of emmetropia (**emm_def**) following the order specified above. |
| **scaled_macula_intensity** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **scaled_macula_intensity_R** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **scaled_macula_intensity_G** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **scaled_macula_intensity_B** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **median_intensity_R** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **median_intensity_G** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **median_intensity_B** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **median_intensity** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **orientation** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **disc_x** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **disc_y** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **macula_x** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **macula_y** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **sex** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **VA** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **cataract** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **diabetic_eye** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **md** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **ethnicity** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **glaucoma** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |
| **eye** | Numeric | Standard deviation of spherical equivalent refraction, SER (specific to emmetropes, if provided).   |


