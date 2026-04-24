# PGSFusion <br>

## Overview
We develop PGSFusion to perform polygenic score (PGS) construction and epidemiological application. PGSFusion accepts many different formats of summary statistics and implements a wide range of PGS methods, including single-trait, multiple-trait, annotation-based, and cross-ancestry methods. PGSFusion provided further epidemiological analysis: i) prediction performance; and ii) joint analysis. <br>
The user can use the server to visit the website: [http://www.pgsfusion.net/#/](http://www.pgsfusion.net/#/)

## PGSFusion pipeline
### UKBB data
+ We select 138 quantitative traits (i.e. 61 Anthropometric traits, 65 Blood and Unrine Test traits and 12 Other traits) and 63 binary traits (i.e. 14 Cancers, 44 Non-cancer Diseases, and 5 Other traits). <br>
+ We use three variables at baseline to assess the SES of each participant at individual level: including family income level (data field: p738_i0), education qualification (data field: p6138_i0), employment status (data field: p6142_i0).
+ We include information on four healthy lifestyle factors collected at baseline, including “no current smoking”, “regular physical activity”, “healthy diet”, and “no alcohol consumption”.
+ We used UKBB data to build the reference panel, validation set and test set. For EUR reference panels, we used 2,000 EUR individuals in UKBB. For AFR and EAS, we used corresponding individuals in 1000 Genome Project. For EUR test set, we selected 50,000 individuals in UKBB. For AFR and EAS, we selected all the corresponding samples in UKBB. <br>
### reference panel
+ We make refenence panel for different methods using 2,000 EUR individuals from UKBB. <br>
+ We make the reference panel for EUR-AFR and EUR-EAS for SDPRX. <br>
### PGS method
+ We integrated 17 PGS methods, including 11 single-trait, 1 multiple-trait, 2 annotation-based and 3 cross-ancestry based methods. <br>
+ We constructed the reference panel for each method.
### PGS application
+ We output the linkage for the selected methods
+ We provide two further applications: i) performance evaluation; ii) joint analysis. <br>

## Update log
### Version 4.1 (June 2025)
+ Add watermark for each figure
+ Add notice for users
### Version 4.0 (Apr 2025)
+ Create new Conda environment <em>pgscalc2</em>. <br>
+ Use top 18 genetic PC to correct population structure. <br>
+ Change the EAS test set with 1419 individuals. <br>
+ Add the linkage of each method in Manual page. <br>
### Version 3.0 (Jan 2025)
+ Add the number of quantitative and binary traits <br>
+ Update [DBSLMM](https://github.com/biostat0903/DBSLMM) (v1.0) <br>
+ Use 11 strata variables <br>
+ Add SBayesRC method <br>
+ Add the comparison of prediction performance for selected methods <br>
+ Use 2,000 EUR individuals to build EUR reference panel <br>
### Version 2.0 (May 2024)
+ Update the PGSFusion in formating summary statistics and performing epidemiological applications
### Version 1.0 (Dec 2023)
+ Build the PGSFusion to streamline the PGS construction

## Citation
If you used our server, please cite the our paper: <br>
<em><strong>Yang S</strong></em>\#\*, Ye X\#, Ji X\#, Li Z, Tian M, Huang P, Cao C. PGSFusion streamlines polygenic score construction and epidemiological applications in biobank-scale cohorts. <em>Genome Med</em>. 2025 Jul 14;17(1):77. [doi: 10.1186/s13073-025-01505-w](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-025-01505-w). 



