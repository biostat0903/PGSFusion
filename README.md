# PGSFusion <br>

## Overview
We develop PGSFusion to perform polygenic score (PGS) construction and epidemiolgical application. PGSFusion accepted many different formats of summary statistics and implements a wide range of PGS methods, including single-trait, multiple-trait, annotation-based, and cross-ancestry methods. PGSFusion provided further epidemiological analysis: i) prediction performance evaluation; and ii) interaction identification. <br>
The user can use the server to visit the website: [http://www.pgsfusion.net/#/](http://www.pgsfusion.net/#/)

## PGSFusion pipeline
### UKBB data
+ We select 138 quantitative traits (i.e. 61 Anthropometric traits, 65 Blood and Unrine Test traits and 12 Other traits) and 63 binary traits (i.e. 14 Cancers, 44 Non-cancer Diseases, and 5 Other traits). <br>
+ We used UKBB data to build the reference panel, validation set and test set. For EUR reference panels, we used 2,000 EUR individuals in UKBB. For AFR and EAS, we used corresponding individuals in 1000 Genome Project. For EUR test set, we selected 50,000 individuals in UKBB. For AFR and EAS, we selected all the corresponding samples in UKBB. <br>
### PGS method
+ We integrated 17 PGS methods, including 11 single-trait, 1 multiple-trait, 2 annotation-based and 3 cross-ancestry based methods. <br>
+ We constructed the reference panel for each method.
### PGS application
+ We output the linkage for the selected methods
+ We provide two further applications: i) prediction performance evaluation; ii) interaction analysis. <br>

## Update log
### Version 3.0 (Jan 2025)
+ Add the number of quantitative and binary traits <br>
+ Add SBayesRC method <br>
+ Add the comparison of prediction performance for selected methods <br>
+ Use 2,000 EUR individuals to build EUR reference panel <br>
### Version 2.0 (May 2024)
+ Update the PGSFusion in formating summary statistics and performing epidemiological applications
### Version 1.0 (Dec 2023)
+ Build the PGSFusion to streamline the PGS construction

## Citation
If you used our server, please cite the our paper: <br>
<em><strong>Yang S</strong></em>\#\*, Ye X\#, Ji X\#, Li Z, Tian M, Huang P\*, Cao C\*. PGSFusion streamlines polygenic score construction and epidemiological applications in biobank-scale cohorts. <em>bioRxiv</em>. 
