This R project implements a test for checking whether the simplifying assumption is reasonable for a given dataset.\\
In the R folder the code for the project can be found. \\
The data folder contains the datasets required for the rest of the script. \\
The models folder contains models, that have been trained on the data, which include (simplified) vine copula models and the models used for classifying the real versus the simplified data. \\
The scripts in the R folder have the following functionality: 
 - simulate_non_simplified_vine.R contains the functions, that allow to simulate from a non simplified vine copula.
 - simulate_data.R contains the code, that actually simulates non-simulated vine copula data.
 - classifier.R trains a Neural Network for classification to distinguish between
 - quantile_reg.R fits a D-vine quantile regression model to the values $r_{i}^n := ln(c_{NP}(u_i) / c_{\hat{RV}}(u_i))$ and performs the test for the simplifying assumption.
