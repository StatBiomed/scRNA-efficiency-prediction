# Prediction of RNA splicing kinetics in scRNA-seq data

This repository contains 
1) scripts for feature extraction on introns and genes, 
2) notebooks for estimating splicing and degradation rates, 
3) notebooks for prediction experiments
4) processed data, including features and estimated splicing rates, which can be directly used as validation for machinie learning model development.


## Project Organization

    ├── README.md          <- Top-level README.md
    ├── data
    │   ├── featuers       <- Featreus at both intron and gene levels
    │   ├── estimated      <- Estimated splicing and degradation rates
    │   ├── predicted      <- Predicted splicing and degradation rates
    │   └── raw            <- Raw data (optional)
    │
    ├── notebooks          <- Jupyter notebooks
    │   └── figures        <- Figure generation for the paper
    │   ├── estimate       <- Estimating splicing and degradation rates
    │   └── predicte       <- Predicting splicing and degradation rates
    │
    └── src                <- Source code for use in this project
        ├── preprocess     <- Scripts to preprocess the data (incl. feature extraction and combination)
        └── models         <- Scripts to train models (optional)