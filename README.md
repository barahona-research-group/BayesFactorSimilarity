Similarity measure for sparse time course data with Gaussian processes
================
Zijing Liu
2018-11-08

Introduction
------------
This repository contains MATLAB functions for modelling time course data with
Gaussian processes (GP) and computing a pair-wise similarity measure in the form of
the Bayes factor

*  BF_onehyp.m - a function for computing the pair-wise similarity matrix, where the hyperparameters are optimised for the whole dataset.
*  BF_twohyp.m - a function for computing the pair-wise similarity matrix, where the hyperparameters are optimised for each pair of time courses.
*  gpRegression.m - a function to do 1-D GP regression.
*  gpml-matlab-version-number/ - the GPML toolbox from http://www.gaussianprocess.org/gpml/code/matlab/doc/

