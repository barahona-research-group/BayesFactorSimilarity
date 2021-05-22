Similarity measure for sparse time course data with Gaussian processes
================
Zijing Liu
2021-05-22

Introduction
------------
This repository contains MATLAB functions for modelling time course data with
Gaussian processes (GP) and computing a pair-wise similarity measure in the form of
a Bayes factor. It uses the GPML toolbox (http://www.gaussianprocess.org/gpml/code/matlab/doc/).

* BF_onehyp.m - a function for computing the pair-wise similarity matrix, where the hyperparameters are optimised for the whole dataset.
* BF_twohyp.m - a function for computing the pair-wise similarity matrix, where the hyperparameters are optimised for each pair of time courses.
* BF_async.m - a function for computing the pair-wise similarity matrix, where the hyperparameters are optimised for the whole dataset and the time courses are asynchronous.
* gpRegression.m - a function to do 1-D GP regression.
* cluster_dist_matrix.m - a function to do clustering given a distance matrix.

* lib/ - the required Matlab packages including:
	* GPML toolbox for Gaussian process
	* Ncut and ZPclustering for spectral clustering
	* InfoTheory toolbox for NMI

* R/ - contains the R code to compute the BHI z-score.
* gene_data.mat - it is the Matlab data file containing the gene expression data.

References
------------
Liu, Zijing, and Mauricio Barahona. "Similarity measure for sparse time course data based on Gaussian processes." arXiv preprint arXiv:2102.12342 (2021).