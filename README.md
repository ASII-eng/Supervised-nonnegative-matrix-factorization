Thank you for downloading this MATLAB package for NMTF classification. 

This folder contains Matlab code related to the following paper:

Asieh Amousoltani Arani, Mohammadreza Sehhati and Mohammad Amin Tabatabaiefar.  "Genetic Variant Effect Prediction by Supervised Nonnegative Matrix Tri-factorization". 

The code has been tested on MATLAB 7.0.a.


Please to execute the sNMTF algorithm, run the  `NMTF_Classification_Main.m' . Please send bug reports, comments, or questions to soltani.asieh@gmail.com.

Please pay attention when you run the code, you need more than 16GB of RAM.

=================================
DESCRIPTION
=================================

The folder contains: 

Datasets: 

training data: TRAIN.txt
testing data: TEST.txt
variant-disease train network:VD_tr.txt
variant-disease test network:VD_ts.txt
disease-disease network:Disease_Disease.txt


Functions:

NMTF_Classification_Main.m: this code, is the Main function of our algorithm, supervised nonnegative matrix factorization(sNMTF) which get the data and implement the algorithm on it.   
 

factorization_Classification.m: this function optimizes the sNMTF  objective function and give back the factorization matrix of train data.

Variant_Varient_Network.m: construct the Variant-varinat network 

block_matrices.m: construct the relation and intratype matrices from input data(.txt) 

matrix_initialization.m: initilize the factorized matrices. 


generic_random_forests.m : implement the random forest classification 

 









