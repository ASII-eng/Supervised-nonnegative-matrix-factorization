clc;
clear all;
close all;
k = [7, 15];
max_iter = 120;
initialization = 'random_acol'; %nnmf,random_acol,random,kmeans
gamma1=.7; 
gamma2=.9; 


%% Select the prediction scores
Score_names={'SIFT_score', 'Polyphen2_HDIV_score','Polyphen2_HVAR_score','MutationAssessor_score','PROVEAN_score','GERP++_score'...
    'phyloP20way_mammalian','phastCons20way_mammalian','SiPhy_29way_logOdds'};

%% reading data

 [num_train,~,row_train] = xlsread('TEST.xlsx');
 [num_test,~,row_test] = xlsread('TRAIN.xlsx');
 
 TRUE_Class_train=row_train(2:end,end);

%% Constructing Variant-Variant Network 


[VAR_VAR_train]=Variant_Varient_Network(num_train,row_train,Score_names);
[VAR_VAR_test]=Variant_Varient_Network(num_test,row_test,Score_names);

%% Input data: relation matrix(variant-disease), variant-variant and disease-disease networks

 adj_list_train = {VAR_VAR_train,'Disease_Disease.txt'};
 rel_file_train = {'VD_tr.txt'};

 adj_list_test = {VAR_VAR_test,'Disease_Disease.txt'};
 rel_file_test = {'VD_ts.txt'};

%% create block matrices

[R_train, A_train, label_list_train] = block_matrices(adj_list_train, rel_file_train);
[R_test, A_test, label_list_test] = block_matrices(adj_list_test, rel_file_test);

clear VAR_VAR_train VAR_VAR_test num_train num_test row_train 
%% run sNMTF and export low rank matrices

[S,G,RE] = factorization_Classification(R_train, A_train,k,max_iter,initialization,gamma1,gamma2);

%% construct feature map & 

Feature_Map_train=R_train{1,2}*G{2}*S'; 

% Mdl = generic_random_forests(full(Feature_Map_train), TRUE_Class_train,200,'classification'); % RF 
% Mdl= generic_random_forests(full(Feature_Map_train), TRUE_Class_train,40,'classification'); % RF
Mdl = fitcensemble(full(Feature_Map_train), TRUE_Class_train,'Method','Bag','FResample',.7,'NumLearningCycles',100);


%% predicting train data

[Pred_Class_train,Scores_pred_train] = predict(Mdl,full(Feature_Map_train)); 

%% predicting test samples from classifier RF

 Feature_Map_test= R_test{1,2}*full(G{2}*S'); 

[Pred_Class_test,Scores_pred_test] = predict(Mdl,full(Feature_Map_test)); 


TRUE_Class_test=row_test(2:end,end);
 
% acc=ismember(TRUE_Class_test,Pred_Class_test); %TP+TN
% ACCURACY=sum(acc==1)/(length(acc)); 


