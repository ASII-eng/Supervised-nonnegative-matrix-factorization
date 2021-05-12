function [VAR_VAR]=Variant_Varient_Network(num,row,Score_names)

 % this function construct variant-variant network 
%------------------------------------------------------------ 

Label=string(zeros(length(num),1)); % pre allocate
%% set the Variant labels
for i=1:length(num)
Label(i)=strcat(string(num(i,1)),row{i+1,2},row{i+1,3});
end

%% set the NAN values to zero for Scores that dont exist in file(missing data)
nans= isnan(num);
num(nans)=0;

%% Selecting the desired SCORES 
% applied a linear transformation to adjust scores according to two criteria: 1) a score should be in
%the interval of [0, 1] and 2) the larger a score, the stronger the evidence of functionally damaging. 
% For SIFT and LRT, we convert them to 1-SIFT and 1-LRT. 
% Then, we transform All scores with formula (f- fmin)/(fmax-fmin), 

num(:,4)=1-num(:,4); %SIFT=1-SIFT
num(:,8)=1-num(:,7); %LRT=1-LRT
%
labels=row(1,:);
[~,col_num]=ismember(Score_names,labels);
Score_SELECTED= num(:,col_num); %select the Scores
FMAX=[1,1,1,1,1,6.49,10.64,14,1,3,1,1,18.301497,35,1,1,1 ,1,1,6.17,10.003,1.199,1,1,37.9718,1]; % eigen & phyloP20way_mammalian%1*29
FMIN=[0,0,0,0,0,-5.17,-16.13,-14,0,-2,0,0,-6.458163,0,0,0,0,0,0,-12.3,-20,-13.282,0,0,0,0]; %1*26
FMAX=FMAX(col_num-3);
FMIN=FMIN(col_num-3);
Score_SELECTED2=(Score_SELECTED-FMIN)./(FMAX-FMIN);
 VAR=corrcoef(Score_SELECTED2');

VAR=(VAR+1)./2; %to transform [-1,1] values to [0,1]


%% costruct variant-variant network
VAR_VAR={cellstr(Label),sparse(VAR)};

end