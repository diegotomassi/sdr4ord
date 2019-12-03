% This script computes predictions scores for method referred to as 
% PFCord in TABLE 3 of the manuscript 
% L. Forzani, R. Garcia Arancibia, P. Llop and D. Tomassi, 
% "Sufficient dimension reduction for ordinal predictors" (submitted). 
% 
% Only continuous response is considered in the script.
% ========================================================================
rng(160480);
GBA = load('traindata_GBA.txt');
[N,p] = size(GBA);
dim = 1;
r = 2;
lambda = 500;
Y = GBA(:,1)/1000; X = GBA(:,2:p);
[erPFCord_GBA,sdPFCord_GBA] = cvPFCord(Y,X,dim,r,lambda,'linear');
    

GBA = load('traindata_Pampeana.txt');
[N,p] = size(GBA);
dim = 1;
r = 3;
lambda = 500;
Y = GBA(:,1)/1000; X = GBA(:,2:p);
[erPFCord_Pampeana,sdPFCord_Pampeana] = cvPFCord(Y,X,dim,r,lambda,'linear');


GBA = load('traindata_NOA.txt');
[N,p] = size(GBA);
dim = 1;
r = 3;
lambda = 500;
Y = GBA(:,1)/1000; X = GBA(:,2:p);
[erPFCord_NOA,sdPFCord_NOA] = cvPFCord(Y,X,dim,r,lambda,'linear');


GBA = load('traindata_NEA.txt');
[N,p] = size(GBA);
dim = 1;
r = 3;
lambda = 500;
Y = GBA(:,1)/1000; X = GBA(:,2:p);
[erPFCord_NEA,sdPFCord_NEA] = cvPFCord(Y,X,dim,r,lambda,'linear');


GBA = load('traindata_Patagonia.txt');
[N,p] = size(GBA);
dim = 1;
r = 3;
lambda = 500;
Y = GBA(:,1)/1000; X = GBA(:,2:p);
[erPFCord_Patagonia,sdPFCord_Patagonia] = cvPFCord(Y,X,dim,r,lambda,'linear');

results_er = [erPFCord_GBA,erPFCord_Pampeana,erPFCord_NOA,erPFCord_NEA,erPFCord_Patagonia];
results_sd = [sdPFCord_GBA,sdPFCord_Pampeana,sdPFCord_NOA,sdPFCord_NEA,sdPFCord_Patagonia];
disp('Averaged MSE over 10-fold CV is:')
disp(results_er)



