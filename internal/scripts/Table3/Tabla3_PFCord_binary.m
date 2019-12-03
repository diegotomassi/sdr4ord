% This script computes predictions scores for method referred to as 
% reg-PFCord in TABLE 3 of the manuscript 
% L. Forzani, R. Garcia Arancibia, P. Llop and D. Tomassi, 
% "Sufficient dimension reduction for ordinal predictors" (submitted). 
% 
% Only binary response is considered in the script.
% ========================================================================
rng(160480);
GBA = load('traindata_binGBA.txt');
[N,p] = size(GBA);
dim = 1;
r = 3;
lambda = 500;
Y = GBA(:,1); X = GBA(:,2:p);
[erPFCord_GBA,sdPFCord_GBA] = cvPFCord(Y,X,dim,r,lambda,'logit');
    
%%
Pampeana = load('traindata_binPampeana.txt');
[N,p] = size(Pampeana);
dim = 1;
r = 3;
lambda = 500;
Y = Pampeana(:,1); X = Pampeana(:,2:p);
[erPFCord_Pampeana,sdPFCord_Pampeana] = cvPFCord(Y,X,dim,r,lambda,'logit');

%%
NOA = load('traindata_binNOA.txt');
[N,p] = size(NOA);
dim = 1;
r = 3;
lambda = 500;
Y = NOA(:,1); X = NOA(:,2:p);
[erPFCord_NOA,sdPFCord_NOA] = cvPFCord(Y,X,dim,r,lambda,'logit');

%%
NEA = load('traindata_binNEA.txt');
[N,p] = size(NEA);
dim = 1;
r = 3;
lambda = 500;
Y = NEA(:,1); X = NEA(:,2:p);
[erPFCord_NEA,sdPFCord_NEA] = cvPFCord(Y,X,dim,r,lambda,'logit');

%%
Patagonia = load('traindata_binPatagonia.txt');
[N,p] = size(Patagonia);
dim = 1;
r = 3;
lambda = 500;
Y = Patagonia(:,1); X = Patagonia(:,2:p);
[erPFCord_Patagonia,sdPFCord_Patagonia] = cvPFCord(Y,X,dim,r,lambda,'logit');

results_er = [erPFCord_GBA,erPFCord_Pampeana,erPFCord_NOA,erPFCord_NEA,erPFCord_Patagonia];
results_sd = [sdPFCord_GBA,sdPFCord_Pampeana,sdPFCord_NOA,sdPFCord_NEA,sdPFCord_Patagonia];

disp('Averaged MSE over 10-fold CV is:')
disp(results_er)


