% This script estimates prediction error in regression for the EPH
% database, using the strategy named FULL in the manuscript.
% Compare results with the row corresponding to FULL in TABLE 3, for
% INCOME PER CAPITA response.
%
%------------------------------------------------------------------------

rng(160480);
data = load('traindata_GBA.txt');
Y = data(:,1)/1000;
X = data(:,2:end);
nfold = 10;
[erORIG1_GBA,sdORIG1_GBA] = cvOrig1(Y,X,'linear',nfold);


data = load('traindata_Pampeana.txt');
Y = data(:,1)/1000;
X = data(:,2:end);
[erORIG1_Pampeana,sdORIG1_Pampeana] = cvOrig1(Y,X,'linear',nfold);


data = load('traindata_NOA.txt');
Y = data(:,1)/1000;
X = data(:,2:end);
[erORIG1_NOA,sdORIG1_NOA] = cvOrig1(Y,X,'linear',nfold);


data = load('traindata_NEA.txt');
Y = data(:,1)/1000;
X = data(:,2:end);
[erORIG1_NEA,sdORIG1_NEA] = cvOrig1(Y,X,'linear',nfold);


data = load('traindata_Patagonia.txt');
Y = data(:,1)/1000;
X = data(:,2:end);
[erORIG1_Patagonia,sdORIG1_Patagonia] = cvOrig1(Y,X,'linear',nfold);

results_er = [erORIG1_GBA,erORIG1_Pampeana,erORIG1_NOA,erORIG1_NEA,erORIG1_Patagonia];
disp('Averaged MSE over 10-fold CV is:')
disp(results_er)
results_sd = [sdORIG1_GBA,sdORIG1_Pampeana,sdORIG1_NOA,sdORIG1_NEA,sdORIG1_Patagonia];
disp('Standard deviation for the 10-fold CV procedure is:')
disp(results_sd)



%%
results_er = [erORIG1_GBA,erORIG1_Pampeana,erORIG1_NOA,erORIG1_NEA,erORIG1_Patagonia];
disp('Averaged Prediction Error over 10-fold CV is:')
disp(results_er)

