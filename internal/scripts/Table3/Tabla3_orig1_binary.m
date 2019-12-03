% This script estimates prediction error in classification for the EPH
% database, using the strategy named FULL in the manuscript.
% Compare results with the row corresponding to FULL in TABLE 3, for
% POVERTY response.
%
%------------------------------------------------------------------------

rng(160480);
nfold = 10;

% =======================================================================
data = load('traindata_binGBA.txt');
Y = data(:,1);
X = data(:,2:end);
nfold = 10;
[erORIG1_GBA,sdORIG1_GBA] = cvOrig1(Y,X,'logit',nfold);


data = load('traindata_binPampeana.txt');
Y = data(:,1);
X = data(:,2:end);
[erORIG1_Pampeana,sdORIG1_Pampeana] = cvOrig1(Y,X,'logit',nfold);


data = load('traindata_binNOA.txt');
Y = data(:,1);
X = data(:,2:end);
[erORIG1_NOA,sdORIG1_NOA] = cvOrig1(Y,X,'logit',nfold);


data = load('traindata_binNEA.txt');
Y = data(:,1);
X = data(:,2:end);
[erORIG1_NEA,sdORIG1_NEA] = cvOrig1(Y,X,'logit',nfold);


data = load('traindata_binPatagonia.txt');
Y = data(:,1);
X = data(:,2:end);
[erORIG1_Patagonia,sdORIG1_Patagonia] = cvOrig1(Y,X,'logit',nfold);


results_er = [erORIG1_GBA,erORIG1_Pampeana,erORIG1_NOA,erORIG1_NEA,erORIG1_Patagonia];
results_sd = [sdORIG1_GBA,sdORIG1_Pampeana,sdORIG1_NOA,sdORIG1_NEA,sdORIG1_Patagonia];
disp('Averaged Prediction Error over 10-fold CV is:')
disp(results_er)

