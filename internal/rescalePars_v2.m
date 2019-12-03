function [alphaout,Deltaout,sig] = rescalePars_v2(alpha,Delta)
% [alphaout,Deltaout,sig] = rescalePars_v2(alpha,Delta)
% 
% This function rescales alpha and Delta in order for Delta to have unit 
% diagonal elements (Delta)_ii = 1.
% 
% - INPUTS: 
% alpha: basis matrix for the dimesnion reduction subspace
% Delta: residual covariance matrix.
% 
% - OUTPUTS:
% alphaout: rescaled basis matrix for the dimension reduction subspace.
% Deltaout: rescaled residual covariance matrix.
% sig: vector of standard deviations of predictors.

[Deltaout,sig] = corrcov(Delta);
alphaout=diag(sig)*alpha;

