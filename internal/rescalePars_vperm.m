function [alphaout,Deltaout,sig] = rescalePars_vperm(alpha,Delta)
% [alphaout,Deltaout] = rescalePars_v2(alpha,Delta)
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

prod = Delta*alpha;
[Delta,sig] = corrcov(Delta);
Deltaout = .5*(Delta + Delta');
alphaout = invsqrtm(Deltaout)*prod;


