function [fn,Wn] = pfc4perm(Y,X,u,r)
%
% This function implements the Principal Fitted Components (EPFC) 
% model for Dimension Reduction in Regression (Cook & Forzani 2009).
% USAGE:
% - outputs:
%     fn: value of the loss function at the optimal point;
%     Wn: generating vectors for the central subspace;
%  - inputs:
%     Y: Response vector. 
%     X: Data matrix. Each row is an observation. It is assumed
%        that rows relate with the corresponding rows in Y, so that Y(k) is 
%	    the response due to X(k,:). 
%     u: Dimension of the sufficient subspace. It must be a 
%        natural greater than 1 and smaller than the number of columns in X.
%     r: order of polynomial basis for th einverse regression of X on Y.
% --------------------------- REFERENCES -----------------------------------
% Cook, R. D. and Forzani, L. (2009). Principal fitted components in regression.
% Statistcal Science 23 (4), 485-501.
%
% =======================================================

[n,p] = size(X);
epsi = 1e-15; epsimtx = epsi*eye(p);
Fy = get_fy(Y,r);
SIGMAfit = get_fitted_cov(Y,X,Fy);
SIGAMfit= SIGMAfit+epsimtx;
SIGMA = get_cov(X)+epsimtx;

%--- optimization .........................................................
aux = invsqrtm(SIGMA)*SIGMAfit*invsqrtm(SIGMA);
[Wn,vals] = firsteigs(aux,u);
Wn = orth(invsqrtm(SIGMA)*Wn);
fn = -sum(log(vals));
