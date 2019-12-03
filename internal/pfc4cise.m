function [fn,Wn,st] = pfc4cise(Y,X,u,r,lmax)
% [Wn,fn,fp] = pfc(Y,X,u)
%
% This function implements the Principal Fitted Components (EPFC) 
% model for Dimension Reduction in Regression (Cook & Forzani 2009).
% USAGE:
% - outputs:
%     Wn: generating vectors for the central subspace;
%     fn: value of the loss function at the optimal point;
%     fp: value of the loss function for the original predictors;
%  - inputs:
%     Y: Response vector. 
%     X: Data matrix. Each row is an observation. It is assumed
%        that rows relate with the corresponding rows in Y, so that Y(k) is 
%	    the response due to X(k,:). 
%     u: Dimension of the sufficient subspace. It must be a 
%        natural greater than 1 and smaller than the number of columns in X.
% --------------------------- REFERENCES -----------------------------------
% Cook, R. D. and Forzani, L. (2009). Principal fitted components in regression.
% Statistcal Science 23 (4), 485-501.
%
% =======================================================
if nargin < 5,
    lmax=0;
end
[n,p] = size(X);
[X,idx]=unique(X,'rows');
Y = Y(idx);
epsi = 1e-16; epsimtx = epsi*eye(p);
Fy = get_fy(Y,r);
SIGMAfit = get_fitted_cov(Y,X,Fy);
% SIGAMfit= SIGMAfit+epsimtx;
SIGMA = get_cov(X);
SIGMA = SIGMA+epsimtx;
% SIGMAres = SIGMA - SIGMAfit;

%--- optimization .........................................................
aux = invsqrtm(SIGMA)*SIGMAfit*invsqrtm(SIGMA);
[Wn,vals] = firsteigs(aux,u);
Wn = orth(invsqrtm(SIGMA)*Wn);
parameters.S=SIGMA;
parameters.Sfit=SIGMAfit;
[fn,beta,st] = mycise(Y,X,u,lmax,parameters);
Wn=orth(beta);