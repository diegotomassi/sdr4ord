function [alpha,Delta] = get_initval(Y,X,r)
% [alpha,Delta] = get_initval(Y,X,r)
% 
% This function computes an initial value for the EM algorithm.
% 
% - INPUTS:
% Y: response (column) vector;
% X: predictor matrix, each row is an observation.
% r: order for the polynomial basis for the inverse regresion of Y on X.
% 
% - OUTPUTS:
% alpha: initial estimate of the basis matrix for the reduction.
% Delta: initial estimate of the residual covariance matrix.

[n,p] = size(X);
X = X - repmat(mean(X),n,1);
if isinZ(Y),
    Fy = get_fyZ(Y);
else
    Fy = get_fy(Y,r);
end
aux = inv(Fy'*Fy);
Qf = eye(n)-Fy*aux*Fy';
Delta = (X'*Qf*X)/n;
alpha = inv(Delta)*X'*Fy*aux;
