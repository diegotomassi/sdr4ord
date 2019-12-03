function fval = ordLogLik(Y,X,hvec,Delta,S,M)
% ordLogLik.m
% 
% This function computes the log-likelihhod of the estimated model 
% for the observed data.

[n,p] = size(M);
invDelta = inv(Delta);
apperror = Y - X*hvec;
%traza = trace(Delta\S) + apperror'*apperror - trace((Delta\M')*M);
traza = trace(invDelta*S) + apperror'*apperror - trace(invDelta*M'*M);
fval = .5*n*logdet(invDelta) - .5*traza;
