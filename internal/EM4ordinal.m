function [alphaout,deltaout,Theta,S,M] = EM4ordinal(X,Y,alpha0,Delta0,tol)
% [alpha,delta,Theta,S,M] = EM4ordinal(X,Y,alpha0,Delta0,tol)
%
% This function estimates a full-rank basis matrix for a dimension reduction subspace 
% for ordinal data, using an approximate EM algorithm.
%
% INPUTS:
% X: predictor matrix, each row is an observation and each column a predictor variable.
% Y: response vector.
% alpha0: initial estimate of the basis matrix for dimension reduction
% Delta0: initial estimate of the residuals covariance matrix.
% tol: tolerance for convergence check.
%
% OUTPUTS:
% alpha: estimated basis matrix for the dimension reduction subspace.
% delta: estimated residual covariance matrix.
% Theta: estimated set of thresholds for the latent variables to generate
% the observed data.
% S: second moment estimates for the latent variables, given the observed
% data.
% M: first moment estimates for the latent variables, given the observed
% data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    error('Not enough input arguments for function EM4ordinal');
end
if nargin<5,
    tol = 1e-3;
end
[p,r] = size(alpha0);
n=length(Y);

if isinZ(Y),
    Fy = get_fyZ(Y);
else
    Fy = get_fy(Y,r);
end

% Initialization
if nargin<4,
  Delta0 = eye(p);
end
alpha = alpha0;
Delta = Delta0;

% ---- MAIN LOOP-----------
history = -1e6;
StopNotMet = 1;
iter = 0;

while (StopNotMet && iter<10)
    iter = iter+1;
    alphaout = alpha; 
    deltaout = Delta; 
        
    % update thresholds
    Theta = findThresolds(X,Fy,Delta,alpha);
    
    % -----------Start EM algorithm-------------
    % E-step
    [S,M] = computeEstep(X,Fy,Delta,alpha,Theta,'levina');
    
    % M-step
    Delta = (S - M'*Fy*inv(Fy'*Fy)*Fy'*M)/n;
    Delta = .5*(Delta+Delta');
    
    [YY,XX,opts] = prepare4glasso(X,Fy,Delta,M);
    hvec = inv(XX'*XX)*XX'*YY;
    alpha = vec2mtx(hvec,p,r);
    [alpha,Delta] = rescalePars_v2(alpha,Delta);
    Delta = .5*(Delta+Delta');
 
    % check convergence criteria... TBA
    [StopNotMet,history] = checkConvergenceOut(hvec,XX,YY,Delta,S,M,history,tol);
end
Theta = findThresolds(X,Fy,deltaout,alphaout);
    

