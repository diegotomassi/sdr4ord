function [alphaout,deltaout,Theta,S,M,Aout,Bout] = EM4ordinal_vrr(X,Y,dim,alpha0,Delta0,A0,B0,lambda,tol)
% 
%
% [alpha,delta,Theta,S,M,A,B,chosen] = EM4ordinal_vrr(X,Y,dim,alpha0,Delta0,A0,B0,lambda,tol,thr)
%
% This function estimates a reduced-rank basis matrix for a dimension reduction subspace 
% for ordinal data, using an approximate EM algorithm. THIS VERSION DOES NOT INCLUDE VARIABLE SELECTION.
%
% INPUTS:
% X: predictor matrix, each row is an observation and each column a predictor variable.
% Y: response vector.
% dim: dimension of the dimension reduction subspace.
% alpha0: initial estimate of the basis matrix for dimension reduction
% Delta0: initial estimate of the residuals covariance matrix.
% A0: initial estimate for A (alpha in the final version of the draft)
% B0: initial estimate for B (varepsilon in the final version of the draft)
% tol: tolerance for convergence check.
%
% OUTPUTS:
% alpha: estimated basis matrix for the full-rank dimension reduction subspace.
% delta: estimated residual covariance matrix.
% Theta: estimated set of thresholds for the latent variables to generate
% the observed data.
% S: second moment estimates for the latent variables, given the observed
% data.
% M: first moment estimates for the latent variables, given the observed
% data.
% A: estimate for the reduced rank reduction (alpha in the final version of the draft).
% B: estimate for varepsilon in the draft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 7
    error('Not enough input arguments for function EM4ordinal');
end

if nargin<9,
    tol = 1e-3;
end

if nargin<8,
    lambda = .001;
end

[p,r] = size(alpha0);
n=length(Y);

if isinZ(Y),
    Fy = get_fyZ(Y);
else
    Fy = get_fy(Y,r);
end

% Initialization
if nargin<5,
  Delta0 = eye(p);
end
alpha = alpha0;
A = A0;
B = B0;
Delta = Delta0;

% ---- FIRST STEP NOW INCLUDED IN MAIN LOOP-----------
history = -1e6;
StopNotMet = 1;
iter = 0;
while (StopNotMet && iter<5)
    iter = iter+1;
    alphaout = alpha; 
    deltaout = Delta; 
    Aout = A;
    Bout = B;
        
    % update thresholds
    Theta = findThresolds(X,Fy,Delta,alpha);
    
    % -----------Start EM algorithm-------------
    % E-step
    [S,M] = computeEstep(X,Fy,Delta,alpha,Theta,'levina');
    
    % M-step

    % update of A (alpha in the draft) 
    %----- option 1: without LASSO
    iShalf = invsqrtm(S/n);
    Sfit = M'*Fy*inv(Fy'*Fy)*Fy'*M/n;
    [v,~] = firsteigs(iShalf*Sfit*iShalf,dim);
    A = orth(iShalf*v);
    
    % update of DELTA
%     a0 = orthcomp(A);
    invS = inv(S/n);
%     Aaux = a0*a0' + a0*inv(a0'*invS*a0)*a0'*invS*A*A';
    Sres = S/n - Sfit;
    Delta = Sres ;%+ Aaux*Sfit + Sfit*Aaux' - Aaux*Sfit*Aaux';
    Delta = .5*(Delta+Delta');
%     isposdef(Delta)    
%     diag(Delta)
    
%     % rescaling ===
     [A,Delta] = rescalePars_v2(A,Delta);
     Delta = .5*(Delta+Delta');

    % update of B (psi in the draft)
    B = inv(A'*Delta*A)*A'*M'*Fy*inv(Fy'*Fy);

    alpha = A*B;
    
    % check convergence criteria... TBA
    [YY,XX,Groups] = prepare4glasso(X,Fy,Delta,M);
    hvec = alpha'; hvec = hvec(:);
    [StopNotMet,history] = checkConvergenceOut(hvec,XX,YY,Delta,S,M,history,tol);
%     pause
end
