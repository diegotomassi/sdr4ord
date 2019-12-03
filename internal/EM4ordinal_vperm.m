function [alphaout,deltaout,Theta,S,M,Aout,Bout,f] = EM4ordinal(X,Y,dim,alpha0,Delta0,A0,B0,lambda,tol)
% THIS VERSION TRIES TO SPEED-UP THE PERMUTATION TESTING PROCEDURE
%
% [alpha,delta,Theta,S,M] = EM4ordinal(X,Y,alpha0,Delta0,tol)
%
% This function estimates a reduced-rank (d<r) basis matrix for a dimension reduction subspace 
% for ordinal data, using an approximate EM algorithm.
%
% INPUTS:
% X: predictor matrix, each row is an observation and each column a predictor variable.
% Y: response vector.
% dim: dimesnion of the characteristic subspace.
% alpha0: initial estimate of (alpha*varepsilom)
% Delta0: initial estimate of the residuals covariance matrix.
% A0: initial estimate of the basis matrix for the reduction subspace (alpha).
% B0: initial estimate of epsilon in model (4)
% lambda: (optional) regularization parameter
% tol: (optional) tolerance for convergence check.
%
% OUTPUTS:
% alpha: estimate of (alpha*epsilon) in model (4).
% delta: estimated residual covariance matrix.
% Theta: estimated set of thresholds for the latent variables to generate
% the observed data.
% S: second moment estimates for the latent variables, given the observed
% data.
% M: first moment estimates for the latent variables, given the observed
% data.
% Aout: estimated basis matrix for the dimension reduction subspace (alpha)
% Bout: estimate of epsilon in model (4)
%
%
% For details, see 
% L. Forzani, R. Garcia Arancibia, P. Llop and D. Tomassi, 
% "Sufficient dimension reduction for ordinal predictors" (submitted).
% Compare results with Figure 1-a) in the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ========================================================================

if nargin < 7
    error('Not enough input arguments for function EM4ordinal');
end

if nargin<9,
    tol = 1e-3;
end

if nargin<8,
    lambda = 0;
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
    % ----------------E-step
    [S,M] = computeEstep(X,Fy,Delta,alpha,Theta,'levina');
    
    % ----------------M-step

    % update of A (alpha in the draft) 
    iShalf = invsqrtm(S/n);
    Sfit = M'*Fy*inv(Fy'*Fy)*Fy'*M/n;
    [v,vals] = firsteigs(iShalf*Sfit*iShalf,dim);
    A = orth(iShalf*v);
    f = -sum(log(vals));
    
    % update of DELTA
    invS = inv(S/n);
    Sres = S/n - Sfit;
    Delta = Sres; % + Aaux*Sfit + Sfit*Aaux' - Aaux*Sfit*Aaux';
    Delta = .5*(Delta+Delta');
    
     % rescaling ===
      [A,Delta] = rescalePars_vperm(A,Delta);
      Delta = .5*(Delta+Delta');

    % update of B (psi in the draft)
    B = inv(A'*Delta*A)*A'*M'*Fy*inv(Fy'*Fy);

    alpha = A*B;
    
    % check convergence criteria
    [YY,XX,Groups] = prepare4glasso(X,Fy,Delta,M);
    hvec = alpha'; hvec = hvec(:);
    [StopNotMet,history] = checkConvergenceOut(hvec,XX,YY,Delta,S,M,history,tol);
end
