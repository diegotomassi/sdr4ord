function [StopNotMet,history] = checkConvergenceOut(hvec,XX,YY,Delta,S,M,history,tol)
if nargin < 8,
    tol = 1e-3;
end
StopNotMet = 1;
funval = ordLogLik(YY,XX,hvec,Delta,S,M);
if ((funval < history)||( abs(history-funval)/abs(funval) < tol)),
    StopNotMet = 0;
end
history = [funval];


