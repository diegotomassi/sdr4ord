function [StopNotMet,history] = checkConvergenceIn(S,Sold,history,tol)
if nargin<4,
    tol = 1e-3;
end
StopNotMet = 1;
fro = norm(S-Sold,'fro');%/norm(S,'fro');
if sum(abs(history-fro)/fro)<tol,
    StopNotMet = 0;
end
history = [fro history]; history(end) = [];
