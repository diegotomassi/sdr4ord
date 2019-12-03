function [alpha,Delta,A,B] = get_initval(Y,X,dim,r)

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
% alpha = inv(Delta)*X'*Fy*aux;
% nueva version alpha=AB
    ideltahalf = invsqrtm(Delta);
    [v,~] = firsteigs(ideltahalf*cov(X)*ideltahalf,dim);
    A = ideltahalf*v;
    B = inv(A'*Delta*A)*A'*X'*Fy*aux;
    alpha = A*B;


    