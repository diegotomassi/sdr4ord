function [YY,XX,Groups] = prepare4glasso(X,Fy,Delta,M)
XX = kron(sqrtm(Delta),Fy);
YY = M*invsqrtm(Delta); YY=YY(:);
p = size(Delta,2);
r = size(Fy,2);
Groups.ind = zeros(3,p);
for i=1:p,
    ini = ((i-1)*r+1);
    fin = i*r;
    Groups.ind(:,i) = [ini;fin;1];
end