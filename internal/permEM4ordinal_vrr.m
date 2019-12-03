function [d] = permEM4ordinal(Y,X,r,nsample)
%
% function [Wmin,d,f] = permLAD(Y,X,morph,parameters);
% 
% This function estimates the dimension of the reduced subspace that best
% describes the data under the LAD model using a permutation test.  
% 
% USAGE:
%  - outputs:
%    - Wmin: generating vectors for the central subespace of estimated
%    dimension.
%    - d: estimated dimension under LRT.
%    - f: value of the optimized function for dimension d. (perhaps this is useless)
%  - inputs: 
%    - Y: response vector;
%    - X: matrix of predictors;
% =========================================================================
if nargin < 4,
    nsample=999;
end
alphalevel=0.01;
%----checking type of response and slicing if needed.......................
if isinZ(Y),
    morph = 'disc';
     Fy = get_fyZ(Y);
     Delta = get_deltaZ(X,Y);
else 
    morph = 'cont';
    Fy = get_fy(Y,r);
    Delta = get_delta(X,Y,Fy);
end

[n,p] = size(X);

% likelihood con todos los predictores
    % estimador inicial
    [alpha0,Delta0,A0,B0] = get_initval_v2(Y,X,r,r);
    [A0,Delta0] = rescalePars_v2(A0,Delta0);
    %
    % Reduccion ordinal
    [alpha,Delta,~,S,M,A,B,f] = EM4ordinal_vperm(X,Y,r,alpha0,Delta0,A0,B0);

fpo = f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for u=1:p-1,
    % estimador inicial
    [alpha0,Delta0,A0,B0] = get_initval_v2(Y,X,u,r);
    [A0,Delta0] = rescalePars_vperm(A0,Delta0);
    %
    % Reduccion ordinal
    [alpha,Delta,~,S,M,A,B,f] = EM4ordinal_vperm(X,Y,u,alpha0,Delta0,A0,B0);
    % likelihood ordinal
    fno = f;
    
    % estadistico observado
    T = 2*(fno-fpo);
    
    %%%%%%%%% PERMUTACIONES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPLEMENTO ORTOGONAL
    Wu = alpha;
    [Q R] = qr(Wu); 
    Wu0 = Q(:,(u+1):p);
    Xnew = X*[Wu Wu0];
    Xrsp = zeros(n,p);
    fns = 1:nsample;
    fps = 1:nsample;
    for i=1:nsample,
        % fijamos las primeras u columnas y permutamos el resto
        Yperm = randperm(n); % ver implementaci??n alternativa en denboot.m
        Xrsp(:,1:u) = Xnew(:,1:u);
        Xrsp(:,(u+1):p) = Xnew(Yperm(:),(u+1):p);

        % logLik sin reduccion para datos permutados
        % estimador inicial
        fps(i) = pfc4perm(Y,Xrsp,r,r);
        fns(i) = pfc4perm(Y,Xrsp,u,r);
    end
    fnmfp = 2*(fns - fps);
    prop = sum(fnmfp > T);
    aux = prop/nsample;
    if (aux > (alphalevel))||(u==min(r,p)),
        d = u;
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
