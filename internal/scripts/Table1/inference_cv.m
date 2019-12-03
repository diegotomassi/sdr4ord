% inference_cv.m
%
% This script reproduces results for inference about the dimension of the
% reduced subspace using cross-validation. Compare results with the
% corresponding column in TABLE 1 of the manuscript.

rng(1);
nreps=500;

p = 10;
d = 2;
r = 4;
n1 = 300;
n2 = 200;
dimset = 1:r;
nfolds=10;
delta = eye(p);
alpha= [ones(1,p);sign(randn(1,p))]';
epsilon = randn(d,r);
alpha=orth(alpha);
%== reduction
beta = delta*alpha*epsilon;

%== Define thresholds on latent variables
Theta = cell(p,1);
Theta{1}=[-1,0,1,2];
Theta{2}=[-1,0,1];
Theta{3}=[0,1];
Theta{4}=[-1,1];
Theta{5}=[-1,0];
Theta{6}=[-1,0,1,2];
Theta{7}=[-1,0,1];
Theta{8}=[0,1];
Theta{9}=[-1,1];
Theta{10}=[-1,0];
Theta{11}=[-1,0,1,2];
Theta{12}=[-1,0,1];
Theta{13}=[0,1];
Theta{14}=[-1,1];
Theta{15}=[-1,0];
Theta{16}=[-1,0,1,2];
Theta{17}=[-1,0,1];
Theta{18}=[0,1];
Theta{19}=[-1,1];
Theta{20}=[-1,0];

%%
scores_n1 = zeros(size(dimset));
scores_n2 = zeros(size(dimset));

for rep = 1:nreps,
    disp(['running repetition ',int2str(rep),' out of ',int2str(nreps)]);
    %=== response
    y = randn(n1,1);
    %=== latent variables
    Fy = get_fy(y,r);
    Z = mvnrnd(Fy*beta',delta);
    %=== observed variables
    X = zeros(n1,p);
    for j = 1:p,
        thetaj = Theta{j};
        Kj = length(thetaj) + 1;
        Zj = Z(:,j);
        Xj = Kj*ones(n1,1);
        for k = (Kj-1):-1:1,
            idx = find(Zj < thetaj(k));
            Xj(idx) = k;
        end
        X(:,j) = Xj;
    end
    Xn2 = X(1:n2,:); yn2=y(1:n2);
    
    %=== START CV procedurenor n=n1
    blksize = floor(n1/nfolds);
    ercv =zeros(nfolds,length(dimset));
    for fold=1:nfolds,
        if fold < nfolds,
            idx = ((fold-1)*blksize+1):(fold*blksize);
        else
            idx = ((fold-1)*blksize+1):n1;
        end
        xtr = X; xtr(idx,:)=[]; ytr = y; ytr(idx) = [];
        xts = X(idx,:);         yts = y(idx);
        
        
    %=== proposed method (ordinalPFC)
    ytshat = zeros(size(yts));
    for kkk=dimset,
        [alpha0b,Delta0b,A0,B0] = get_initval_v2(ytr,xtr,kkk,r);
        [A0,Delta0] = rescalePars_v2(A0,Delta0b);
        [~,~,~,~,~,alphaout,~] = EM4ordinal_vrr(xtr,ytr,kkk,alpha0b,Delta0,A0,B0);
        ID = knnsearch(xtr*alphaout,xts*alphaout,'K',6);
        
        for jj=1:length(yts),
            ytshat(jj)=mean(ytr(ID(jj,:)));
        end
        ercv(fold,kkk)=mean(dot(ytshat-yts,ytshat-yts));
    end
    end
    [~,kmin]=min(mean(ercv));
    scores_n1(kmin)=scores_n1(kmin)+1;
    
    
    %=== START CV procedure for n=n2
    blksize = floor(n2/nfolds);
    ercvn2 =zeros(nfolds,length(dimset));
    for fold=1:nfolds,
        if fold < nfolds,
            idx = ((fold-1)*blksize+1):(fold*blksize);
        else
            idx = ((fold-1)*blksize+1):n2;
        end
        xtr = Xn2; xtr(idx,:)=[]; ytr = yn2; ytr(idx) = [];
        xts = Xn2(idx,:);         yts = yn2(idx);
        
        
    %=== proposed method (ordinalPFC)
    ytshat = zeros(size(yts));
    for kkk=dimset,
        [alpha0b,Delta0b,A0,B0] = get_initval_v2(ytr,xtr,kkk,r);
        [A0,Delta0] = rescalePars_v2(A0,Delta0b);
        [~,~,~,~,~,alphaout,~] = EM4ordinal_vrr(xtr,ytr,kkk,alpha0b,Delta0,A0,B0);
        ID = knnsearch(xtr*alphaout,xts*alphaout,'K',6);
        
        for jj=1:length(yts),
            ytshat(jj)=mean(ytr(ID(jj,:)));
        end
        ercvn2(fold,kkk)=mean(dot(ytshat-yts,ytshat-yts));
    end
    end
    [~,kmin]=min(mean(ercvn2));
    scores_n2(kmin)=scores_n2(kmin)+1;
end
%%
