% inference_perm.m
%
% This script reproduces results for inference about the dimension of the
% reduced subspace using sequential permutation testing. 
% Compare results with the corresponding column in TABLE 1 of the manuscript.

rng(1);
nreps=500;

p = 10;
d = 2;
r = 4;
n1 = 300;
n2 = 200;
dimset = 1:r;
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

dim_n1 = zeros(nreps,1);
dim_n2 = zeros(nreps,1);

N=n1;
parfor rep = 1:nreps,
    disp(['running repetition ', int2str(rep)])
    y = randn(N,1);
    ruido = randn(N,p);
    Fy = get_fy(y,r);
    Z = Fy*beta' + ruido;

    X = zeros(N,p);
    for j = 1:p,
        thetaj = Theta{j};
        Kj = length(thetaj) + 1;
        Zj = Z(:,j);
        Xj = Kj*ones(N,1);
        for k = (Kj-1):-1:1,
            idx = find(Zj < thetaj(k));
            Xj(idx) = k;
        end
        X(:,j) = Xj;
    end
    dim_n1(rep) = permEM4ordinal_vrr(y,X,r);
    dim_n2(rep) = permEM4ordinal_vrr(y(1:n2),X(1:n2,:),r);
end

p1_n1 = mean(dim_n1==1);
p2_n1 = mean(dim_n1==2);
p3_n1 = mean(dim_n1==3);
p4_n1 = mean(dim_n1==4);

p1_n2 = mean(dim_n2==1);
p2_n2 = mean(dim_n2==2);
p3_n2 = mean(dim_n2==3);
p4_n2 = mean(dim_n2==4);
