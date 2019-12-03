%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation_figure2.m
%
% This simulation illustrates the performance of the ordinal-PFC
% method including regularization, as described in Section 5.4 of 
% L. Forzani, R. Garcia Arancibia, P. Llop and D. Tomassi, 
% "Sufficient dimension reduction for ordinal predictors" (submitted).
% Compare results with Figure 2 in the paper.
%
% Requires (included): 
% - LDR package, by RD Cook,L. Forzani and D. Tomassi:
% http://www.jstatsoft.org/article/view/v039i03
% - SLEP package by J. Liu, S. Ji, and J. Ye. 
% http://www.public.asu.edu/~jye02/ Software/SLEP.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rng(2015);

N = 500;
nreps=100;

p = 10;
r = 2;
aux = rand(p); 
delta = 4*eye(p) + 0.1*aux*aux';
alpha= [1 -1;
        1  1;
        1  1;
        1 -1;
        zeros(p-4,2)];
alpha=orth(alpha);

beta = delta*alpha;

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

%%
lambda = 100; % regularization parameter (not optimized)
angPFC = ones(nreps,1);
angORD = ones(nreps,1);
angORDlasso = ones(nreps,1);

betas = cell(nreps,1);
parfor rep = 1:nreps,
    disp(['running repetition ',int2str(rep),' out of ',int2str(nreps)]);
    y = randn(N,1);
    ruido = randn(N,p);
    Fy = get_fy(y,r);
    Z = Fy*beta' + ruido*sqrtm(delta);
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

    %---- PFC comun
    Fy = get_fy(y,r);
    [~,alpha0] = ldr(y,X,'pfc','cont',r,'Fy',Fy);
    angPFC(rep) = subspace(alpha0,alpha);
    
    %--- PFC ordinal
    [alpha0b,Delta0b] = get_initval(y,X,r);
    [alpha0,Delta0] = rescalePars_v2(alpha0b,Delta0b);
    alphaout = EM4ordinal_wlasso(X,y,alpha0,Delta0,lambda);
    angORDlasso(rep) = subspace(alphaout,alpha);
    
    [alphaout] = EM4ordinal(X,y,alpha0,Delta0);
    angORD(rep) = subspace(alphaout,alpha);
end
%%
resultados = [angPFC,angORD,angORDlasso]*180/pi;



%%
figure;
boxplot(resultados,'labels',{'standard PFC', 'ordinal PFC', 'ordinal+lasso'})
ylabel('angle')