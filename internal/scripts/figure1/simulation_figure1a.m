%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This simulation illustrates the performance of the ordinal-PFC
% method for normal latent variables, as described in Section 5.2 of 
% L. Forzani, R. Garcia Arancibia, P. Llop and D. Tomassi, 
% "Sufficient dimension reduction for ordinal predictors" (submitted).
% Compare results with Figure 1-a) in the paper.
%
% Requires: 
% - LDR package, by RD Cook,L. Forzani and D. Tomassi:
% http://www.jstatsoft.org/article/view/v039i03

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Starting simulation. This can take some time...']);

%% ============== PERFORMANCE FOR NORMAL LATENT VARIABLES ===============

%===== Simulation parameters
N = 500;
nreps = 100;
rng(160480);

%===== Model parameters for latent-variables
p = 20;
r = 2;
dim = 2;
%== covariance matrix
aux = rand(p); 
delta = 4*eye(p) + 0.3*aux*aux';
alpha= [ones(1,p);sign(randn(1,p))]';
alpha=orth(alpha);
%== reduction
beta = delta*alpha;

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

%===== Initialization
angPFC = ones(nreps,1); % angle between PFC estimate and true reduction
angORD = ones(nreps,1); % angle between proposed estimate and true reduction
%%
parfor rep = 1:nreps,
    disp(['running repetition ',int2str(rep),' out of ',int2str(nreps)]);
    %=== response
    y = randn(N,1);
    %=== latent variables
    Fy = get_fy(y,r);
    Z = mvnrnd(Fy*beta',delta);
    %=== observed variables
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

    %=== standard PFC estimate
    Fy = get_fy(y,r);
    [~,alpha0] = ldr(y,X,'pfc','cont',r,'Fy',Fy);
    angPFC(rep) = subspace(alpha0,alpha);
    
    %=== proposed method (ordinalPFC)
    [alpha0b,Delta0b] = get_initval(y,X,r);
    [alpha0,Delta0] = rescalePars_v2(alpha0b,Delta0b);
    alphaout = EM4ordinal_wlasso(X,y,alpha0,Delta0,0);
    angORD(rep) = subspace(alphaout,alpha);

end

%% FIGURE 1-A)
figure;
boxplot([angPFC,angORD]*180/pi,'labels',{'standard PFC', 'ordinal PFC'})
ylabel('angle');
ylim([10 90])






