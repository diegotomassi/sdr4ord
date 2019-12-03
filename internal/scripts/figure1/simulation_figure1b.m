%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation_figura1b.m
%
% This simulation illustrates the performance of the ordinal-PFC
% method for nonnormal latent variables, as described in Section 5.2 of 
% L. Forzani, R. Garcia Arancibia, P. Llop and D. Tomassi, 
% "Sufficient dimension reduction for ordinal predictors" (submitted).
% Compare results with Figure 1-b) in the paper.
%
% Requires (included): 
% - LDR package, by RD Cook,L. Forzani and D. Tomassi:
% http://www.jstatsoft.org/article/view/v039i03

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ======= PERFORMANCE FOR LATENT VARIABLES DEVIATING FROM NORMALITY ======

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
%== reduction
alpha= [ones(1,p);sign(randn(1,p))]';
alpha=orth(alpha);
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

%%
%===== Initialization
angPFCnonnormal = ones(nreps,1); % angle between PFC estimate and true reduction
angORDnonnormal = ones(nreps,1); % angle between proposed estimate and true reduction

%===== main loop
parfor rep = 1:nreps,
    disp(['running repetition ',int2str(rep),' out of ',int2str(nreps)]);
    %=== response
    y = randn(N,1);
    %=== latent variables
    ruido = chi2rnd(5,[N,p]);
    scale = diag(1./sqrt(diag(get_cov(ruido))));
    ruido = (ruido - repmat(mean(ruido),N,1))*scale*sqrtm(delta);
    Fy = get_fy(y,r);
    Z = Fy*beta' + ruido;

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
    [~,alpha0] = ldr(y,X,'pfc','cont',r,'Fy',Fy);
    angPFCnonnormal(rep) = subspace(alpha0,alpha);
    
    %=== proposed method (ordinalPFC)
    [alpha0b,Delta0b] = get_initval(y,X,r);
    [alpha0,Delta0] = rescalePars_v2(alpha0b,Delta0b);
    alphaout = EM4ordinal_wlasso(X,y,alpha0,Delta0,0);
    angORDnonnormal(rep) = subspace(alphaout,alpha);

end

%% =========================================================================%
% FIGURE 1-B)
figure;
boxplot([angPFCnonnormal,angORDnonnormal]*180/pi,'labels',{'standard PFC', 'ordinal PFC'})
ylabel('angle');
ylim([10 90])
