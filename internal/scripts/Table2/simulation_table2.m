%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation52a.m
%
% This simulation illustrates the performance of the ordinal-PFC
% method for normal latent variables, as described in Section 5.4 of 
% L. Forzani, R. Garcia Arancibia, P. Llop and D. Tomassi, 
% "Sufficient dimension reduction for ordinal predictors" (submitted).
% Compare results with TABLE 2 in the paper.
%
% Requires: 
% - (PROVIDED) LDR package, by RD Cook,L. Forzani and D. Tomassi:
% http://www.jstatsoft.org/article/view/v039i03

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Starting simulation. This can take some time...']);

%% ============== PERFORMANCE FOR NORMAL LATENT VARIABLES ===============
dcorcase='high';
%===== Simulation parameters
N = 200;
nreps = 100;
rng(160480);

%===== Model parameters for latent-variables
p = 20;
r = 2;
dim = 2;
%== covariance matrix
aux = rand(p); 
% DCOR approx 0.18
switch dcorcase,
    case 'independent'
        delta = eye(p);
    case 'low'
    % DCOR approx 0.30
    delta = 4*eye(p) + 0.2*aux*aux';
    case 'moderate'
    % DCOR approx 0.50
    delta = 4*eye(p) + 0.3*aux*aux';
    case 'high'
    % DCOR approx 0.7
    delta = 4*eye(p) + 0.5*aux*aux';
end

alpha= [ones(1,p);sign(randn(1,p))]';
alpha(5:p,:)=0;
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
selection_pfc = zeros(p,nreps);
selection = zeros(p,nreps);

%%
parfor rep = 1:nreps,
    disp(['running repetition ',int2str(rep),' out of ',int2str(nreps)]);
    %=== response
    y = randn(N,1);
    %=== latent variables
    Fy = get_fy(y,r);
    Z = mvnrnd(Fy*beta',delta);
    corre(rep)=distcorr(Z(:,1:4),Z(:,5:p))
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
    [fn,Wn,st] = pfc4cise(y,X,dim,r,100);
    selection_pfc(:,rep)=st;    
    %=== proposed method (ordinalPFC)
    [alpha0b,Delta0b,A0,B0] = get_initval_v2(y,X,dim,r);
     [A0,Delta0] = rescalePars_v2(A0,Delta0b);
    [~,~,~,~,~,alphaout,~,st] = EM4ordinal_wcise(X,y,dim,alpha0b,Delta0,A0,B0,100);
     selection(:,rep)=st;
end


%%
cardinality_pfc = mean(sum(selection_pfc));
cardinality_ord = mean(sum(selection));

screening_pfc = mean(sum(selection_pfc(1:4,:))==4);
screening_ord = mean(sum(selection(1:4,:))==4);

insertion_pfc = mean(sum(selection_pfc(5:p,:))~=0);
insertion_ord = mean(sum(selection(5:p,:))~=0);







