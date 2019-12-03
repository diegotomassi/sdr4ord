function [Eznew_ij,Ezznew_ijj] = updateEzANDEzz(Ez,Ezz,Theta,X,Fy,Delta,alpha,idx)

i = idx(1); j= idx(2);
Xij = X(i,j);
thetaj = Theta{j};
Kj = length(thetaj)+1;

beta = Delta*alpha;

Sigma_jmj = Delta(j,:); Sigma_jmj(j) = [];
Sigma_mjmj = Delta; Sigma_mjmj(j,:) = []; Sigma_mjmj(:,j) = [];

sigmatilde_ij = sqrt(Delta(j,j) - Sigma_jmj*inv(Sigma_mjmj)*Sigma_jmj');

fyi = Fy(i,:)';
Mumj = beta*fyi; 
Muj = Mumj(j); Mumj(j)=[];
Ez_imj = Ez(i,:); Ez_imj(j) = [];
Ezz_imjmj = squeeze(Ezz(i,:,:)); Ezz_imjmj(j,:)=[]; Ezz_imjmj(:,j)=[];
Emutilde_ij = Sigma_jmj*inv(Sigma_mjmj)*(Ez_imj' - Mumj) + Muj; 


if Xij==1,
    deltatilde_ij = (thetaj(Xij) - Emutilde_ij)/sigmatilde_ij;    
    coef = -1*normpdf(deltatilde_ij)/normcdf(deltatilde_ij);
    if ~isfinite(coef),
%         disp('INF coef found')        
        coef = 0;
    end
    
    % actualizamos E(Zij | todo)
    Eznew_ij = Emutilde_ij + coef * sigmatilde_ij;

    % actualizamos E(Zij*Zij | todo)     
    coef2 = (0 - deltatilde_ij*normpdf(deltatilde_ij))/(normcdf(deltatilde_ij) - 0);
    if ~isfinite(coef),
%         disp('INF coef2 found')        
        coef2 = 0;
    end
    Ezznew_ijj = sigmatilde_ij^2 + Emutilde_ij^2 + ...
                 2*coef*Emutilde_ij*sigmatilde_ij + ...
                 coef2*sigmatilde_ij^2;
elseif Xij==Kj,
    deltatilde_ijm = (thetaj(Xij-1) - Emutilde_ij)/sigmatilde_ij;
    coef = normpdf(deltatilde_ijm)/(1-normcdf(deltatilde_ijm));
    if ~isfinite(coef),
%         disp('INF coef found')        
        coef = 0;
    end
    
    % actualizamos E(Zij | todo)
    Eznew_ij = Emutilde_ij + coef * sigmatilde_ij;

    % actualizamos E(Zij*Zij | todo)         
    coef2 = (deltatilde_ijm*normpdf(deltatilde_ijm) - 0)/(1 - normcdf(deltatilde_ijm));
    if ~isfinite(coef2),
%         disp('INF coef2 found')        
        coef2=0;
    end
    Ezznew_ijj = sigmatilde_ij^2 + Emutilde_ij^2 + ...
                 2*coef*Emutilde_ij*sigmatilde_ij + ...
                 coef2*sigmatilde_ij^2;
    
else
    deltatilde_ijm = (thetaj(Xij-1) - Emutilde_ij)/sigmatilde_ij;
    deltatilde_ij = (thetaj(Xij) - Emutilde_ij)/sigmatilde_ij;
    coef = (normpdf(deltatilde_ijm) - normpdf(deltatilde_ij)) / (normcdf(deltatilde_ij) - normcdf(deltatilde_ijm));
    if ~isfinite(coef),
%         disp('INF coef found')
        coef=0;
    end
    % actualizamos E(Zij | todo)
    Eznew_ij = Emutilde_ij + coef * sigmatilde_ij;

    % actualizamos E(Zij*Zij | todo)         
    coef2 = (deltatilde_ijm*normpdf(deltatilde_ijm) - deltatilde_ij*normpdf(deltatilde_ij))/(normcdf(deltatilde_ij) - normcdf(deltatilde_ijm));
    if ~isfinite(coef2),
%         disp('INF coef2 found')
        coef2 = 0;
    end
    Ezznew_ijj = sigmatilde_ij^2 + Emutilde_ij^2 + ...
                 2*coef*Emutilde_ij*sigmatilde_ij + ...
                 coef2*sigmatilde_ij^2;
        
end
        
