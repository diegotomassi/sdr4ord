function [Eznew_ij,Ezznew_ijj] = updateEzANDEzz(Ez,Ezz,Theta,X,Fy,Delta,alpha,idx)

i = idx(1); j= idx(2);
Xij = X(i,j);
thetaj = Theta{j};

beta = Delta*alpha;

Sigma_jmj = Delta(j,:); Sigma_jmj(j) = [];
Sigma_mjmj = Delta; Sigma_mjmj(j,:) = []; Sigma_mjmj(:,j) = [];

sigmatilde_ij = Delta(i,j) - Sigma_jmj*inv(Sigma_mjmj)*Sigma_jmj';

fyi = Fy(i,:)';
Mumj = beta*fyi; 
Muj = Mumj(j); Mumj(j)=[];
Ez_imj = Ez(i,:); Ez_imj(j) = [];
Emutilde_ij = Sigma_jmj*inv(Sigma_mjmj)*(Ez_imj - Mumj) + Muj; 

deltatilde_ijm = (thetaj(Xij-1) - Emutilde_ij)/sigmatilde_ij;
deltatilde_ij = (thetaj(Xij-1) - Emutilde_ij)/sigmatilde_ij;
coef = (normpdf(deltatilde_ijm) - normpdf(deltatilde_ij)) / (normcdf(deltatilde_ij) - normcdf(deltatilde_ijm));

% actualizamos E(Zij | todo)
Eznew_ij = Emutilde_jj + coef * sigmatilde_ij;


% actualizamos E(Zij*Zij | todo)
Ezznew_ijj = Sigma_jmj*inv(Sigma_mjmj)*Ezz(i,j,j)*inv(Sigma_mjmj)*Sigma_jmj' + sigmatilde_ij + ...
      2*coef*Emutilde_ij*sigmatilde_ij + ...
      (deltatilde_ijm*normpdf(deltatilde_ijm) - deltatilde_ij*normpdf(deltatilde_ij)/(normcdf(deltatilde_ij) - normcdf(deltatilde_ijm))*sigmatilde_ij;
        
        

