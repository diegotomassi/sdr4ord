function [S,M] = computeEstep(X,Fy,Delta, alpha, Theta,type)
if nargin < 6,
    type = 'levina';
end

switch type,
    case 'levina'
        [S,M] = levinaEstep(X,Fy,Delta,alpha,Theta);
    otherwise
        error('Selected method for E-step is not available')
end


%---------------------------------------------
function [S,M] = levinaEstep(X,Fy,Delta,alpha,Theta)

[n,p] = size(X);
Ez = zeros(n,p);
Ezz = zeros(n,p,p);

Eznew = zeros(n,p);
Ezznew = zeros(n,p,p);

% Initialization
StopNotMet = 1;
history = ones(1,5);
Sold = eye(size(Delta));
iter = 0;

while(StopNotMet && iter<5),
    iter = iter + 1;
    for i=1:n,
        for j=1:p,
            [Eznew(i,j),Ezznew(i,j,j)] = updateEzANDEzz(Ez,Ezz,Theta,X,Fy,Delta,alpha,[i,j]);
        end
        for j=1:p,
            for jprima = 1:p,
                if j~=jprima,
                    Ezznew(i,j,jprima) = Eznew(i,j)*Eznew(i,jprima);
                end
            end
        end
    end
    % Update S and M
    for j = 1:p,
        for jprima = 1:p,
            S(j,jprima) = sum(Ezznew(:,j,jprima));
        end
    end
    M = Eznew;
    
    Ezz = Ezznew;
    Ez = Eznew;
    
    [StopNotMet,history] = checkConvergenceIn(S,Sold,history);
    Sold = S;
end
            
            
            
