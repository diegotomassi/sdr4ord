function [Theta] = findThresolds(X,Fy,Delta, alpha)
[n,p] = size(X);
unos = ones(n,1);
beta = Delta*alpha;
Theta = cell(p,1); % structure to store the thresholds
for j=1:p,
%     disp(j)
    dj = Delta(j,j);
    betaj = beta(j,:);
    Xj = X(:,j);
    Kj = max(Xj);
    theta_j = zeros(1,Kj-1);
    medias = Fy*betaj';
    for k=1:Kj-1,
        njk = length(find(Xj <= k));
        objfun = @(x)  (njk - sum(normcdf((x*unos - medias)/dj))); % see page 6.
        % We need a starting value or an interval within to search...
%         x0 = dj*norminv(njk) + mean(Fy*betaj');
        x0 = [min(medias)-5*dj , max(medias)+5*dj];
        theta_j(k) = fzero(objfun,x0);
    end
    Theta{j} = theta_j;
end
    