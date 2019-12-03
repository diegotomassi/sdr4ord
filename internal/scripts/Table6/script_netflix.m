% This script illustrates results reported in Section 6.2 of:
% Forzani, Llop, Tomassi and Arancibia: Sufficient Dimension Reduction for
% Ordinal Predictors.
%
% Requires: 
% - (PROVIDED) LDR package, by RD Cook,L. Forzani and D. Tomassi:
% http://www.jstatsoft.org/article/view/v039i03
% - LIBSVM: https://www.csie.ntu.edu.tw/~cjlin/libsvm/
%
% IMPORTANT: this script can take several weeks to complete
% ========================================================================
disp('WARNING: this script can take several weeks to complete');

[Y,X] = getDATA('pelis20_netflix.txt');
data = [Y,X];
[n,p] = size(data);

%% Results referred to as MEAN
res = zeros(n,p);
for j=1:p,
    Xj = data(:,j);
    for i=1:n,
        xx = Xj; xx(i)=[];
        res(i,j) = Xj(i) - mean(xx);
    end
end
res = res(:);
err_mean = res'*res/(n*p)


%% Results referred to as FULLlin
res_FULLlin = zeros(n,p);
for j = 1:p,
    j
    Y =  data(:,j);
    X = data; X(:,j) = [];
    X = X - repmat(mean(X),n,1);
    for i=1:n,
        yy = Y; yy(i) = []; ytest = Y(i);
        xx = X; xx(i,:) = [];
        xtest = X(i,:);
        %=== linear model for prediction ====
        m = glmfit(xx,yy,'normal');
        yhat = glmval(m,xtest,'identity');
        res_FULLlin(i,j) = ytest-yhat;
    end
end
res = res_FULLlin(:);
err_FULLlin = res'*res/(n*p)



%% Results referred to as PFClin
beta_pfc = cell(n,p);
res_PFClin = zeros(n,p,4);
for j = 1:p,
    Y =  data(:,j);
    X = data; X(:,j) = [];
    X = X - repmat(mean(X),n,1);
    for i=1:n,
        i
        yy = Y; yy(i) = []; ytest = Y(i);
        xx = X; xx(i,:) = [];
        xtest = X(i,:);
        [~,W] = ldr(yy,xx,'PFC','disc',4);
        beta_pfc{i,j}=W;
        for dim=1:4,
            %=== linear model for prediction ====
            m = glmfit(xx*W(:,1:dim),yy,'normal');
            yhat = glmval(m,xtest*W(:,1:dim),'identity');
            res_PFClin(i,j,dim) = ytest-yhat;
        end
    end
end
for dim=1:4,
    res = res_PFClin(:,:,dim); res = res(:);
    err_PFClin(dim) = res'*res/(n*p);
end




%% Results referred to as PFCsvr
beta_pfc = cell(n,p);
res_PFClin = zeros(n,p,4);
for j = 1:p,
    Y =  data(:,j);
    X = data; X(:,j) = [];
    X = X - repmat(mean(X),n,1);
    for i=1:n,
        yy = Y; yy(i) = []; ytest = Y(i);
        xx = X; xx(i,:) = [];
        xtest = X(i,:);
        [~,W] = ldr(yy,xx,'PFC','disc',4);
        beta_pfc{i,j}=W;
        for dim=1:4,
        %== SVR for prediction (non-optimized)
            m = svmtrain(yy,xx*W(:,1:dim),'-s 3 -t 2 -c 4 -g 1 -p .1');
            yhat = svmpredict(ytest,xtest*W(:,1:dim),m);
            res_PFCsvr(i,j,dim) = ytest-yhat;
        end
    end
end
for dim=1:4,
    res = res_PFCsvr(:,:,dim); res = res(:);
    err_PFCsvr(dim) = res'*res/(n*p);
end





%% Results referred to as ORDlin
beta_ord = cell(n,p,4);
res_ORDlin = zeros(n,p,4);

for j=1:p,
    Y =  data(:,j);
    X = data; X(:,j) = [];
  for i=1:n,
    yy = Y; yy(i) = []; ytest = Y(i);
    xx = X; xx(i,:) = [];
    xtest = X(i,:);
    Fy = get_fy(yy,4);
    for dim=1:4,
        [alpha0b,Delta0b,A0,B0] = get_initval_v2(yy,xx,dim,r);
        [A0,Delta0] = rescalePars_v2(A0,Delta0b);
        [~,~,~,~,~,W,~] = EM4ordinal_wcise(xx,yy,dim,alpha0b,Delta0,A0,B0);
        beta_ord{i,j,dim}=W;
        %=== linear model for prediction ====
        m = glmfit(xx*W,yy,'normal');
        yhat = glmval(m,xtest*W,'identity');
        res_ORDlin(ij,dim) = ytest-yhat;
    end
  end
end
  
for dim=1:4,
    res = res_ORDlin(:,:,dim); res = res(:);
    err_ORDlin(dim) = res'*res/(n*p);
end



%% Results referred to as ORDlin
res_ORDsvr = zeros(n,p,4);

for j=1:p,
    Y =  data(:,j);
    X = data; X(:,j) = [];
  for i=1:n,
    yy = Y; yy(i) = []; ytest = Y(i);
    xx = X; xx(i,:) = [];
    xtest = X(i,:);
    Fy = get_fy(yy,4);
    for dim=1:4,
        W = beta_ord{i,j,dim};
        %== SVR for prediction (non-optimized)
        m = svmtrain(yy,xx*W,'-s 3 -t 2 -c 4 -g 1 -p .1');
        yhat = svmpredict(ytest,xtest*W,m);
        res_ORDsvr(i,j,dim) = ytest-yhat;
    end
  end
end
for dim=1:4,
    res = res_ORDsvr(:,:,dim); res = res(:);
    err_ORDsvr(dim) = res'*res/(n*p);
end


