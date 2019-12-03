function [mse,sd] = cvPFCord(Y,X,dim,r,lambda,model,nfold)
if nargin<7,
    nfold=10;
end
disp('== starting CV procedure for PFCord (this can take several minutes)===')
N = length(Y);
blksize = floor(N/nfold);
predictions = zeros(1,nfold);
shuf = randperm(N);
Y = Y(shuf); X = X(shuf,:);
 
if strcmpi(model,'linear'),   % fit a linear model
  parfor k=1:nfold,
    % set training and testing partitions
    if k < nfold,
        idx = ((k-1)*blksize + 1):(k*blksize);
    else
        idx = ((k-1)*blksize + 1):N;
    end
    xts = X(idx,:); yts = Y(idx,:);
    xtr = X; xtr(idx,:) = [];   ytr = Y; ytr(idx)=[];
    
    % estimate the reduction
    [alpha0b,Delta0b,A0,B0] = get_initval_v2(ytr,xtr,dim,r);
    [A0,Delta0] = rescalePars_v2(A0,Delta0b);
    [~,~,~,~,~,alphaout,~] = EM4ordinal_wcise(xtr,ytr,dim,alpha0b,Delta0,A0,B0,lambda);

     m = glmfit(xtr*alphaout,ytr,'normal') ;
     yhat = glmval(m,xts*alphaout,'identity');
     er = yts-yhat;
     predictions(k)= er'*er/length(yts);
  end
elseif strcmpi(model,'logit'),
  parfor k=1:nfold,
    % set training and testing partitions
    if k < nfold,
        idx = ((k-1)*blksize + 1):(k*blksize);
    else
        idx = ((k-1)*blksize + 1):N;
    end
    xts = X(idx,:); yts = Y(idx,:);
    xtr = X; xtr(idx,:) = [];   ytr = Y; ytr(idx)=[];
    % estimate the reduction
    [alpha0b,Delta0b,A0,B0] = get_initval_v2(ytr,xtr,dim,r);
    [A0,Delta0] = rescalePars_v2(A0,Delta0b);
    [~,~,~,~,~,alphaout,~] = EM4ordinal_wcise(xtr,ytr,dim,alpha0b,Delta0,A0,B0,lambda);
     m = glmfit(xtr*alphaout,ytr,'binomial');
     yhat = glmval(m,xts*alphaout,'logit') > 0.5;
     predictions(k) = mean(yts~=yhat);
  end
end
mse = mean(predictions);
sd = std(predictions);
disp('Estimated prediction error from CV procedure is:')
disp(mse)

