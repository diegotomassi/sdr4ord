function [mse,sd] = cvOrig1(Y,X,model,nfold)
if nargin<3,
    nfold=10;
end
N = length(Y);
blksize = floor(N/nfold);
predictions = zeros(1,nfold);

shuf = randperm(N);
Y = Y(shuf); X = X(shuf,:);
parfor k=1:nfold,
    % set training and testing partitions
    if k < nfold,
        idx = ((k-1)*blksize + 1):(k*blksize);
    else
        idx = ((k-1)*blksize + 1):N;
    end
    xts = X(idx,:); yts = Y(idx,:);
    xtr = X; xtr(idx,:) = [];   ytr = Y; ytr(idx)=[];

    if strcmpi(model,'linear'),   % fil a linear model
        m = glmfit(xtr,ytr,'normal') ;
        yhat = glmval(m,xts,'identity');
        er = yts-yhat;
        predictions(k)= er'*er/length(yts);
    elseif strcmpi(model,'logit'),
        m = glmfit(xtr,ytr,'binomial');
        yhat = glmval(m,xts,'logit') > 0.5;
        predictions(k) = mean(yts~=yhat);
    end
end
mse = mean(predictions)
sd = std(predictions)