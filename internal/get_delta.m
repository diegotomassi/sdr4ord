function Delta = get_delta(X,Y,Fy)
    SIGMAfit = get_fitted_cov(Y,X,Fy);
    SIGMA = get_cov(X);
    Delta = SIGMA-SIGMAfit;
    Delta = .5*(Delta + Delta');