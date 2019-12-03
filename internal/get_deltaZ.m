function Delta = get_delta(X,Y)
    Y = grp2idx(Y);
    data_parameters = setdatapars(Y,X,max(Y));
    SIGMAfit = get_average_cov(Y,X,Fy);
    SIGMA = get_cov(X);
    Delta = SIGMA-SIGMAfit;
    Delta = .5*(Delta + Delta');