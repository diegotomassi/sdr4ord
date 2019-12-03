warning off;
GBA = load('traindata_binGBA.txt');
y = GBA(:,1);
X = GBA(:,2:end);
r=3;
dim=1;
[alpha0b,Delta0b,A0,B0] = get_initval_v2(y,X,dim,r);
[A0,Delta0] = rescalePars_v2(A0,Delta0b);
[~,~,~,~,~,alphaoutGBA,~,st] = EM4ordinal_wcise(X,y,dim,alpha0b,Delta0,A0,B0,500,1e-3);
[alphaoutGBA,st] = newThresh(alphaoutGBA,0.05);

%%
Pampeana = load('traindata_binPampeana.txt');
y = Pampeana(:,1);
X = Pampeana(:,2:end);
r=3;
dim=1;
[alpha0b,Delta0b,A0,B0] = get_initval_v2(y,X,dim,r);
[A0,Delta0] = rescalePars_v2(A0,Delta0b);
[~,~,~,~,~,alphaoutPAMPA,~,st] = EM4ordinal_wcise(X,y,dim,alpha0b,Delta0,A0,B0,500,1e-3);
[alphaoutPAMPA,st] = newThresh(alphaoutPAMPA,0.05);

%%
NOA = load('traindata_binNOA.txt');
y = NOA(:,1);
X = NOA(:,2:end);
r=3;
dim=1;
[alpha0b,Delta0b,A0,B0] = get_initval_v2(y,X,dim,r);
[A0,Delta0] = rescalePars_v2(A0,Delta0b);
[~,~,~,~,~,alphaoutNOA,~,st] = EM4ordinal_wcise(X,y,dim,alpha0b,Delta0,A0,B0,500,1e-3);
[alphaoutNOA,st] = newThresh(alphaoutNOA,0.05);


%%
NEA = load('traindata_binNEA.txt');
y = NEA(:,1);
X = NEA(:,2:end);
r=3;
dim=1;
[alpha0b,Delta0b,A0,B0] = get_initval_v2(y,X,dim,r);
[A0,Delta0] = rescalePars_v2(A0,Delta0b);
[~,~,~,~,~,alphaoutNEA,~,st] = EM4ordinal_wcise(X,y,dim,alpha0b,Delta0,A0,B0,500,1e-3);
[alphaoutNEA,st] = newThresh(alphaoutNEA,0.05);

%%
PATAGON = load('traindata_binPatagonia.txt');
y = PATAGON(:,1);
X = PATAGON(:,2:end);
r=3;
dim=1;
[alpha0b,Delta0b,A0,B0] = get_initval_v2(y,X,dim,r);
[A0,Delta0] = rescalePars_v2(A0,Delta0b);
[~,~,~,~,~,alphaoutPATAGON,~,st] = EM4ordinal_wcise(X,y,dim,alpha0b,Delta0,A0,B0,500,1e-3);
[alphaoutPATAGON,st] = newThresh(alphaoutPATAGON,0.05);

%%
loadings = [alphaoutGBA,alphaoutPAMPA,alphaoutNOA,alphaoutNEA,alphaoutPATAGON]


%%
warning on;