function [dcor,a,b] = distcorr(x,y)

% This function calculates the distance correlation between x and y.

% dcor: value of distance correlation between x and y 
% a: distances between all points of x
% b: distances between all points of y

% ============================================================
%
% LA IMPLEMENTACION SE CORRESPONDE CON EL PAPER DE SZEKELY, RIZZO Y BAKIROV
%
% ============================================================

% Check if the sizes of the inputs match
%
% ==============================================================================
% QUERES MEDIR SI HAY DEPENDENCIA ENTRE X E Y, ENTONCES SE MUESTREARON JUNTAS.
% LUEGO, SI O SI NECESITAN TENER EL MISMO TAMAÑO.
% ==============================================================================
if size(x,1) ~= size(y,1);
    error('Inputs must have the same number of rows')
end

% ==============================================================================
% ELIMINA DATOS FALTANTES PARA QUE NO JOROBEN AL HACER LAS CUENTAS
% ==============================================================================
% Delete rows containing unobserved values
N = any([isnan(x) isnan(y)],2);
x(N,:) = [];
y(N,:) = [];


% ==============================================================================
% AQUI COMIENZA EL COMPUTO. LA FORMULITA QUE QUIERE IMPLEMENTAR ES:
%
% V(X,Y) = SUM_J SUM_K (A_jk*B_jk)
%
% con A_jk = |X_j-X_k| - .... (ver formulas 2.7 a 2.9 en 
% http://projecteuclid.org/download/pdfview_1/euclid.aos/1201012979)
% ==============================================================================
% Calculate doubly centered distance matrices for x and y
a = pdist2(x,x);
mcol = mean(a);
mrow = mean(a,2);
ajbar = ones(size(mrow))*mcol;
akbar = mrow*ones(size(mcol));
abar = mean(mean(a))*ones(size(a));
A = a - ajbar - akbar + abar;

% ==============================================================================
% ACA REPITE PARA Y LO QUE HIZO PARA X 
% ==============================================================================
b = pdist2(y,y);
mcol = mean(b);
mrow = mean(b,2);
bjbar = ones(size(mrow))*mcol;
bkbar = mrow*ones(size(mcol));
bbar = mean(mean(b))*ones(size(b));
B = b - bjbar - bkbar + bbar;


% Calculate squared sample distance covariance and variances
dcov = sum(sum(A.*B))/(size(mrow,1)^2);

dvarx = sum(sum(A.*A))/(size(mrow,1)^2);
dvary = sum(sum(B.*B))/(size(mrow,1)^2);

% Calculate the distance correlation
dcor = sqrt(dcov/sqrt(dvarx*dvary));