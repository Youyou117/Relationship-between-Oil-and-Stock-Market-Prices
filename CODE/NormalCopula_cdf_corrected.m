% Computes the value of the Normal copula at a specified point
% 
% INPUTS:	U,   a Tx1 vector (or a scalar) of F(X[t])
%				V,   a Tx1 vector (or a scalat) of G(Y[t])
%				RHO, a Tx1 vector (or a scalar) of correlation coefficients
%
% Monday, 13 Nov, 2000
%
% Andrew Patton

function out1 = NormalCopula(data,RHO)

U = data(:,1);
V = data(:,2);
T = max([size(U,1),size(V,1),size(RHO,1)]);
N = max([size(U,2),size(V,2),size(RHO,2)]);
out1 = -999.99*ones(T,N);

% stretching the input vectors to match
if size(U,1)<T;
   U = ones(T,1)*U(1,:);
end
if size(V,1)<T;
   V = ones(T,1)*V(1,:);
end
if size(RHO,1)<T;
   RHO = ones(T,1)*RHO(1,:);
end
if size(U,2)<N;
   U = U(:,1)*ones(1,N);
end
if size(V,2)<N;
   V = V(:,1)*ones(1,N);
end
if size(RHO,2)<N;
   RHO = RHO(:,1)*ones(1,N);
end

X = norminv(U,0,1);
Y = norminv(V,0,1);

out1 = bivarnormcdf(X,Y,RHO);
