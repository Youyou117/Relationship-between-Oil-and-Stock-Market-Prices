% Computes the value of the Normal copula at a specified point
% 
% INPUTS:	U,   a Tx1 vector (or a scalar) of F(X[t])
%				V,   a Tx1 vector (or a scalat) of G(Y[t])
%				RHO, a Tx1 vector (or a scalar) of correlation coefficients
%
% Monday, 13 Nov, 2000
%
% Andrew Patton

function out1 = NormalCopula_pdf(U,V,RHO)

T = max([size(U,1),size(V,1),size(RHO,1)]);
out1 = -999.99*ones(T,1);

% stretching the input vectors to match
if size(U,1)<T;
   U = U*ones(T,1);
end
if size(V,1)<T;
   V = V*ones(T,1);
end
if size(RHO,1)<T;
   RHO = RHO*ones(T,1);
end

X = norminv(U,0,1);
Y = norminv(V,0,1);

out1 = (1-RHO.^2).^(-0.5).*exp(-0.5*(1-RHO.^2).^(-1).*(X.^2+Y.^2-2*RHO.*X.*Y)).*exp(0.5*(X.^2+Y.^2));
