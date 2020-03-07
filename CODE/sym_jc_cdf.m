% Computes the value of the symmetrised Joe-Clayton
% copula at a specified point
% 
% INPUTS:	U, a Tx1 vector (or a scalar) of F(X[t])
%				V, a Tx1 vector (or a scalat) of G(Y[t])
%				tauU, a Tx1 vector (or a scalar) of upper tail dependence measures
%				tauL, a Tx1 vector (or a scalar) of lower tail dependence measures
%
% Monday, 28 January, 2002
%
% Andrew Patton

% Published in:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence,
% International Economic Review, 47(2), 527-556. 
%
% See also:
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of
% Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
%
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and
% Asymmetric Dependence for Asset Allocation, Journal of Financial
% Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton


function out1 = sym_jc_cdf(U,V,tauU,tauL)

T = max([size(U,1),size(V,1),size(tauU,1),size(tauL,1)]);
out1 = -999.99*ones(T,1);

% stretching the input vectors to match
if size(U,1)<T;
   U = ones(T,1)*U(1);
end
if size(V,1)<T;
   V = ones(T,1)*V(1);
end
if size(tauU,1)<T;
   tauU = ones(T,1)*tauU(1);
end
if size(tauL,1)<T;
   tauL = ones(T,1)*tauL(1);
end

K =  1./log2(2-tauU);
G = -1./log2(tauL);
out1 = 1-((1-(((1-((1-U).^K)).^(-G))+((1-((1-V).^K)).^(-G))-1).^(-1./G)).^(1./K));
K =  1./log2(2-tauL); % switching the upper and lower measures
G = -1./log2(tauU);
U = 1-U;
V = 1-V;
out2 = (1-U) + (1-V) - 1 + 1-((1-(((1-((1-U).^K)).^(-G))+((1-((1-V).^K)).^(-G))-1).^(-1./G)).^(1./K));
out1 = 0.5*(out1+out2);

