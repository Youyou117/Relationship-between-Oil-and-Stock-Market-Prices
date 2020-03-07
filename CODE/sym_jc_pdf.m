function CL = sym_jc_pdf(u,v,tauU,tauL)
% function CL = sym_jc_pdf(u,v,tauU,tauL)
%
% Computes the value of the symmetrised Joe-Clayton
% copula pdf at a specified point
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


T = max([size(u,1),size(v,1),size(tauU,1),size(tauL,1)]);
out1 = -999.99*ones(T,1);

% stretching the input vectors to match
if size(u,1)<T;
   u = u*ones(T,1);
end
if size(v,1)<T;
   v = v*ones(T,1);
end
if size(tauU,1)<T;
   tauU = ones(T,1)*tauU(1);
end
if size(tauL,1)<T;
   tauL = ones(T,1)*tauL(1);
end

k1 =  1./log2(2-tauU);
k2 = -1./log2(tauL);
CL1 = ((1 - (1 - u).^k1).^(k2 - 1).* (1 - u).^(k1 - 1).*(-1 + k1.*(k2.* (-1 + (-1 + (1 - (1 - u).^k1).^(-k2) + (1 - (1 - v).^k1).^(-k2)).^(k2.^(-1))) + (-1 + (1 - (1 - u).^k1).^(-k2) + (1 - (1 - v).^k1).^(-k2)).^(k2.^(-1)))).* (1 - (-1 + (1 - (1 - u).^k1).^(-k2) + (1 - (1 - v).^k1).^(-k2)).^(-k2.^(-1))).^(k1.^(-1)).* (1 - (1 - v).^k1).^(k2 - 1).* (1 - v).^(k1 - 1));
CL2 = (((-1 + (-1 + (1 - (1 - u).^k1).^(-k2) + (1 - (1 - v).^k1).^(-k2)).^(k2.^(-1))).^2) .* ((1 - (1 - u).^k1).^k2 + (1 - (1 - v).^k1).^k2 - (1 - (1 - u).^k1).^k2.* (1 - (1 - v).^k1).^k2).^2);
CL1 = CL1./CL2;

k1 =  1./log2(2-tauL);
k2 = -1./log2(tauU);
u  = 1-u;
v  = 1-v;
CL3 = ((1 - (1 - u).^k1).^(k2 - 1).* (1 - u).^(k1 - 1).*(-1 + k1.*(k2.* (-1 + (-1 + (1 - (1 - u).^k1).^(-k2) + (1 - (1 - v).^k1).^(-k2)).^(k2.^(-1))) + (-1 + (1 - (1 - u).^k1).^(-k2) + (1 - (1 - v).^k1).^(-k2)).^(k2.^(-1)))).* (1 - (-1 + (1 - (1 - u).^k1).^(-k2) + (1 - (1 - v).^k1).^(-k2)).^(-k2.^(-1))).^(k1.^(-1)).* (1 - (1 - v).^k1).^(k2 - 1).* (1 - v).^(k1 - 1));
CL4 = (((-1 + (-1 + (1 - (1 - u).^k1).^(-k2) + (1 - (1 - v).^k1).^(-k2)).^(k2.^(-1))).^2) .* ((1 - (1 - u).^k1).^k2 + (1 - (1 - v).^k1).^k2 - (1 - (1 - u).^k1).^k2.* (1 - (1 - v).^k1).^k2).^2);
CL3 = CL3./CL4;
CL  = 0.5*(CL1+CL3);


