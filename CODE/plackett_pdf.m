function pdf = plackett_pdf(u,v,k)
% function pdf = plackett_pdf(u,v,k)
%
% Computes the value of the Plackett copula pdf at a specified point
% 
% INPUTS:	U, a Tx1 vector (or a scalar) of F(X[t])
%				V, a Tx1 vector (or a scalat) of G(Y[t])
%				K, a Tx1 vector (or a scalar) of kappas
%
% Thursday, 24 May, 2001.
%
% Andrew Patton


% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton


T = max([size(u,1),size(v,1),size(k,1)]);

% stretching the input vectors to match
if size(u,1)<T;
   u = u*ones(T,1);
end
if size(v,1)<T;
   v = v*ones(T,1);
end
if size(k,1)<T;
   k = k*ones(T,1);
end

pdf = k.*(1+(u - 2*u.*v + v).*(k-1))./((((1+(k-1).*(u+v)).^2) - 4*u.*v.*k.*(k-1)).^(3/2));