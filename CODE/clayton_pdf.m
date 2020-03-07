function pdf = clayton_pdf(u,v,k1)
% function pdf = clayton_pdf(u,v,k1)
%
% Computes the value of the Clayton copula pdf at a specified point
% 
% INPUTS:	U, a Tx1 vector (or a scalar) of F(X[t])
%				V, a Tx1 vector (or a scalat) of G(Y[t])
%				K, a Tx1 vector (or a scalar) of kappas
%
% Monday, 7 May, 2001.
%
% Andrew Patton
%
%  NOTE:  only works numerically for k1<=34 or so. 
%
% FURTHER NOTE: when k1=0 this copula converges to the "independence
% copula", which needs to be computed separately. (this copula pdf takes
% the value of 1 for all (u,v) in [0,1]x[0,1]

% lower tail dep = 2^(-1/K) for Clayton's copula
% upper tail dep = 0


% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton

T = max([size(u,1),size(v,1),size(k1,1)]);

% stretching the input vectors to match
if size(u,1)<T;
   u = u*ones(T,1);
end
if size(v,1)<T;
   v = v*ones(T,1);
end
if size(k1,1)<T;
   k1 = k1*ones(T,1);
end


pdf = (1+k1).*(u.*v).^(-k1-1).*(u.^(-k1)+v.^(-k1)-1).^(-2-k1.^(-1));


