function out1 = BB7UgivenV_t(u,v,t,k,g)
% function out1 = BB7UgivenV_t(u,v,t,k,g)
% 
% Computes the value of the conditional BB7 (U given V)
% copula cdf at a specified point
% 
% INPUTS:	u, a scalar of F(X[t])
%				v, a scalar of G(Y[t])
%				t, the value of C(u|v)
%				k, a scalar of kappa
%				g, a scalar of gamma
%
% Tuesday, 14 Nov, 2000
%
% Andrew Patton

% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton


out1 = (1-(-1+(1-(1-u).^k).^(-g)+(1-(1-v).^k).^(-g)).^(-1./g)).^(-1+1./k).*(-1+(1-(1-u).^k).^...
   (-g)+(1-(1-v).^k).^(-g)).^(-1-1./g).*(1-(1-v).^k).^(-1-g).*(1-v).^(-1+k);

%out1 = (1-(-1+(1-(1-u)^k)^(-g)+(1-(1-v)^k)^(-g))^(-1/g))^(-1+1/k)*(-1+(1-(1-u)^k)^(-g)+(1-(1-v)^k)^(-g))^(-1-1/g)*(1-(1-v)^k)^(-1-g)*(1-v)^(-1+k);
out1 = out1 - t;
