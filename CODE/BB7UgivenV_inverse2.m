function u = BB7UgivenV_inverse2(v,t,k,g)
% function u = BB7UgivenV_inverse2(v,t,k,g)
%
% Computes the ivnerse value of the conditional BB7 (U given V)
% copula cdf at a specified point
% 
% INPUTS:	t, a scalar uniform r.v.
%				v, a scalar r.v. G(Y[t])
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

info=struct('maxit',{100},'tol',{10^-6},'pflag',0);
[u,junk] = bisect2('BB7UgivenV_t',10^-4,1-10^-4,info,v,t,k,g);


% dealing with the extreme values that don't converge
if junk==-999.99 		% then the zero of the function is very near to 0
   [u,junk] = bisect2('BB7UgivenV_t',10^-18,10^-4,info,v,t,k,g);
elseif junk==999.99		% then the zero of the funciton is very near to 1
   [u,junk] = bisect2('BB7UgivenV_t',1-10^-4,1-10^-18,info,v,t,k,g);
end
