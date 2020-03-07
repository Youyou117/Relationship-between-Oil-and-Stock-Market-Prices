function [CL, rhohat] = bivnorm_tvp1_CL(theta,Zdata,rhobar)
% function [CL, rhohat] = bivnorm_tvp1_CL(theta,Zdata,rhobar)
%
% The negative copula log-likelihood of a 
% bivariate Normal distribution
% with time-varying correlation
%
% Monday, 4 Sep, 2000
%
% Andrew Patton
%
% INPUTS: theta ;
%				data = [U V];
%				rhobar = parameter of the Normal copula without time-variation

% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton



% Below I use an average of the last 10 lags rather than just the last lag.
% This can of course be changed

T = size(Zdata,1);
x = norminv(Zdata(:,1),0,1);
y = norminv(Zdata(:,2),0,1);

kappa = -999.99*ones(T,1);
kappa(1) = rhobar;			% this is the MLE of kappa in the time-invariant version of this model
for jj = 2:T
    if jj<=10
        psi = theta(1) + theta(2)*mean(x(1:jj-1).*y(1:jj-1)) + theta(3)*kappa(jj-1);
    else
        psi = theta(1) + theta(2)*mean(x(jj-10:jj-1).*y(jj-10:jj-1)) + theta(3)*kappa(jj-1);
    end
   kappa(jj) = 1.998/(1+exp(-psi))-0.999;		% a modified logistic transformation
end
rhohat = kappa;  % time-path of conditional copula parameter

CL = -1*(2*(1-kappa.^2)).^(-1).*(x.^2+y.^2-2*kappa.*x.*y);
CL = CL + 0.5*(x.^2+y.^2);  
CL = CL - 0.5*log(1-kappa.^2);
CL = sum(CL);
CL = -CL;
