function [CL, kappa] = Gumbel_tvp1_CL(theta,data,kappabar)
%function [CL, kappa] = Gumbel_tvp1_CL(theta,data,kappabar)
% 
% The negative copula log-likelihood of Gumbel's copula with time-varying parameter 
%
% Sunday, 5 Aug, 2001.
%
% Andrew Patton
%
% INPUTS: theta ;
%				data = [U V], a Tx2 matrix of data
%				thetabar = parameter of the Clayton copula without time-variation
%				exog, a Txk matrix of exogenous conditional correlation regressors

% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton


T = size(data,1);
u = data(:,1);
v = data(:,2);

kappa = -999.99*ones(T,1);
kappa(1) = kappabar;			
for jj = 2:T;
   if jj<=10
      psi1 = theta(1)+ theta(2)*kappa(jj-1) + theta(3)*(mean(abs(u(1:jj-1)-v(1:jj-1))));
   else
      psi1 = theta(1)+ theta(2)*kappa(jj-1) + theta(3)*(mean(abs(u(jj-10:jj-1)-v(jj-10:jj-1))));
   end
   kappa(jj) = 1.0001 + psi1^2;
end

ut = -log(data(:,1));
vt = -log(data(:,2));

CL = log(Gumbel_cdf(u,v,kappa)) - log(u) - log(v);
CL = CL + (kappa-1).*(log(ut)+log(vt)) - (2-1./kappa).*(log(ut.^kappa+vt.^kappa));
CL = CL + log((ut.^kappa + vt.^kappa).^(1./kappa) + kappa - 1);
CL = sum(CL);
CL = -CL;

if isreal(theta)==0;
   CL = 1e6;
elseif isreal(CL)==0;
   CL = 1e7;
elseif isnan(CL)
   CL = 1e8;
elseif isinf(CL)
   CL = 1e9;
end