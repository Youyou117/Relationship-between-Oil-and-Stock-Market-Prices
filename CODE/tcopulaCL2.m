function CL = tcopulaCL2(theta,data)
% function CL = tcopulaCL2(theta,data)
% 
% The negative copula log-likelihood of a Student's t copula
% Differs from tcopulaCL.m only in that second parameter is *inverse*
% degree of freedom (so lies inside (0,0.5) not (2,Inf)
%
% From my own calculations - Roncalli's seems to have some mistakes in it
%
% Sunday, 22 April, 2001, Modified: 1 April 2011
%
%
% Andrew Patton
%
% INPUTS: theta = [rho, nuinv] = [ correlation coeff, *inverse* degree of freedom parameter]
%				data = [U V];
%
%


% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton


rho = theta(1);
nu  = 1/theta(2);

if nu>=100;			% for nu>100, normal is good approximation, and MUCH quicker.
   x = norminv(data(:,1));
   y = norminv(data(:,2));
else
   x = tinv(data(:,1),nu);
   y = tinv(data(:,2),nu);
end

CL = gammaln((nu+2)/2) + gammaln(nu/2) - 2*gammaln((nu+1)/2) - 0.5*log((1-rho^2));
CL = CL - (nu+2)/2*log(1+(x.^2 + y.^2 - 2*rho*x.*y)/(nu*(1-rho^2)));
CL = CL + (nu+1)/2*log(1+x.^2/nu) + (nu+1)/2*log(1+y.^2/nu);
CL = sum(CL);
CL = -CL;