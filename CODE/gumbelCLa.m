function CLa = gumbelCL(theta,data)
% function CLa = gumbelCL(theta,data)
%
% The negative copula log-likelihood of a  member of Gumbel's family
% Taken from Joe (1997), p142.
%
% Thursday, 27 July, 2000
%
% Andrew Patton
%
% INPUTS: theta ;
%				data = [U V];

% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556.
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168.
%
% http://fmg.lse.ac.uk/~patton



ut = -log(data(:,1));
vt = -log(data(:,2));

CL = log(gumbel_cdf(data(:,1),data(:,2),theta)) - log(data(:,1)) - log(data(:,2));
CL = CL + (theta-1)*(log(ut)+log(vt)) - (2-1/theta)*(log(ut.^theta+vt.^theta));
CL = CL + log((ut.^theta + vt.^theta).^(1/theta) + theta - 1);
CLa = -CL;
