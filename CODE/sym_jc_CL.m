function CL = bb7CL(theta,data)
% function CL = bb7CL(theta,data)
%
% The negative copula log-likelihood of a 
% member of the symmetrised Joe-Clayton family
% From Joe(1997), p153
%
% Monday, 28 January, 2002
%
% Andrew Patton
%
% INPUTS: theta = [tauU,tauL]
%				data = [U V];

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


CL = sym_jc_pdf(data(:,1),data(:,2),theta(1),theta(2));
CL = log(CL);		% lazy way of obtaining the log-likelihood.


CL = sum(CL);
CL = -CL;
%[CL]
if isnan(CL)
   CL = 1e6;
end
if isreal(CL)==0
   CL = 1e7;
end
if isreal(theta)==0
   CL = 1e8;
end