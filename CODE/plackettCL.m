function CL = plackettCL(theta,data)
% function CL = plackettCL(theta,data)
%
% The negative copula log-likelihood of a 
% Plackett copula
%
% Thursday, 24 May, 2001.
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


if theta<=0
   CL = 1e6;
elseif isreal(theta)==0
   CL=1e7;
else
   CL = plackett_pdf(data(:,1),data(:,2),theta);
   CL = log(CL);  % lazy way to get a log-likelihood!
   CL = -sum(CL);
end

