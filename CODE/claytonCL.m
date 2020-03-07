function CL = claytonCL(theta,data)
% function CL = claytonCL(theta,data)
%
% The negative copula log-likelihood of a 
% member of Clayton's family
% From Joe(1997), p141
%
% Friday, 29 Sep, 2000
%
% Andrew Patton

% INPUTS: theta 
%				data = [U V];


% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton

if theta~=0;
    CL = log(data(:,1).*data(:,2))*(1+theta);
    CL = CL + (2+1/theta).*log( (data(:,1).^(-theta)) + (data(:,2).^(-theta)) -1);
    CL = log(1+theta) - CL;
    CL = sum(CL);
    CL = -CL;
else
    CL = 0;  % under independence the copula pdf=1, so the log-likelihood=0
end

