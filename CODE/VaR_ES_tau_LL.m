function Eloss = VaR_ES_tau_LL(theta,Y,alpha,tau)
% function Eloss = VaR_ES_tau_LL(theta,Y,alpha,tau)
%
% Function to return the average loss using the VaR_ES_loss function,
% imposing some smoothing (via parameter tau). This is for the *constant*
% case.
% 
% INPUTS:   theta = [VaR,ES]
%           Y, a Tx1 vector, returns on the target variable
%           alpha, a scalar in (0,1), the probability level
%           tau, a scalar, a positive number if want to smooth the hit variable (eg, tau=10), use -1 (default) if want NO smoothing and want to use the true objective function
%
% OUTPUTS:  Eloss, a scalar, the average loss using these parameter values
%
%  Andrew Patton
%
%  20 Feb 2016
%
% This code was used in: 
% Patton, A.J., J.F. Ziegel, and R. Chen, 2018, Dynamic Semiparametric
% Models for Expected Shortfall (and Value-at-Risk), Journal of Econometrics, forthcoming. 


T = length(Y);
VEhat = nan(T,2);
if nargin<4 || isempty(tau)
    tau = -1;  % default option is to NOT use any smoothing of the indicator function
end

if tau==-1
    hitS = (Y<theta(1));  % then use the actual indicator function
else
    hitS = 1./(1+exp( tau*(Y-theta(1)) ));
end

loss = -1/alpha/theta(2)*hitS.*(theta(1)-Y) - 1/theta(2)*(theta(2)-theta(1)) + log(-theta(2));
Eloss = mean(loss);
