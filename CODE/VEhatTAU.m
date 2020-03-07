function VEhat = VEhatTAU(data,alpha,tau)
% function VEhat = VEhatTAU(data,alpha,tau);
%
% Function to estimate the sample "VaR" and "ES" using the *smoothed" FZ0
% loss function (so these quantities are not really VaR and ES, unless
% tau->+inf)
%
%  INPUTS:  data, a Txk matrix of data, organized in columns
%           alpha, a scalar or 1xk vector, in (0,1), the probability level
%           tau, a scalar or 1xk vector, the smoothing parameter (set to -1 if want no smoothing)
%
% OUTPUTS:  VEhat, a kx2 matrix of estimated VaR and ES
%
%  Andrew Patton
%
%  20 Feb 2016
%
% This code was used in: 
% Patton, A.J., J.F. Ziegel, and R. Chen, 2018, Dynamic Semiparametric
% Models for Expected Shortfall (and Value-at-Risk), Journal of Econometrics, forthcoming. 


T = size(data,1);
k = max([size(data,2),length(alpha),length(tau)]);
if size(data,2)<k
    data = data*ones(1,k);  % making shape conform
end
if length(alpha)<k
    alpha = alpha(1)*ones(1,k);
end
if length(tau)<k
    tau = tau(1)*ones(1,k);
end


options = optimset('Display','off','TolCon',10^-6,'TolFun',10^-6,'TolX',10^-4,'DiffMaxChange',Inf,'DiffMinChange',0);


VEhat = nan(k,2);
for kk=1:k
    VEhat0 = sample_VE(data(:,kk),alpha(kk));  % this is the sample VaR and ES (no smoothing). Use as starting value
    warning off;  % get a warning about not providing a gradient - turning this off.
    thetahat1 = fminunc('VaR_ES_tau_LL',VEhat0,options,data(:,kk),alpha,tau(kk));
    VEhat(kk,:) = thetahat1';
end

