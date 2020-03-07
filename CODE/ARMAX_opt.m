function [theta,sig2,vcv,order,resids,yhat] = ARMAX_opt(depvar,maxp,maxq,crit)
%function [theta,sig2,vcv,order,resids] = ARMAX_opt(depvar,maxp,maxq,crit);
%
% This program determines the optimal AR and MA lags, according to a specified information criterion.
%
% INPUTS:	depvar, a Tx1 vector of observations of the variable under analysis
%           maxp, a scalar >=0, the maximum AR lag length
%			maxq, a scalar >=0, the maximum MA lag length
%			crit, a string, either 'AIC' or 'BIC' or 'HQ ' corresponding
%               to which of the Akaike or Schwarz or Hannan-Quinn Information	
%				criteria is to be used.
%
% OUTPUTS:	theta, the (ar+ma+1) vector of parameter estimates from the optimal model
%           sig2 ,	the esimated variance of the residuals from the optimal model
%			vcv, the variance-covarance matrix of the parameter estimates from the optimal model
%			order, the order [p_opt, q_opt] of the optimal ARMA model
%           resids, the vector of resids from the optimal ARMA model
%           yhat, a scalar, the one-step forecast of the depvar using the optimal ARMA model
%
% Tuesday, 26 December, 2000
%
%  Andrew Patton

T = size(depvar,1);
theta = mean(depvar);
sig2 = cov(depvar);
vcv = sig2/T;
order = [0;0];
crit1 = log(sig2);  % with p=q=0 both criteria take the same value
resids = depvar-mean(depvar);
yhat = mean(depvar);

if nargin==3
    crit = 'BIC';
end

for ii = maxp:-1:0
    for jj = maxq:-1:0
        [theta_1,sig2_1,vcv_1,resids_1,yhat_1] = ARMAX_2(depvar,ii,jj);
        sig2_1 = mean(resids_1.^2);  % AJP, 31mar11: using usual definition of MSE here, rather than "armax" built-in
        if strcmp('AIC',crit)
            crit2 = log(sig2_1) + (ii+jj+1)*2/T;
        elseif strcmp('HQ',crit)
            crit2 = log(sig2_1) + (ii+jj+1)*2*log(log(T))/T;  
        else	% crit default is the BIC
            crit2 = log(sig2_1) + (ii+jj+1)*log(T)/T;
        end
        if crit2<crit1
            theta = theta_1;
            sig2 = sig2_1;
            vcv = vcv_1;
            order = [ii;jj];
            crit1 = crit2;
            resids = resids_1;
            yhat = yhat_1;
        end
    end
end

