function [theta,sig2,vcv,resids,yhat] = ARMAX_2(depvar,ar,ma)
%function [theta,sig2,vcv,resids,yhat] = ARMAX_2(depvar,ar,ma);
%
% This program is just a nicer front-end to the existing ARMAX function.
%
% INPUTS:	depvar, a Tx1 vector of data
%				ar,		number of AR lags to use
%				ma, 		number of MA lags to use
%
% OUTPUTS:	theta, the (1+ar+ma) vector of parameter estimates
%					sig2 ,	the esimated variance of the residuals
%					vcv, the variance-covarance matrix of the parameter estimates
%
% Friday, 20 October, 2000
%  -- updated 25 June 2004
%
%  Andrew Patton
%
% see also: ARMAX_opt.m, armax

T = size(depvar,1);
resids=-999.99;
if ar==0 && ma==0
    theta = mean(depvar);
    sig2 = cov(depvar);
    vcv = sig2/size(depvar,1);
    resids = depvar-theta;
    yhat = mean(depvar);
else
    meandepvar = mean(depvar);
    depvar = depvar - meandepvar;
    M = armax(depvar,[ar,ma]);
    
    theta = -999.99*ones(1+ar+ma,1);
    theta(1) = meandepvar;
    
    
    if ar>0
        theta(2:ar+1) = -M.A(2:end);  % need to take negative for AR coefficients
    end
    if ma>0
        theta(ar+2:ar+ma+1) = M.C(2:end);
    end
    sig2 = M.NoiseVariance;
    vcv = [];        % new version of armax.m does not seem to return the covariance matrix for the parameter estimates
    
    %new section
    if nargout>3
        if ma==0
            resids = depvar(ar+1:T);
            for jj = 1:ar
                resids = resids - theta(1+jj)*depvar(ar+1-jj:T-jj);
            end
        else
            resids = zeros(T,1);
            for tt = max(ar,ma)+1:T
                resids(tt) = depvar(tt);
                for jj = 1:ar
                    resids(tt) = resids(tt) - theta(1+jj)*depvar(tt-jj);
                end
                for jj = 1:ma
                    resids(tt) = resids(tt) - theta(1+ar+jj)*resids(tt-jj);
                end
            end
            resids = resids(max(ar,ma)+1:T);
        end
    end
end
yhat = [mean(depvar)*ones(T-length(resids),1);depvar(T-length(resids)+1:end)-resids];   % making sure yhat vector is the right length, using sample mean for first few values if needed
resids = [zeros(T-length(resids),1);resids];  % making sure resids vector is the right length

