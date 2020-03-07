function [LL,muhat,hhat,uhat] = skewtdis_ARMA_GJRGARCH_LL(theta,data,ARp,MAq,archP,tarchO,garchQ);
% function skewtdis_ARMA_GJRGARCH_LL(theta,ar,ma,arch,tarch,garch,data);
%
%  Function to return the log-likelihood of the standardized residuals from
%  an ARMA-GJRGARCH model with Hansen's skew t residuals. 
%
%  This function is used to construct standard errors for the skew t dis parameters,
%  controlling for the estimation error from the mean and vol estimation
%
%  INPUTS:  theta, a 1+ARp+MAq+archP+tarchO+garchQ+2 vector, containing all parameters
%           data, a Tx1 vector, the data
%           ARp, a scalar, the order of the AR part
%           MAq, a scalar, the order of the MA part
%           archP, a scalar, the order of the ARCH part
%           tarchO, a scalar, the order of the threshold ARCH part
%           garchQ, a scalar, the order of the GARCH part
%
% OUTPUTS:  LL, the negative log-likelihood of this model
%           muhat, a Tx1 vector, the fitted conditional mean
%           hhat, a Tx1 vector, the fitted conditional variance
%           uhat, a Tx1 vector, the estimated prob integral transforms (needs to be FOURTH output to work in two-stage std error code)
%
%  Andrew Patton
%
%  31 March 2011

theta1 = theta(1:1+ARp+MAq);
theta2 = theta(1+ARp+MAq+ 1:1+ARp+MAq+ 1+archP+tarchO+garchQ);
theta3 = theta(end-1:end);

muhat = ARMA_mean(theta1,ARp,MAq,data);
resids = data - muhat;
hhat = GJR_GARCH_vol(theta2,archP,tarchO,garchQ,resids);
stdresids = resids./sqrt(hhat);

LL = skewtdis_LL(theta3, stdresids);

if nargout>3
    uhat = skewtdis_cdf(stdresids,theta3(1),theta3(2));
end