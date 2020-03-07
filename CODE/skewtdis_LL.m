function LL = skewtdis_LL(theta, x)
% PURPOSE: returns the log-likelihood at x of Hansen's (1994) 'skewed t' distribution
%---------------------------------------------------
% USAGE: LL = skewtdis_LL(x,nu,lambda)
% where: x  = a vector of data
%        theta = [nu;lambda]
%					nu = degrees of freedom parameter 
%			  		lambda = skewness parameter 
%---------------------------------------------------
% RETURNS:
%        a matrix of log-likelihood at each element of x
% --------------------------------------------------
% SEE ALSO: tdis_cdf, tdis_rnd, tdis_inv, tdis_prb, skewtdis_pdf
%---------------------------------------------------
%
%  Andrew Patton
%
%  25 June, 2001
%
% This code was used in: 
%
% Patton, Andrew J., 2004, "On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for
% Asset Allocation", Journal of Financial Econometrics, 2(1), 130-168.

nu = theta(1);
lambda = theta(2);

T = size(x,1);
nu = nu(1)*ones(T,1);			
lambda = lambda(1)*ones(T,1);	
   
c = gamma((nu+1)/2)./(sqrt(pi*(nu-2)).*gamma(nu/2));
a = 4*lambda.*c.*((nu-2)./(nu-1));
b = sqrt(1 + 3*lambda.^2 - a.^2);

logc = gammaln((nu+1)/2) - gammaln(nu/2) - 0.5*log(pi*(nu-2));
logb = 0.5*log(1 + 3*lambda.^2 - a.^2);

find1 = (x<(-a./b));
find2 = (x>=(-a./b));
LL1   = logb + logc - (nu+1)/2.*log(1+1./(nu-2).*((b.*x+a)./(1-lambda)).^2);
LL2   = logb + logc - (nu+1)/2.*log(1+1./(nu-2).*((b.*x+a)./(1+lambda)).^2);
LL    = sum(LL1(find1)) + sum(LL2(find2));
LL    = -LL;