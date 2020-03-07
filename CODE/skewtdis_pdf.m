function pdf = skewtdis_pdf(x, nu, lambda)
% PURPOSE: returns the pdf at x of Hansen's (1994) 'skewed t' distribution
%---------------------------------------------------
% USAGE: pdf = skewtdis_pdf(x,nu,lambda)
% where: x  = a matrix, vector or scalar 
%        nu = a matrix or scalar degrees of freedom parameter 
%			  lambda = a maxtrix or scalar skewness parameter 
%---------------------------------------------------
% RETURNS:
%        a matrix of pdf at each element of x       
% --------------------------------------------------
% SEE ALSO: tdis_cdf, tdis_rnd, tdis_inv, tdis_prb, skewtdis_cdf
%---------------------------------------------------
%
% Based on tdis_pdf.m from the "Spatial Econometrics"
% toolbox of James P. LeSage
% http://www.spatial-econometrics.com/
%
%  Andrew Patton
%
%  25 June, 2001
%
% This code was used in: 
%
% Patton, Andrew J., 2004, "On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for
% Asset Allocation", Journal of Financial Econometrics, 2(1), 130-168.



T = size(x,1);
if size(nu,1)<T
   nu = nu(1)*ones(T,1);
end
if size(lambda,1)<T
   lambda = lambda(1)*ones(T,1);
end
c = gamma((nu+1)/2)./(sqrt(pi*(nu-2)).*gamma(nu/2));
a = 4*lambda.*c.*((nu-2)./(nu-1));
b = sqrt(1 + 3*lambda.^2 - a.^2);

pdf1 = b.*c.*(1 + 1./(nu-2).*((b.*x+a)./(1-lambda)).^2).^(-(nu+1)/2);
pdf2 = b.*c.*(1 + 1./(nu-2).*((b.*x+a)./(1+lambda)).^2).^(-(nu+1)/2);
pdf  = pdf1.*(x<(-a./b)) + pdf2.*(x>=(-a./b));