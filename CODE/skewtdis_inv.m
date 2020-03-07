function inv = skewtdis_inv(u, nu, lambda)
% PURPOSE: returns the inverse cdf at u of Hansen's (1994) 'skewed t' distribution
%---------------------------------------------------
% USAGE: inv = skewtdis_inv(u,nu,lambda)
% where: U  = a matrix, vector or scalar in the unit interval
%        nu = a matrix or scalar degrees of freedom parameter 
%			  lambda = a maxtrix or scalar skewness parameter 
%---------------------------------------------------
% RETURNS:
%        a matrix of the inverse cdf at each element of u 
% --------------------------------------------------
% SEE ALSO: tdis_cdf, tdis_rnd, tdis_inv, tdis_prb, skewtdis_pdf
%---------------------------------------------------
%
% Based on tdis_inv.m from the "Spatial Econometrics"
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


T = size(u,1);
if size(nu,1)<T
   nu = nu(1)*ones(T,1);
end
if size(lambda,1)<T
   lambda = lambda(1)*ones(T,1);
end
c = gamma((nu+1)/2)./(sqrt(pi*(nu-2)).*gamma(nu/2));
a = 4*lambda.*c.*((nu-2)./(nu-1));
b = sqrt(1 + 3*lambda.^2 - a.^2);

f1 = find(u<(1-lambda)/2);
f2 = find(u>=(1-lambda)/2);

inv1 = (1-lambda(f1))./b(f1).*sqrt((nu(f1)-2)./nu(f1)).*tdis_inv(u(f1)./(1-lambda(f1)),nu(f1))-a(f1)./b(f1);
inv2 = (1+lambda(f2))./b(f2).*sqrt((nu(f2)-2)./nu(f2)).*tdis_inv(0.5+1./(1+lambda(f2)).*(u(f2)-(1-lambda(f2))./2),nu(f2))-a(f2)./b(f2);
inv = -999.99*ones(T,1);
inv(f1) = inv1;
inv(f2) = inv2;