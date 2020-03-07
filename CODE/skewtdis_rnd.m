function out1 = skewtdis_rnd(nu,lambda,T,state)
% PURPOSE: returns random draws from Hansen's (1994) 'skewed t' distribution
%---------------------------------------------------
% USAGE: out1 = skewtdis_rnd(nu,lambda,T)
% where:  nu = a matrix or scalar degrees of freedom parameter 
%			  lambda = a maxtrix or scalar skewness parameter 
%				T = number of draws
%				state, an integer to use to seed the random number generator
%---------------------------------------------------
% RETURNS:
%        a Tx1 vector of draws from the dist'n  
% --------------------------------------------------
% SEE ALSO: tdis_cdf, tdis_rnd, tdis_inv, tdis_prb
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



if nargin<4
   rand('state',sum(1234*clock));	% setting RNG to new seed according to computer clock time.
else
   rand('state',state);
end
if size(nu,1)<T
   nu = nu(1)*ones(T,1);
end
if size(lambda,1)<T
   lambda = lambda(1)*ones(T,1);
end
u = rand(T,1);
out1 = skewtdis_inv(u,nu,lambda);
