function out1 = kappa2tau(kappa);
% Takes in parameters of the Joe-Clayton
% copula and returns the implied
% upper and lower tail dependence measures
%
%	INPUTS:		kappa, a 2x1 vector of copula parameters
%
%	OUTPUS:		out1 = [tauU, tauL] the implied upper and lower tail dependence measures
%
%  Andrew Patton
%
%  4 May, 2002

out1 = nines(size(kappa,1),size(kappa,2));
out1(1) = 2-2^(1/kappa(1));
out1(2) = 2^(-1/kappa(2));