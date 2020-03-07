function out1 = tau2kappa(tau);
% Takes in upper and lower tail dependence measures
% and returns the required parameters of the Joe-Clayton copula
%
%	INPUTS:		tau = [tauU, tauL], a 2xT vector of upper and lower tail dependence measures
%
%	OUTPUS:		kappa, a 2xT vector of required Joe-Clayton copula parameters
%
%  Andrew Patton
%
%  4 May, 2002

out1 = nines(size(tau,1),size(tau,2));
out1(:,1) = 1./log2(2-tau(:,1));
out1(:,2) = -1./log2(tau(:,2));
