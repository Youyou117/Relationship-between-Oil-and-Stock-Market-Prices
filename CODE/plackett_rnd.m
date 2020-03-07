function [out1,t] = plackett_rnd(kappa,n,state)
%function out1 = plackett_rnd(kappa,n)
% Generates n (bivariate) random numbers from the Plackett copula
% with parameter kappa.
%
% INPUTS:	kappa, a scalar, the parameter of the Plackett copula, must be greater than 0 (correl=0.5 at about kappa=5.5)
%				n , a scalar, the number of draws required
%				state, an integer to use to seed the random number generator
%
% OUTPUTS: out1, an nx2 matrix of random numbers 
%
%  Andrew Patton
%
%  Saturday, 12 October, 2002.
%
% Following Nelsen (1998), page 83, Exercise 3.37.

% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton


if nargin==1
   n = 1;	
end
if nargin<3
   rand('state',sum(1234*clock));	% setting RNG to new seed according to computer clock time.
else
   rand('state',state);
end

U = rand(n,1);
t = rand(n,1);		% interim variable

a = t.*(1-t);
b = kappa + a.*((kappa-1)^2);
c = 2.*a.*(U.*(kappa^2)+1-U)+kappa*(1-2*a);
d = sqrt(kappa)*sqrt(kappa + 4*a.*U.*(1-U).*((1-kappa)^2));
V = (c-(1-2*t).*d)./(2*b);
out1 = [U,V];