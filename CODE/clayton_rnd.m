function [out1,out2] = clayton_rnd(kappa,n,state)
%function out1 = clayton_rnd(kappa,n)
% Generates n (bivariate) random numbers from the Clayton copula
% with parameter kappa.
%
% INPUTS:	kappa, a scalar or an nx1 vector, the parameter of Clayton's copula
%				n , a scalar, the number of draws required
%				state, an integer to use to seed the random number generator
%
% OUTPUTS:	out1, an nx2 matrix of random numbers from the bivariate Clayton copula, OR
%				 	if nargout=2, out1=nx1 and out2=nx1 vectors of individual Unif(0,1) r.v.'s
%
%  Andrew Patton
%
%  11 May, 2001.


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
if size(kappa,1)==1
   kappa = kappa*ones(n,1);		% stretching vector
end
if nargin<3
   rand('state',sum(1234*clock));	% setting RNG to new seed according to computer clock time.
else
   rand('state',state);
end
U = rand(n,1);
t = rand(n,1);		% interim variable
V = (1-U.^(-kappa)+(t.*(U.^(1+kappa))).^(-kappa./(1+kappa))).^(-1./kappa);
out1 = [U,V];
if nargout==2
   out2 = out1(:,2);
   out1 = out1(:,1);
end