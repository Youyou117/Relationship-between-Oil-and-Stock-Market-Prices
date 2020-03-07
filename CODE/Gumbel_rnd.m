function out1 = Gumbel_rnd(kappa,T,state)
% function out1 = Gumbel_rnd(kappa,T,state)
%
% This program generates observations from
% a Gumbel copula
%
%  INPUTS:	kappa, a scalar or a Tx1 vector, the copula parameter
%					T,	a scalar, the number of obs in each sample
%					state, an integer to use to seed the random number generator
%
%  OUTPUT:	out1, a Tx2 matrix of data.
%
% Thursday, 8 November, 2001
%
%	Andrew Patton

% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton


if size(kappa,1)~=T;		% stretching the input vectors to match each other
   kappa = kappa*ones(T,1);
end
if nargin<3
   rand('state',sum(1234*clock));	% setting RNG to new seed according to computer clock time.
else
   rand('state',state);
end


out1 = -999.99*ones(T,2);
out1(1:T,2) = rand(T,1);			% transforms of Y
temp = rand(T,1);							% interim variable
for jj = 1:T;
   out1(jj,1) = GumbelUgivenV_inverse2(out1(jj,2),temp(jj),kappa(jj));
end
