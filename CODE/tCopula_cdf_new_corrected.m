function out1 = tCopula_cdf_new(U,V,theta,nobs)
% function out1 = tCopula_cdf_new(U,theta,nobs)
%
% Computes the value of the Student's t copula at a specified point
% 
% INPUTS:	U,   a Tx1 vector (or a scalar) of F(X[t])
%				V,   a Tx1 vector (or a scalat) of G(Y[t])
%				RHO, a Tx1 vector (or a scalar) of correlation coefficients
%				NU, 	a Tx1 vector (or a scalar) of degrees of freedom coefficients
%				nobs, a scalar: the number of Monte Carlo samples to draw (actually 4 times this effectively)
%						default: 12,500.
%
% Thursday, 26 April, 2001.
%    Modified 2 April 2011 (just changing style of inputs)
%
% Andrew Patton


% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton

% V = U(:,2);
% U = U(:,1);
RHO = theta(1);
NU = 1/theta(2);  % note that second paraemter is INVERSE degrees of freedom

T = size(U,1);

if nargin<=4;
   nobs = 12500;
end
   
out1 = -999.99*ones(T,1);

% stretching the input vectors to match
RHO = ones(T,1)*RHO;
NU = ones(T,1)*NU;

out1 = -999.99*ones(T,1);
for tt = 1:T;
   if NU(tt)>100;		% then just use Normal copula
      out1(tt) = NormalCopula_cdf_mc([U(tt),V(tt)],RHO(tt),sims);
   else
      x = tinv(U(tt),NU(tt))*sqrt((NU(tt)-2)/NU(tt));  % need to adjust these as the bivartcdfmc.m is for *standardised* t random variables
      y = tinv(V(tt),NU(tt))*sqrt((NU(tt)-2)/NU(tt));
      out1(tt) = bivartcdfmc(x,y,RHO(tt),NU(tt),nobs);
   end
end

