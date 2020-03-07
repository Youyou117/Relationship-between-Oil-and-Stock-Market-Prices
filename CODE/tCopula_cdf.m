function out1 = tCopula(U,V,RHO,NU,nobs)
% function out1 = tCopula(U,V,RHO,NU,nobs)
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
%
% Andrew Patton


% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton


T = max([size(U,1),size(V,1),size(RHO,1),size(NU,1)]);

if nargin<=4;
   nobs = 12500;
end
   
out1 = -999.99*ones(T,1);

% stretching the input vectors to match
if size(U,1)<T;
   U = ones(T,1)*U(1);
end
if size(V,1)<T;
   V = ones(T,1)*V(1);
end
if size(RHO,1)<T;
   RHO = ones(T,1)*RHO(1);
end
if size(NU,1)<T;
   NU = ones(T,1)*NU(1);
end

out1 = -999.99*ones(T,1);
for tt = 1:T;
   if NU(tt)>100;		% then just use Normal copula
      out1(tt) = NormalCopula(U(tt),V(tt),RHO(tt));
   else
      x = tdis_inv(U(tt),NU(tt))*sqrt((NU(tt)-2)/NU(tt));  % need to adjust these as the bivartcdfmc.m is for *standardised* t random variables
      y = tdis_inv(V(tt),NU(tt))*sqrt((NU(tt)-2)/NU(tt));
      out1(tt) = bivartcdfmc(x,y,RHO(tt),NU(tt),nobs);
   end
end

