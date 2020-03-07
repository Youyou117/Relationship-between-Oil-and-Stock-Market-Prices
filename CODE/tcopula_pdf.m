function out1 = tcopula_pdf(U,V,RHO,NU)
% function out1 = tcopula_pdf(U,V,RHO,NU)
%
% Computes the value of the bivariate Student's t copula pdf at a specified point
% 
% INPUTS:	U,  a Tx1 vector or a scalar
%				V,  a Tx1 vector or a scalar
%				RHO, a Tx1 vector or a scaler, correlation coefficient
%				NU, 	a Tx1 vector or a scalar, the degrees of freedom parameter
%
% Sunday, 29 July, 2001.
%
% Andrew Patton

% Written for the following papers:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton


X1 = tdis_inv(U,NU);
X2 = tdis_inv(V,NU);

T = max(size(X1,1),max(size(X2,1),max(size(RHO,1),size(NU,1))));
if size(X1,1)~=T;
   X1 = ones(T,1)*X1(1);
end
if size(X2,1)~=T;
   X2 = ones(T,1)*X2(1);
end
if size(RHO,1)~=T;
   RHO = ones(T,1)*RHO(1);
end
if size(NU,1)~=T;
   NU = ones(T,1)*NU(1);
end
out1 = -999.99*ones(T,1);
for tt = 1:T;
   out1(tt) = gamma((NU(tt)+2)/2)/gamma(NU(tt)/2)*((NU(tt)*pi)^(-1))/sqrt(1-RHO(tt)^2)*((1 + (X1(tt)^2+X2(tt)^2-2*RHO(tt)*X1(tt)*X2(tt))/(NU(tt)*(1-RHO(tt)^2)))^(-(NU(tt)+2)/2));
   out1(tt) = out1(tt)/(tdis_pdf(X1(tt),NU(tt))*tdis_pdf(X2(tt),NU(tt)));
end
   