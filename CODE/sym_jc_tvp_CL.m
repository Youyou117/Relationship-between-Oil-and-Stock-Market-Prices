function [CL, TAU1, TAU2] = sym_jc_tvp_CL(theta,Z,thetabar)
% function [CL, TAU1, TAU2] = sym_jc_tvp_CL(theta,Z,thetabar)
%
% The negative copula log-likelihood of the symmetrised Joe-Clayton copula 
% with time-varying tail dependence
%
% Wed, 10 May, 2003.
%
% Andrew Patton
%
% INPUTS: theta = [kappa,gamma]
%				Z = [U V];
%				thetabar = [kappabar,gammabar], results from estimation of unconditional model

% Published in:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence,
% International Economic Review, 47(2), 527-556. 
%
% See also:
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of
% Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
%
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and
% Asymmetric Dependence for Asset Allocation, Journal of Financial
% Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton


w1 = theta(1);
a1 = theta(2);
b1 = theta(3);
w2 = theta(4);
a2 = theta(5);
b2 = theta(6);

u = Z(:,1);
v = Z(:,2);
T  = size(Z,1);

TAU1    = -999.99*ones(T,1);
TAU2    = -999.99*ones(T,1);
TAU1(1) = thetabar(1);		% unconditional tau1 and tau2, based on time-invariant version of this model
TAU2(1) = thetabar(2);
psi(1,:)=zeros(1,2);
for jj = 2:T;
   if jj<=10
      psi1 = w1 + b1*TAU1(jj-1) + a1*(mean(abs(u(1:jj-1)-v(1:jj-1))));
      psi2 = w2 + b2*TAU2(jj-1) + a2*(mean(abs(u(1:jj-1)-v(1:jj-1))));
   else
      psi1 = w1 + b1*TAU1(jj-1) + a1*(mean(abs(u(jj-10:jj-1)-v(jj-10:jj-1))));
      psi2 = w2 + b2*TAU2(jj-1) + a2*(mean(abs(u(jj-10:jj-1)-v(jj-10:jj-1))));
   end
   psi(jj,:) = [psi1,psi2];
   TAU1(jj) = 0.998/(1+exp(-psi1)) + 0.001;		% tail dependence parameters
   TAU2(jj) = 0.998/(1+exp(-psi2)) + 0.001;
end

CL = sym_jc_pdf(u,v,TAU1,TAU2);

CL = log(CL);
CL = sum(CL);
CL = -CL;

if isnan(CL)
   CL = 1e6;
end
if isreal(CL)==0
   CL = 1e7;
end
if isreal(theta)==0
   CL = 1e8;
end