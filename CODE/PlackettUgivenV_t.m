% Computes the value of the conditional Plackett (U given V)
% copula cdf at a specified point
% 
% INPUTS:	u, a scalar of F(X[t])
%				v, a scalar of G(Y[t])
%				t, the value of C(u|v)
%				k, a scalar of kappa
%
% Monday, 13 August, 2001
%
% Andrew Patton

function out1 = PlackettUgivenV_t(u,v,t,k)

eta = k-1;
out1 = 0.5 - 0.5*(eta.*v + 1 - (eta+2).*u)./sqrt((1 + eta.*(u+v)).^2 - 4*k.*eta.*u.*v);
out1 = out1 - t;
