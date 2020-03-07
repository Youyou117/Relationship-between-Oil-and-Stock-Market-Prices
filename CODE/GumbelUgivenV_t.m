function out1 = GumbelUgivenV_t(u,v,t,k)
% function out1 = GumbelUgivenV_t(u,v,t,k)
%
% Computes the value of the conditional Gumbel (U given V)
% copula cdf at a specified point
% 
% INPUTS:	u, a scalar of F(X[t])
%				v, a scalar of G(Y[t])
%				t, the value of C(u|v)
%				k, a scalar of kappa
%
% Saturday, 28 July, 2001
%
% Andrew Patton


ut = -log(u);
vt = -log(v);

out1 = (v.^(-1)).*exp(-(ut.^k+vt.^k).^(1/k)).*(1+(ut./vt).^k).^(1/k-1);
out1 = out1 - t;
