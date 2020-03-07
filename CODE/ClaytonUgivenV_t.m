function out1 = ClaytonUgivenV_t(u,v,t,k)
% function out1 = ClaytonUgivenV_t(u,v,t,k)
%
% Computes the value of the conditional Clayton (U given V)
% copula cdf at a specified point
% 
% INPUTS:	u, a scalar of F(X[t])
%				v, a scalar of G(Y[t])
%				t, the value of C(u|v)
%				k, a scalar of kappa
%
% 17 Oct 2011 (based on "GumbelUgivenV_t.m")
%
% Andrew Patton


if nargin<4 || isempty(v)  % then data was given as a Tx2 matrix, so need to re-label the inputs
    k=t;
    t=v;
    v = u(:,2);
    u = u(:,1);
end

out1 = (v.^(-1-k)) .* ( ( -1 + (u.^(-k)) + (v.^(-k)) ) .^ (-1-1/k) );  % from Mathematica

out1 = out1 - t;
