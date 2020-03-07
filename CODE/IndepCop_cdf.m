% Computes the value of the Independence copula cdf at a specified point
% 
% INPUTS:	U, a Tx1 vector (or a scalar) of F(X[t])
%				V, a Tx1 vector (or a scalat) of G(Y[t])
%
% Monday, 2 July, 2001.
%
% Andrew Patton

function cdf = IndepCop_cdf(u,v)


if nargin<2 || isempty(v)
    v = u(:,2);
    u = u(:,1);
end

cdf = (u>=0).*(u<=1).*(v>=0).*(v<=1).*u.*v + (u>=0).*(u<=1).*(v>1).*u + ...
   		(v>=0).*(v<=1).*(u>1).*v + (u>1).*(v>1);