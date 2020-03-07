function out1 = bivarnormcdf_arg(X,b2,rho);
% Used in 'bivarnormcdf.m'
%
% INPUTS:	rho, a scalar, the correlation coefficient of the joint density
%				X, a Tx1 vector, the value of the first (outer) variable
%				b1, a scalar, the upper bound of the inner integral
%
% OUTPUTS: out1, the value of the argument at the given co-ordinate
%
%  Andrew Patton
%
% Thursday, 14 December, 2000


out1 = normcdf(b2,rho*X,sqrt(1-rho^2)).*normpdf(X,0,1);

