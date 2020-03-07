function out1 = isposdef(X);
%function out1 = isposdef(X);
%
% Returns a 1 if the matrix X is positive definite
% and a 0 otherwise
%
%  INPUTS:	X, a kxk matrix
%
%  OUTPUTS: out1, a scalar, =1 if X is pos def, =0 else
%
%  Andrew Patton
%
%  Wednesday, 3 July, 2002.

eigs = eig(X);
out1 = (sum(eigs>0)==length(X))*isreal(eigs);  % want all eigenvalues to be positive real numbers, else return zero.

