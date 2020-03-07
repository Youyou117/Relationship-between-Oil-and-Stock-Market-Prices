function r = ranks(x)
%function r = ranks(x)
% Finds the ranks of each of the arguments of the 
% input matrix
%
% r = ranks(x);
%
% INPUTS: x, a Txk matrix of data, with each
%  variable in a different column 
%
% OUTPUT: r, a Txk matrix of the ranks of
%  each observation for each variable
%
% Monday, 25 Sep, 2000
%
% Andrew Patton
%
% Modified 13dec08


[T,k] = size(x);
r = -999.99*ones(T,k);
for jj=1:k;
    temp = sortrows([(1:T)',sortrows([(1:T)',x(:,jj)],2)],2);
    r(:,jj) = temp(:,1);
end
