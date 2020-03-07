function out1 = empirical_copula(data,points)
% function out1 = empirical_copula(data,points)
%
%  Function to compute the empirical copula (CDF) using a given data set, and then evaluate the empirical copula at the set of given points.
%
%  INPUTS:  data, a Txk matrix of Unif(0,1) data to use to estimate the copula CDF
%           points, a Sxk matrix of values inside [0,1]^k at which to evaluate the empirical copula. (If empty, will set points=data)
%
%  OUTPUTS: out1, a Tx1 vector, the value of the empirical copula evaluated at the matrix "points"
%
%  Andrew Patton
%
%  30 September 2011

[T,k] = size(data);

if nargin<2 || isempty(points)
    points = data;
end

S = size(points,1);

out1 = nan(S,1);

for ss=1:S;
    out1(ss) = 1/(T+1)*sum( prod(1*(data<=(ones(T,1)*points(ss,:))),2) );  % note that I multiply by "1" here as the function "prod" does not like logical inputs (0/1). Multiplying by 1 converts to a number.
end
