function out1 = weighted_chi2cdf(x,weights);
% function out1 = weighted_chi2cdf(x,weights);
%
%  Function to compute a weighted sum of chi-squared cdf's
%
%  I assume that the sum involves chi-squared variables with degrees of
%  freedom ranging from zero to K
%
% INPUTS:   x, a px1 vector of values at which to evaluate the CDF
%           weights, a (K+1)x1 vector of weights to attach to each chi-squared CDF
%
% OUPUTS:   out1, a px1 vector, the value of the CDF at the specified points
%
%  Andrew Patton
%
%  4 May 2008

K = length(weights)-1;

out1 = weights(0+1)*(x>=0);  % a chi-squared with zero degrees of freedom equals 0 with probability 1. so its CDF is zero until 0, and one for values greater or equal than 0
for ii=1:K;
    out1 = out1 + weights(ii+1)*chi2cdf(x,ii);
end
