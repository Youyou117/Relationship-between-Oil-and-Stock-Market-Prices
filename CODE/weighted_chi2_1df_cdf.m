function out1 = weighted_chi2_1df_cdf(x,weights,sims)
% function out1 = weighted_chi2_1df_cdf(x,weights,sims)
%
%  Function to compute the CDF of a variable equal to a weighted sum of
%  independent chi-squared variables each with one degree of freedom.
%
%   Using simulations to obtain this distribution
%
% INPUTS:   x, a px1 vector of values at which to evaluate the CDF
%           weights, a kx1 vector of weights to attach to each chi-squared variable
%
% OUPUTS:   out1, a px1 vector, the value of the CDF at the specified points
%
%  Andrew Patton
%
%  4 April 2011


if nargin<3 || isempty(sims)
    sims = 10000;
end

K = length(x);

simdata = chi2rnd(1,sims,length(weights));  % K independent chi2 variables with 1 DoF
simdata = simdata*weights;    % the simulated weighted average variable

out1 = nan(K,1);
for ii=1:K;
    out1(ii) = mean(simdata<=x(ii));
end
