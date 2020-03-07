function pval = AR_test_rank_correl(data,ARp,bootreps)
% function pval = AR_test_rank_correl(data,ARp,bootreps)a
%
%  Function to test for a time-varying copula through an AR(p) model for rank correl
%
%  * Constraining AR coeffs to be the same (so essentially predicting
%  future dep measure using a rolloing avg of past measures)
%
% INPUTS:   data, a Tx2 matrix of (approx) Unif(0,1) data
%           ARp, a scalar, the order of the AR model to try (default=5)
%           bootreps, a scalar, the number of bootstrap replications (default=1000)
%
% OUTPUTS:  pval, a scalar, the p-value for the test of time-varying rank correlation
%
%  Andrew Patton
%
%  7 April 2011

T = size(data,1);

if nargin<2 || isempty(ARp)
    ARp=5;
end
if nargin<3 || isempty(bootreps)
    bootreps=1000;
end

UV = data(:,1).*data(:,2);
X = [ones(T,1),mean(mlag(UV,ARp),2)];
X = X(ARp+1:end,:);
Y = UV(ARp+1:end);
xpxi = (X'*X)\eye(1+1);
beta = xpxi*(X'*Y);
e = Y-X*beta;
vcv = (e'*e)/(T-ARp-1)*xpxi;
teststat = beta(2:end)'*(vcv(2:end,2:end)\beta(2:end));

bootdates = randi(T,T,bootreps);  % these are the "dates" to use in the iid boostrap
bootstats = nan(bootreps,1);
for bb = 1:bootreps;
    UV = data(bootdates(:,bb),1).*data(bootdates(:,bb),2);
    X = [ones(T,1),mean(mlag(UV,ARp),2)];
    X = X(ARp+1:end,:);
    Y = UV(ARp+1:end);
    xpxi = (X'*X)\eye(1+1);
    beta = xpxi*(X'*Y);
    e = Y-X*beta;
    vcv = (e'*e)/(T-ARp-1)*xpxi;
    bootstats(bb,1) = beta(2:end)'*(vcv(2:end,2:end)\beta(2:end));
end
pval = mean(bootstats>teststat);  % proportion of simulations that generated a test stat bigger than the one observed in the data
