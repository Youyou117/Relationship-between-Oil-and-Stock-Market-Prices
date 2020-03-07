function out1 = MR_test_22_simple(data,bootreps,direction,block_length,rand_state);
%
% Function to compute the "monotonic relationship" tests in: 
%
% "Monotonicity in Asset Returns: New Tests with Applications to the Term Structure, the CAPM and Portfolio Sorts"
% by Andrew J. Patton and Allan Timmermann, forthcoming in the Journal of Financial Economics, 2010.
%
% Paper available at:
% http://www.economics.ox.ac.uk/members/andrew.patton/research.html
%
%  THIS FUNCTION IS FOR VARIABLES SORTED ON A SINGLE FACTOR
%
%  INPUTS:  data, a Txn matrix of returns on the n portfolios
%           bootreps, a scalar, the number of bootstrap replications to use for the MR test (default=1000)
%           direction, =  1 if want to test for an *increasing* relationship, 
%                      = -1 if want to test for a *decreasing* relationship, 
%                      =  0 (default) if want to use the sign of (muhat_n-muhat_1) to determine the sign to look at
%                      =  2, if data are already in differences (eg: for two-way sorts) and want to test for an *increasing* relationship
%                      = -2, if data are already in differences (eg: for two-way sorts) and want to test for an *decreasing* relationship
%           block_length, a scalar, the average length of the block to use in Politis and Romano's stationary bootstrap (default=10)
%           rand_state, a scalar, the "seed" for Matlab's random number generator (default is sum(100*clock) )
%
%  OUTPUTS: out1, a scalar, the MR test p-value from the proposed test
%
%  Andrew Patton and Allan Timmermann
%
%  4 May 2008


out1 = -999.99*ones(4,2);

[T,n] = size(data);

if nargin<2 | isempty(bootreps);
    bootreps=1000;
end
if nargin<3 | isempty(direction)
    direction=0;
end
if nargin<4 | isempty(block_length);
    block_length=10;
end
if nargin<5 | isempty(rand_state)
    rand_state = [];
end
if direction==0;
    direction = sign(mean(data(:,n)-data(:,1)));
end

if abs(direction)<2  % then need to difference the data
    diffdata = data(:,2:end)-data(:,1:end-1);   % differences across the columns of the data
    diffdata = direction*diffdata;              % changing the sign of these if want to look for a *decreasing* pattern rather than an *increasing* pattern
else
    diffdata = data;  % data are already in differences
    diffdata = sgn(direction)*diffdata;  % changing the sign of these if want to look for a *decreasing* pattern rather than an *increasing* pattern
    n = n + 1;      % if data are already differenced, then this is like having one more column in the raw data
end
dmuhat = mean(diffdata)';

% generating the time indices for the bootstrapped data:
bootdates = stat_bootstrap_function_21(T,bootreps,block_length,rand_state);     

% this is the long-run variance of diffdata, according to the stationary bootstrap (obtained analytically using Lemma 1 of PR(94)
% I only have to compute this once, as it does not depend on the bootstrap random draws
SBvariance = SB_variance(diffdata,block_length);  


temp = -999.99*ones(bootreps,n-1);
tempS = -999.99*ones(bootreps,n-1);
for ii=1:n-1;
    temp2 = diffdata(:,ii);
    temp(:,ii) = mean(temp2(bootdates))'-dmuhat(ii);  % the mean of each of the bootstrap shuffles of the original data, minus the mean of the original data (the re-centering bit)
    tempS(:,ii) = (mean(temp2(bootdates))'-dmuhat(ii))/sqrt(SBvariance(ii,ii));  % studentising the difference using the SB estimate of its std deviation
end
tempS = min(tempS,[],2);                            % looking at the minimum STUDENTISED difference in portfolio means, for each of the bootstrapped data sets
JstatS = min(dmuhat./sqrt(diag(SBvariance)));       % the STUDENTISED test statistic
out1 = mean(tempS>JstatS);                       % the p-value associated with the STUDENTISED test statistic
