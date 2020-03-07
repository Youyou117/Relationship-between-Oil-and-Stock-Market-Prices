function out1 = MR_test_22(data,bootreps,direction,block_length,rand_state);
%
% Function to compute the "monotonic relationship" tests in: 
%
% "Monotonicity in Asset Returns: New Tests with Applications to the Term
% Structure, the CAPM and Portfolio Sorts"
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
%  OUTPUTS: out1, a 4x2 vector, (1) the t-statistic associated with a t-test that mu_1=mu_n, 
%                               (2) the p-value associated with a t-test that mu_1=mu_n, 
%                               (3) the MR test p-value from the proposed test (non-studentised)
%                               (4) the MR test p-value from the proposed test, on ALL pair-wise comparisons (non-studentised)
%                   first column is non-studentised, second column is studentised (NOTE: tstats and pvals are not affected by studentisation, so are the same)
%
%  Andrew Patton and Allan Timmermann
%
%  4 May 2008

% This paper is a revised version of the paper previously circulated as:
% 
%  "Portfolio Sorts and Tests of Cross-Sectional Patterns in Expected Returns"
%    by Andrew J. Patton and Allan Timmermann


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
temp = min(temp,[],2);          % looking at the minimum difference in portfolio means, for each of the bootstrapped data sets
Jstat = min(dmuhat);            % the test statistic
out1(3,1) = mean(temp>Jstat);     % the p-value associated with the test statistic

tempS = min(tempS,[],2);                            % looking at the minimum STUDENTISED difference in portfolio means, for each of the bootstrapped data sets
JstatS = min(dmuhat./sqrt(diag(SBvariance)));       % the STUDENTISED test statistic
out1(3,2) = mean(tempS>JstatS);                       % the p-value associated with the STUDENTISED test statistic

if abs(direction)<2  % if data are already in differences, then I cannot do the t-test. this must be done on the actual returns
    % now getting the t-statistic and p-value from the usual t-test
    temp = nwest(direction*(data(:,end)-data(:,1)), ones(T,1), floor(4*((T/100)^(2/9))));  % truncation lag for Newey-West standard errors, see Newey and West (1987, 1994)
    out1(1,:) = temp.tstat*ones(1,2);
    out1(2,:) = normcdf(-temp.tstat)*ones(1,2);
end

% 17apr08: running the test on all possible differences, rather than just the adjacent portfolios
% this matrix will "stretch the diffs in adjacent portfolios to cover all
% possible differences
R = eye(n-1);
for ii=2:n-1;  % looping through the block sizes (2 up to N-1)
    for jj=1:n-1 - ii+1  % looping through the starting column for blocks
        R = [R;[zeros(1,jj-1),ones(1,ii),zeros(1,n-1-ii-jj+1)]];
    end
end

diffdata2 = diffdata*R';
dmuhat2 = mean(diffdata2)';
temp = -999.99*ones(bootreps,n*(n-1)/2);
tempS = -999.99*ones(bootreps,n*(n-1)/2);
SBvariance = SB_variance(diffdata2,block_length);  
for ii=1:n*(n-1)/2;
    temp2 = diffdata2(:,ii);
    temp(:,ii) = mean(temp2(bootdates))'-dmuhat2(ii);      % the mean of each of the bootstrap shuffles of the original data, minus the mean of the original data (the re-centering bit)
    tempS(:,ii) = (mean(temp2(bootdates))'-dmuhat2(ii))/sqrt(SBvariance(ii,ii));  % studentising the difference using the SB estimate of its std deviation
end
temp = min(temp,[],2);          % looking at the minimum difference in portfolio means, for each of the bootstrapped data sets
Jstat = min(dmuhat2);            % the test statistic
out1(4,1) = mean(temp>Jstat);     % the p-value associated with the test statistic

tempS = min(tempS,[],2);                            % looking at the minimum STUDENTISED difference in portfolio means, for each of the bootstrapped data sets
JstatS = min(dmuhat2./sqrt(diag(SBvariance)));       % the STUDENTISED test statistic
out1(4,2) = mean(tempS>JstatS);                       % the p-value associated with the STUDENTISED test statistic

