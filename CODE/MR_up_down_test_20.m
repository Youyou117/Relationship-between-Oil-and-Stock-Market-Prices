function out1 = MR_up_down_test_20(data,bootreps,block_length,rand_state)
%
% Function to compute the up and down tests based on:
%  (1) sum of squared differences for positive diffs and negative diffs
%  (2) sum of absolute differences for positive diffs and negative diffs
%  
% From the paper:
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
%           block_length, a scalar, the average length of the block to use in Politis and Romano's stationary bootstrap (default=10)
%           rand_state, a scalar, the "seed" for Matlab's random number generator (default is sum(100*clock) )
%
%  OUTPUTS: out1, a 4x2 vector, the bootstrap p-values from a test for a  monotonic *increaasing* relationship, and a monotonic *decreasing* relationship. 
%               First two rows for squared diffs, second two rows for abs diffs
%               First col not studentised, second col studentised
%
%  Andrew Patton and Allan Timmermann
%
%  6 April 2009

out1 = -999.99*ones(6,2);

[T,n] = size(data);

if nargin<2 || isempty(bootreps);
    bootreps=1000;
end
if nargin<3 || isempty(block_length);
    block_length=10;
end
if nargin<4 || isempty(rand_state)
    rand_state = [];
end

diffdata = data(:,2:end)-data(:,1:end-1);   % differences across the columns of the data
dmuhat = mean(diffdata)';

indic_pos = (dmuhat>0);
indic_neg = (dmuhat<=0);

Jstatpos2 = sum((dmuhat.^2).*indic_pos);   % test stat for an increasing relationship on real data
Jstatneg2 = sum((dmuhat.^2).*indic_neg);   % test stat for an decreasing relationship on real data
Jstatpos1 = sum(abs(dmuhat).*indic_pos);   % test stat for an increasing relationship on real data
Jstatneg1 = sum(abs(dmuhat).*indic_neg);   % test stat for an decreasing relationship on real data

% now studentising the differences
SBvariance = SB_variance(diffdata,block_length);    % stat boot covariance matrix of all differences.
diffdatastd = diffdata./(ones(T,1)*sqrt(diag(SBvariance)'));  % studentised differences
dmuhatstd = mean(diffdatastd)';

Jstatpos2std = sum((dmuhatstd.^2).*indic_pos);   % test stat for an increasing relationship on real data
Jstatneg2std = sum((dmuhatstd.^2).*indic_neg);   % test stat for an decreasing relationship on real data
Jstatpos1std = sum(abs(dmuhatstd).*indic_pos);   % test stat for an increasing relationship on real data
Jstatneg1std = sum(abs(dmuhatstd).*indic_neg);   % test stat for an decreasing relationship on real data

JstatsALL = [[Jstatpos2;Jstatneg2;Jstatpos1;Jstatneg1],[Jstatpos2std;Jstatneg2std;Jstatpos1std;Jstatneg1std]];

% generating the time indices for the bootstrapped data:
bootdates = stat_bootstrap_function_21(T,bootreps,block_length,rand_state);     

% below I loop across columns of diffdata (which looks weird) rather than
% looping across bootstrap samples (which looks easier) for speed.
temp = zeros(bootreps,4);
tempS = zeros(bootreps,4);
for ii=1:n-1;
    temp2 = diffdata(:,ii);
    dmuhatB = mean(temp2(bootdates))' - dmuhat(ii);  % computing mean using bootstrapped data, then re-centering using actual means
    
    temp(:,1) = temp(:,1) + (dmuhatB.^2).*(dmuhatB>0);
    temp(:,2) = temp(:,2) + (dmuhatB.^2).*(dmuhatB<0);
    temp(:,3) = temp(:,3) + abs(dmuhatB).*(dmuhatB>0);
    temp(:,4) = temp(:,4) + abs(dmuhatB).*(dmuhatB<0);

    temp2 = diffdatastd(:,ii);
    dmuhatstdB = mean(temp2(bootdates))' - dmuhatstd(ii);    % computing mean using bootstrapped data, then re-centering using actual means
    tempS(:,1) = tempS(:,1) + (dmuhatstdB.^2).*(dmuhatstdB>0);
    tempS(:,2) = tempS(:,2) + (dmuhatstdB.^2).*(dmuhatstdB<0);
    tempS(:,3) = tempS(:,3) + abs(dmuhatstdB).*(dmuhatstdB>0);
    tempS(:,4) = tempS(:,4) + abs(dmuhatstdB).*(dmuhatstdB<0);
end

out1(1,1) = mean(temp(:,1)>Jstatpos2);
out1(2,1) = mean(temp(:,2)>Jstatneg2);
out1(3,1) = mean(temp(:,3)>Jstatpos1);
out1(4,1) = mean(temp(:,4)>Jstatneg1);

out1(1,2) = mean(tempS(:,1)>Jstatpos2std);
out1(2,2) = mean(tempS(:,2)>Jstatneg2std);
out1(3,2) = mean(tempS(:,3)>Jstatpos1std);
out1(4,2) = mean(tempS(:,4)>Jstatneg1std);



