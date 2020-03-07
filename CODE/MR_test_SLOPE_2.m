function [out1 out2 ] = MR_test_SLOPE_1(depvars,indepvar,Nindepvar,bootreps,direction,block_length,rand_state);
%
% Function to compute the "monotonic relationship" test on estimated regression coefficients in: 
%
% "Monotonicity in Asset Returns: New Tests with Applications to the Term
% Structure, the CAPM and Portfolio Sorts"
%
% Paper available at:
% http://www.economics.ox.ac.uk/members/andrew.patton/research.html
%
%  THIS FUNCTION TESTS MONOTONICITY OF COEFFICIENTS IN A LINEAR REGRESSION
%  (will test the FIRST coefficient in the regression)
%
%  THIS FUNCTION IS FOR VARIABLES SORTED ON A SINGLE FACTOR
%
%  INPUTS:  depvars, a Txn matrix of returns on the n portfolios (the dependent variables in the regressions)
%           indepvar, a TxQ matrix of independent variables in the regressions 
%           Nindepvar, a nx1 vector of integers, the number of indepvars in each regression (eg: could 
%                   be K*ones(n,1) if all have K regressors, or could have different number per regression.
%                   OR: equals 0 if want to use the same matrix of independent variables for all dependent variables
%           bootreps, a scalar, the number of bootstrap replications to use for the MR test (default=1000)
%           direction, =  1 if want to test for an *increasing* relationship, 
%                      = -1 if want to test for a *decreasing* relationship, 
%                      =  0 (default) if want to use the sign of (muhat_n-muhat_1) to determine the sign to look at
%                      =  2, if data are already in differences (eg: for two-way sorts) and want to test for an *increasing* relationship
%                      = -2, if data are already in differences (eg: for two-way sorts) and want to test for an *decreasing* relationship
%           block_length, a scalar, the average length of the block to use in Politis and Romano's stationary bootstrap (default=10)
%           rand_state, a scalar, the "seed" for Matlab's random number generator (default is sum(100*clock) )
%
%  OUTPUTS: out1, a 5x2 matrix, (1) top minus bottom value, (2) bootstrap "t-test" on top minus bottom, 
%                               (3) the p-value associated with a test that beta1_1=beta1_n (slope equivalent to top-minus-bottom t-test in standard case)
%                               (4) the MR test p-value from the proposed test, (5) the MR test p-value from the proposed test on ALL pair-wise comparisons 
%                               First col on non-studentised betas, second on approximated studentised betas
%           out2, a 1xn vector of the estimated slope coefficients for each portfolio (only the coeff of interest, namely the first one)
%
%  Andrew Patton and Allan Timmermann
%
%  12 April 2009

% NOTE: To simplify, we assume that the independent variables are ordered
% so that the test of monotonicity is applied to the first regression
% coefficient

out1 = -999.99*ones(5,2);
[T,n] = size(depvars);

% creating a cell array to store the indepvars
% (easier for me to use this below, and most Matlab users in finance will
% be more comfortable with big matrices than with cell arrays)
indepvar2 = cell(n,1);
if sum(Nindepvar)>0  % then user has specified the regressors separately for each regression
    counter = 0;
    for ii=1:n;
        indepvar2{ii} = indepvar(:,counter+1:counter+Nindepvar(ii));
        counter = counter + Nindepvar(ii);
    end
elseif Nindepvar==0  % then user wants to use same set of indep vars in all regressions
    for ii=1:n;
        indepvar2{ii} = indepvar;
    end
else 
    'Problem: input ''Nindepvar'' must be a nx1 vector or a scalar equal to zero'
end    


if nargin<4 || isempty(bootreps);
    bootreps=1000;
end
if nargin<5 || isempty(direction)
    direction=0;
end
if nargin<6 || isempty(block_length);
    block_length=10;
end
if nargin<7 || isempty(rand_state)
    rand_state = [];
end
if direction==0;  % then use difference between top and bottom beta1 to determine the direction
    temp1 = ols(depvars(:,1),indepvar2{1});
    tempn = ols(depvars(:,n),indepvar2{n});
    direction = sign(tempn.beta(1)-temp1.beta(1));
end

% generating the time indices for the bootstrapped data:
bootdates = stat_bootstrap_function_21(T,bootreps,block_length,rand_state);     

% estimated beta and studentised beta (ie, t-stat) for actual data
beta0 = -999.99*ones(n,1); 
beta0vcv = -999.99*ones(n,1);
for ii=1:n;
    beta = indepvar2{ii}\depvars(:,ii);  % fast way to obtain ols estimate
    sig2e = cov(depvars(:,ii)-indepvar2{ii}*beta);  % variance of residal
    vcv = sig2e*inv(indepvar2{ii}'*indepvar2{ii});  % standard OLS estimate of covariance matrix of estiamted betas
    beta0(ii) = beta(1);  % test will be based just on the first estimated beta, so only need to store that
    beta0vcv(ii) = vcv(1,1); 
end

out2 = beta0;  % outputting the estimated slope coefficients to users

% this matrix will "stretch the diffs in adjacent portfolios to cover all possible differences
% (used in the MRall test)
R = eye(n-1);
for ii=2:n-1;  % looping through the block sizes (2 up to N-1)
    for jj=1:n-1 - ii+1  % looping through the starting column for blocks
        R = [R;[zeros(1,jj-1),ones(1,ii),zeros(1,n-1-ii-jj+1)]];
    end
end

% now obtaining the bootstrap distributions of the betas
temp   = -999.99*ones(bootreps,n);
tempA  = -999.99*ones(bootreps,n*(n-1)/2);
for bb=1:bootreps;
    for ii=1:n;
        temp2 = indepvar2{ii};
        beta = temp2(bootdates(:,bb),:)\depvars(bootdates(:,bb),ii);  % fast way to obtain ols estimate
        sig2e = cov(depvars(bootdates(:,bb),ii)-temp2(bootdates(:,bb),:)*beta);  % variance of residal
        vcv = sig2e*inv(temp2(bootdates(:,bb),:)'*temp2(bootdates(:,bb),:));  % standard OLS estimate of covariance matrix of estiamted betas
        temp(bb,ii) = beta(1)-beta0(ii);  % first estimated beta, re-centered using actual estiamte of first beta
    end
    tempA(bb,:) = R*( direction*( temp(bb,2:n)-temp(bb,1:n-1) )' );  % stretching the adjacent differences to all possible differences
end


% test 1: testing whether betaN is different from beta1
Jstat1 = direction*(beta0(n)-beta0(1));
out1(1,1) = direction*(beta0(n)-beta0(1));
out1(2,1) = direction*(beta0(n)-beta0(1))/std(direction*(temp(:,n)-temp(:,1)));
out1(3,1) = mean( direction*(temp(:,n)-temp(:,1)) > Jstat1 );  % bootstrap p-value on top minus bottom beta
out1(1:3,2) = out1(1:3,1);  % for a simple pair-wise test, studentisation does not change things

% test 2: testing all adjacent inequality constraints (the MR test)
Jstat2 = min( direction*( beta0(2:n)-beta0(1:n-1) ) );
out1(4,1) = mean( min( direction*( temp(:,2:n)-temp(:,1:n-1) ) ,[], 2)>Jstat2 );  % bootstrap p-value for the MR test 

Jstat2s = min( direction*( (beta0(2:n)-beta0(1:n-1))./sqrt( beta0vcv(2:n)+beta0vcv(1:n-1) ) ) );  % studentised MR test (using naive studentisation: not HAC variance esitmate, and ignoring cross-sectional correlation in beta1hat)
out1(4,2) = mean( min( direction*( (temp(:,2:n)-temp(:,1:n-1)))./(ones(bootreps,1)*sqrt( beta0vcv(2:n)+beta0vcv(1:n-1) )')  ,[],2) > Jstat2s );   % bootstrap p-value for the studentised MR test 

% test 3: testing all possible inequality constraints (the MRall test)
Jstat3 = min( R*(direction*( beta0(2:n)-beta0(1:n-1) ) ) );
out1(5,1) = mean( min( tempA ,[], 2)> Jstat3 );  % bootstrap p-value for the MRall test 

beta0vcvA = R*diag(beta0vcv(2:n)+beta0vcv(1:n-1))*R';   % covariance matrix for all possible comparisons, assuming zero cross-sectional correlation
Jstat3s = min( (R*(direction*( beta0(2:n)-beta0(1:n-1) )))./sqrt(diag(beta0vcvA)) );  
out1(5,2) = mean( min( tempA./(ones(bootreps,1)*sqrt(diag(beta0vcvA))') ,[], 2)> Jstat3s );  % bootstrap p-value for the studentised MRall test 
