function out1 = Bonferroni_MR_test1(data)
% function out1 = Bonferroni_MR_test1(data)
%
% Function to implement a test of weak monotonicity using Bonferroni bounds
%
%          H0: d1>=0,d2>=0,...,dK>=0
%     vs.  H1: dj<0  for some j=1,2,..,K
%
%          H0*: d1<=0,d2<=0,...,dK<=0
%     vs.  H1*: dj>0  for some j=1,2,..,K
%
%
% Used in the paper:
%
% "Monotonicity in Asset Returns: New Tests with Applications to the Term
% Structure, the CAPM and Portfolio Sorts"
% Paper available at:
% http://www.economics.ox.ac.uk/members/andrew.patton/research.html
%
% INPUTS:   data, a TxK matrix of data (return *differences* for the sorting application)
%
% OUTPUTS:  out1, a 2x1 vector, the Bonferroni p-values from the two tests
%
%  Andrew Patton
%
%  28 May 2008

[T,K] = size(data);

muhat = mean(data)';  % unconstrained estimate of mu
omegahat = newey_west(data, floor(4*((T/100)^(2/9))) )/T; % HAC estimator of the covariance matrix of muhat
tstats = muhat./sqrt(diag(omegahat));

% below are "Bonferroni p-values", in the sense that we reject the null if
% they are less than the size of the test. NOTE of course that unlike usual
% p-vals these won't be Unif(0,1) under the null. In fact, they do not even
% have to lie in [0,1] - they could be bigger than 1.
out1(1) = K*normcdf(min(tstats));  % if this is less than alpha (size of test) then reject H0 in favour of H1
out1(2) = K*normcdf(-max(tstats)); % if this is less than alpha (size of test) then reject H0* in favour of H1*
