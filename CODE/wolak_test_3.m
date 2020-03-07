function out1 = wolak_test_3(data,reps)
% function out1 = wolak_test_2(data)
%
% Function to implement Wolak's (1989, JoE) test of
%
%       H0*: d1>=0,d2>=0,...,dK>=0
%
%  vs.  H1*: (d1,d2,...,dK) in R^K  (ie: general alternative)
%
% *and*
%
%       H0**: d1=d2=...=dK=0
%
%  vs.  H1**: d1>0, d2>0, ..., dK>0
%
% 
% Used in the paper:
%
% "Monotonicity in Asset Returns: New Tests with Applications to the Term
% Structure, the CAPM and Portfolio Sorts"
% Paper available at:
% http://www.economics.ox.ac.uk/members/andrew.patton/research.html
% 
%
% INPUTS:   data, a TxK matrix of data (return *differences* for the sorting application)
%           reps, a scalar, the number of simulations to use to estimate the weight function in the weighted-sum of chi2 variables (default=100)
%
% OUTPUTS:  out1, a 2x1 vector, the p-values from Wolak's test 1 (null is monotonic) and test 2 (alternative is monotonic)
%
%  Andrew Patton
%
%  5 May 2008

% working well. slow, but OK.

[T,K] = size(data);

if nargin<2 || isempty(reps);
    reps=100;
end

options = optimset('Display','off','TolCon',10^-8,'TolFun',10^-8,'TolX',10^-8);
Aineq = -eye(K);
Bineq = zeros(K,1);

muhat = mean(data)';  % unconstrained estimate of mu
omegahat = newey_west(data, floor(4*((T/100)^(2/9))) )/T; % HAC estimator of the covariance matrix of muhat
warning off;  % trying to speed up the code - eliminating warning messages printed to the screen
mutilda = fmincon('constrained_mean',muhat,Aineq,Bineq,[],[],[],[],[],options,muhat,omegahat);
IU = (muhat-mutilda)'*inv(omegahat)*(muhat-mutilda);  % the first test stat, equation 16
EI = mutilda'*inv(omegahat)*mutilda;  % the second test statistic, see just after equation 18 of Wolak (1989, JoE)
warning on;


% next: use monte carlo to obtain the weights for the weighted sum of chi-squareds
weights = zeros*ones(1+size(omegahat,1),1);
for jj=1:reps;
    tempdata = mvnrnd(zeros(1,size(omegahat,1)),omegahat,1);  % simulating iid Normal data
    warning off;
    mutilda1 = fmincon('constrained_mean',tempdata',Aineq,Bineq,[],[],[],[],[],options,tempdata',omegahat);
    warning on;
    temp = sum(mutilda1>0);  % counting how many elements of mutilda are greater than zero
    weights(1+temp) = weights(1+temp) + 1/reps;  % adding one more unit of weight to this element of the weight vector
end
out1(1) = 1-weighted_chi2cdf(IU,weights(end:-1:1));  % pvalue from wolak's first test
out1(2) = 1-weighted_chi2cdf(EI,weights);  % pvalue from wolak's second test
