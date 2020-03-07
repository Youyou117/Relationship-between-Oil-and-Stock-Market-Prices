function [teststats,pvals] = chen_fan_PLR_test(theta1,theta2,data,LL1_str,LL2_str,v1,v2,varargin)
% function out1 = chen_fan_PLR_test(theta1,theta2,data,LL1_str,LL2_str);
%
%  Function to implement tests from Chen and Fan (2006, JoE) 
%       1) Test for generalized nestedness (p139 of Chen and Fan)
%       2) PLR test for generalized NON-NESTED models
%       3) PLR test for generalized NESTED models
%
% INPUTS:   theta1, a k1x1 vector, estimated parameter for first copula model
%           theta2, a k2x1 vector, estimated parameter for second copula model
%           data, a Txd matrix of (approx) Unif(0,1) data
%           LL1_str, a string containing name of log-likelihood function for first model that returns the value of the log-likelihood *at each point in time*
%           LL2_str, a string containing name of log-likelihood function for second model that returns the value of the log-likelihood *at each point in time*
%       	v1, a scalar indicating the number of additional arguments to pass to LL1_str
%			v2, a scalar indicating the number of additional arguments to pass to LL2_str
%			varargin, arguments to pass to LL1_str and LL2_str
%
% OUTPUTS:  teststats, a 3x1 vector, [stat1, stat2, stat3], teststats for three tests listed above. (under their respective nulls, stat1 is a weighted chi2, stat2 is N(0,1), and stat3 is weighted chi2)
%           pvals, a 3x1 vector, [pval1,pval2,pval3] pvalues for the three tests listed above
%
%  Andrew Patton
%
%  4 April 2011

[T,d] = size(data);
k1 = length(theta1);
k2 = length(theta2);

if nargin<6 || isempty(v1)
    v1=0;
end
if nargin<7 || isempty(v2)
    v2=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computing matrices that relate to asymptotic covariance of estimated parameters, and that are used in the tests below
scores1 = LLgrad_1(LL1_str,theta1,data,varargin{1:v1});  
scores1 = -scores1;     % taking negative as my Matlab likelihoods are all *negative* likelihoods (so that I can use fmincon)
B1 = cov(scores1);
[sigma1,~,~,scoresS1,Qjt1] = copula_deriv2_chen_fan(LL1_str,theta1,data,varargin{1:v1});   % this is a slow line, around 13 seconds

scores2 = LLgrad_1(LL2_str,theta2,data,varargin{1+v1:v1+v2});
scores2 = -scores2;     % taking negative as my Matlab likelihoods are all *negative* likelihoods (so that I can use fmincon)
B2 = cov(scores2);
[sigma2,~,~,scoresS2,Qjt2] = copula_deriv2_chen_fan(LL2_str,theta2,data,varargin{1+v1:v1+v2});

sigma12 = nan(k1,k2);  % defined on p136, just above eq 4.1
for ii=1:k1;
    for jj=1:k2;
        sigma12(ii,jj) = cov12(scoresS1(:,ii),scoresS2(:,jj));
    end
end

What = [ [sigma2*inv(B2) , -sigma12'*inv(B1)] ; [sigma12*inv(B2) , -sigma1*inv(B1)] ];  % first line of p137
eigW = eig(What);

LL1a = -feval(LL1_str,theta1,data,varargin{1:v1});  % note taking negative here as my likelihoods are all *negative* log-likelihoods
LL2a = -feval(LL2_str,theta2,data,varargin{1+v1:v1+v2});
LR = LL2a-LL1a;   
LRn = mean(LR);     % middle of p135
sigma2a = mean( ( LR - LRn ).^2 );   % eq 4.7 on p138

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the test for generalized nestedness. Reject null of "generalized nestedness" if pval is less than, say, 0.05. (See Theorem 4.5, p139)
stat1 = T*sigma2a;
pval1 = 1-weighted_chi2_1df_cdf(stat1,eigW.^2);  % top of p139. Note the use of squared eigenvalues here


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pseudo-LR test for NON-NESTED case
Qjt1s = sum(Qjt1,2);        % these are the two summation (j=1,..,d) terms that appear in eq 4.4
Qjt2s = sum(Qjt2,2); 
sigma2hat = mean( (LR - LRn + Qjt1s-Qjt2s ).^2 );  % eq 4.4 on p137
stat2 = sqrt(T)*LRn/sqrt(sigma2hat);              % eq 4.5 on p138
pval2 = 2*normcdf(-abs(stat2));                         % two-sided p-value for this statistic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pseudo-LR test for NESTED case
stat3 = 2*T*LRn;  % eq 4.6 on p138
pval3 = 1-weighted_chi2_1df_cdf(stat1,eigW);  % part (2) of Theorem 4.2 on p136

teststats = [stat1;stat2;stat3];
pvals = [pval1;pval2;pval3];



