function CL = tail_Gumbel_copulaCL(theta,data,q)
% function CL = tail_tcopulaCL(theta,data,q)
% 
% The negative copula log-likelihood of a Gumbel copula for the lower (q,q) tail
%
% INPUTS:   theta = [rho, nu] = [ correlation coeff, degree of freedom parameter]
%			data = [U V];
%           q, a scalar in (0,1), the value defining the tail
%
% OUTPUTS:  CL, a scalar, the negative log-likelihood
%
% 1 Sep 2011
%
% Andrew Patton

% seems to be working nicely. It overestimates tail dependence when the true tail dependence is zero, as one might expect. But when tail dependence is
% present it gets it nicely.

data = 1-data;  % want to use the *upper* tail of the Gumbel copula, as this is the tail with non-zero tail deps
% so now we condition on being in the *upper* (q,q) tail

tail_prob = 1-2*q+gumbel_cdf(q,q,theta);  

% using Chen, Fan, Pouzo and Ying (JoE 2010) on copulas for censored data
data1 = [max(data(:,1),1-q),max(data(:,2),1-q)];    % censoring the data at 1-q
delta = [(data(:,1)>=1-q),(data(:,2)>=1-q)];        % indicator for whether censoring is present

CL1 = -gumbelCLa(theta,data1);  % log-likelihood (not negative likelihood)
CL2 = log(GumbelUgivenV_t(data1(:,2),data1(:,1),0,theta));  % log( dC(u,v)/du )
CL3 = log(GumbelUgivenV_t(data1(:,1),data1(:,2),0,theta));  % log( dC(u,v)/dv )

CL = delta(:,1).*delta(:,2).*CL1 + delta(:,1).*(1-delta(:,2)).*CL2 + (1-delta(:,1)).*delta(:,2).*CL3 + (1-delta(:,1)).*(1-delta(:,2))*log(tail_prob);
CL = sum(CL);
CL = -CL;

