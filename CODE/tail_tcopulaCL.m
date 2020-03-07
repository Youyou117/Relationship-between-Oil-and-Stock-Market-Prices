function CL = tail_tcopulaCL(theta,data,q)
% function CL = tail_tcopulaCL(theta,data,q)
%
% The negative copula log-likelihood of a Student's t copula for the lower (q,q) tail
%
% INPUTS:   theta = [rho, nu] = [ correlation coeff, degree of freedom parameter]
%			data = [U V];
%           q, a scalar in (0,1), the value defining the tail
%
% OUTPUTS:  CL, a scalar, the negative log-likelihood
%
% 1 Sep 2011, revised 13apr12
%
% Andrew Patton

rho = theta(1);
nu  = theta(2);

rng('default');  % re-setting RNG seed so that we get the same answer every time we evaluate the likelihood
tail_prob = 1-2*q+tCopula(q,q,rho,nu);

data = 1-data;  % doing this so that this code matches the Gumbel copula code (where we want to use the upper rather than the lower tail - for Stud t this makes no difference)
% so now we condition on being in the *upper* (q,q) tail

% using Chen, Fan, Pouzo and Ying (JoE 2010) on copulas for censored data
data1 = [max(data(:,1),1-q),max(data(:,2),1-q)];    % censoring the data at 1-q
delta = [(data(:,1)>=1-q),(data(:,2)>=1-q)];        % equals 0 if data is censored, -1 if not.


dataT = nan(size(data,1),2);
dataT(delta(:,1)==0,1) = tdis_inv(1-q,nu)*ones(sum(delta(:,1)==0),1);
dataT(delta(:,2)==0,2) = tdis_inv(1-q,nu)*ones(sum(delta(:,2)==0),1);
dataT(delta(:,1)==1,1) = tdis_inv(data(delta(:,1)==1,1),nu);
dataT(delta(:,2)==1,2) = tdis_inv(data(delta(:,2)==1,2),nu);


CL1 = -tcopulaCLa(theta,data);  % log-likelihood (not negative likelihood)
% think the two lines below are correct, but will try 1/c rather than c (CHECKED: the lines below ARE correct.
CL2 = log( tdis_cdf( (dataT(:,2)-rho*dataT(:,1)).*sqrt( 1/(1-rho^2)*(nu+1)./(nu+(dataT(:,1).^2)) ), nu+1) );   % log( dC(u,v)/du )
CL3 = log( tdis_cdf( (dataT(:,1)-rho*dataT(:,2)).*sqrt( 1/(1-rho^2)*(nu+1)./(nu+(dataT(:,2).^2)) ), nu+1) );   % log( dC(u,v)/dv )

%CL2 = log( tdis_cdf( (dataT(:,2)-rho*dataT(:,1)).*sqrt( 1/(1-rho^2)/(nu+1).*(nu+(dataT(:,1).^2)) ), nu+1) );   % log( dC(u,v)/du )
%CL3 = log( tdis_cdf( (dataT(:,1)-rho*dataT(:,2)).*sqrt( 1/(1-rho^2)/(nu+1).*(nu+(dataT(:,2).^2)) ), nu+1) );   % log( dC(u,v)/dv )

CL = delta(:,1).*delta(:,2).*CL1 + delta(:,1).*(1-delta(:,2)).*CL2 + (1-delta(:,1)).*delta(:,2).*CL3 + (1-delta(:,1)).*(1-delta(:,2))*log(tail_prob);
CL = sum(CL);
CL = -CL;








if 0
    
    % THIS IS NOT CORRECT - NEED DIFF FORM FOR LIKELIHOOD. SEE "tail_Gumbel_copulaCL.m" FOR CORRECT CODE WHEN ASSUMING GUMBEL COPULA. COME BACK AND FIX
    % THIS CODE LATER..
    
    rho = theta(1);
    nu  = theta(2);
    
    CLa = -tcopulaCLa(theta,data);  % log-like not negative log-like
    
    
    rng('default');  % re-setting RNG seed so that we get the same answer every time we evaluate the likelihood
    tail_prob = tCopula(q,q,rho,nu);
    
    CL = sum( (CLa - log(tail_prob)).*(data(:,1)<=q).*(data(:,2)<=q) );  % summing likelihood over points that lie in the (q,q) lower tail
    CL = -CL;
end
