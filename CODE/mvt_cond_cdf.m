function out1 = mvt_cond_cdf(X2,X1,mu,S,nu);
% function out1 = mvt_cond_cdf(X2,X1,mu,S,nu);
%
% Conditional CDF of X2|X1, where [X1,X2] ~ Student's t (mu,S,nu). 
%
% Based on pp284-286 of Box, Jenkins and Reinsel (1994, "Time Series Analysis")
%
% INPUTS:   X2, a Txp2 matrix
%           X1, a Txp1 vector
%           mu, a px1 vector (p1+p2=p) of the mean of X=[X1,X2]
%           S, a pxp matrix, the scale matrix of X
%           nu, a scalar, the degress of freedom parameter
%
%  Andrew Patton
%
%  17 October 2011

[T,p2] = size(X2);
p1 = size(X1,2);
p = p1+p2;

S12 = S(1:p1,p1+1:end);  % Matrix S is broken into [[S11,S12];[S12',S22]]
S11 = S(1:p1,1:p1);
S22 = S(1+p1:end,1+p1:end);

beta21 = S12'*inv(S11);
S2211 = S22 - S12'*inv(S11)*S12;  

mu1 = mu(1:p1);
mu2 = mu(1+p1:end);

out1 = nan(T,1);
for tt=1:T;
    cc = (nu+p1)/(nu+ (X1(tt,:)'-mu1)'*inv(S11)*(X1(tt,:)'-mu1));  % constant that appears in the scaling term. Note that it changes with the value of the conditioning variable, unlike the Normal case
    mu21 = mu2 + beta21*(X1(tt,:)'-mu1);  % conditional mean of X2 given X1
    
    V21 = 1/cc*S2211;  % conditional scale matrix of X2 given X1
    
    D = diag(sqrt(diag(V21)));  % diagonal matrix with sqrt of variances on diagonal. Needed because built-in mvtcdf requires the variable to have variances equal to one
    R21 = inv(D)*V21*inv(D);    % conditional correlation matrix of X2 given X1

    out1(tt) = mvtcdf( inv(D)*(X2(tt,:)' - mu21) , R21 , nu+p1);
end
