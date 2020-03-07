% INEQTEST.M
% This Matlab performs a test of
% H0: d>=0_r vs. H1: d in R^r.
% due to Wolak (1989, Journal of Econometrics)
%
% Input:
% dhat: estimate of d (a rx1 column vector)
% V: covariance matrix of dhat
% Output:
% LR: likelihood ratio test statistic
% pval: p-value of the test
%
%
% From Raymond Kan, 16 Sep 2010

function [LR,pval] = ineqtest(dhat,V)
w = waldweight(V);

if all(dhat>=0)
    LR = 0;
    pval = 1-w(end);
else
    %
    %   Use quadratic programming to construct the Kuhn-Tucker test statistic
    %   min_{lam} lam'*dhat+0.5*lam'*V*lam  s.t.  lam>=0_r
    %
    r = length(dhat);
    [lam,obj] = quadprog(V,dhat,-eye(r),zeros(r,1),[],[],[],[],[], optimset('Display','off','Diagnostics','off','LargeScale','off','Algorithm','active-set'));
    pval = gammainc(-obj,[1:r]./2,'upper')*w(r:-1:1);
    LR = -2*obj;
end
