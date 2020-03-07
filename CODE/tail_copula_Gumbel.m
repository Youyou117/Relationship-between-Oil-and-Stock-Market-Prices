function [tauL,tauU,thetahatL,thetahatU,nobs,outboot] = tail_copula_Gumbel(data,q,bootreps);
% function [tauL,tauU,thetahatL,thetahatU] = tail_copula_Gumbel(data,q);
%
%  Function to estimate the upper and lower tail dependence coefficients by approximating the upper and lower tail copulas with the Gumbel copula
%  (estimated separately for each tail).
%
%  INPUTS:  data, a Tx2 matrix of Unif(0,1) data (covering the full support)
%           q, a scalar in (0,1), the quantile to use to define the lower tail (upper tail will be 1-q)
%           bootreps, a scalara, the number of bootstrap replications to use to get confidence intervals (default=0)
%
%  OUTPUTS: tauL, a scalar, the estimated lower tail dependence coefficient
%           tauU, a scalar, the estimated upper tail dependence coefficient
%           thetahatL, a scalar, the estimated parameter of the Gumbel copula for the lower tail
%           thetahatU, a scalar, the estimated parameter of the Gumbel copula for the upper tail
%           nobs, a 1x2 vector, the number of observations used in estimation
%
%  Andrew Patton
%
%  1 Sep 2011


T = size(data,1);

if nargin<2 || isempty(q)
    q=0.1;
end
if nargin<3 || isempty(bootreps) || bootreps==0
    bootreps=0;
    outboot = nan;
end


options = optimset('Display','off','TolCon',10^-12,'TolFun',10^-4,'TolX',10^-6,'Algorithm','interior-point');
lower = 1.0001;
upper = 10;
theta0 = 2;

nobs(1) = sum( (data(:,1)<=q).*(data(:,2)<=q) );
nobs(2) = sum( (data(:,1)>=(1-q)).*(data(:,2)>=(1-q)) );

% lower tail
if nobs(1)>3
    thetahatL  = fmincon('tail_Gumbel_copulaCL',theta0,[],[],[],[],lower,upper,[],options,data,q);  
else
    thetahatL=nan;
end
tauL = 2- (2^(1/thetahatL));

% upper tail
if nobs(2)>3
    thetahatU = fmincon('tail_Gumbel_copulaCL',theta0,[],[],[],[],lower,upper,[],options,1-data,q);
else
    thetahatU =nan;
end
tauU = 2- (2^(1/thetahatU));

if bootreps>0
    bootdates = randint(1,T,T,bootreps);  % a set of dates to use in bootstrap (data are assumed iid so just shuffle dates without using blocks
    outboot = nan(bootreps,4);
    for bb=1:bootreps
        [tauLb,tauUb,thetahatLb,thetahatUb] = tail_copula_Gumbel(data(bootdates(:,bb),:),q,0);
        outboot(bb,:) = [tauLb,tauUb,thetahatLb,thetahatUb];
    end
end    
