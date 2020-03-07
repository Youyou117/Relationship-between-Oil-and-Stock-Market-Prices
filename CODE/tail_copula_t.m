function [tauL,tauU,thetahatL,thetahatU,nobs,outboot] = tail_copula_t(data,q,bootreps);
% function [tauL,tauU,thetahatL,thetahatU] = tail_copula_t(data,q);
%
%  Function to estimate the upper and lower tail dependence coefficients by approximating the upper and lower tail copulas with the Student's t copula
%  (estimated separately for each tail).
%
%  INPUTS:  data, a Tx2 matrix of Unif(0,1) data (covering the full support)
%           q, a scalar in (0,1), the quantile to use to define the lower tail (upper tail will be 1-q)
%
%  OUTPUTS: tauL, a scalar, the estimated lower tail dependence coefficient
%           tauU, a scalar, the estimated upper tail dependence coefficient
%           thetahatL, a 2x1 vector, the estimated parameters of the Student's t copula for the lower tail
%           thetahatU, a 2x1 vector, the estimated parameters of the Student's t copula for the upper tail
%
%  Andrew Patton
%
%  1 Sep 2011

T = size(data,1);

if nargin<2 || isempty(q)
    q=0.1;
end
if nargin<3 || isempty(bootreps)
    bootreps=0;
end


[tauL,tauU,thetahatL,thetahatU] = tail_copula_t_calc(data,q,[0.3;6],[0.3;6]);


nobs(1) = sum( (data(:,1)<=q).*(data(:,2)<=q) );
nobs(2) = sum( (data(:,1)>=(1-q)).*(data(:,2)>=(1-q)) );


if bootreps>0
    bootdates = randi(T,T,bootreps);  % a set of dates to use in bootstrap (data are assumed iid so just shuffle dates without using blocks
    outboot = nan(bootreps,6);
    for bb=1:bootreps
        thetahatLo = thetahatL + rand(2,1)*0.01-0.005;  % want to start estimation off near actual estimate, but not exactly as sometimes the search algorithm does not leave the starting value.
        thetahatUo = thetahatU + rand(2,1)*0.01-0.005;
        [tauLb,tauUb,thetahatLb,thetahatUb] = tail_copula_t_calc(data(bootdates(:,bb),:),q,thetahatLo,thetahatUo);
        outboot(bb,:) = [tauLb,tauUb,thetahatLb',thetahatUb'];
        bb
    end
end


    function [tauLb,tauUb,thetahatLb,thetahatUb] = tail_copula_t_calc(data,q,theta0L,theta0U);
        
        options = optimset('Display','off','TolCon',10^-12,'TolFun',10^-4,'TolX',10^-6,'Algorithm','interior-point');
        lower = [-0.99;2.1];
        upper = [0.99,100];
        
        % lower tail
        thetahatLb  = fmincon('tail_tcopulaCL',theta0L,[],[],[],[],lower,upper,[],options,data,q);
        tauLb = 2*tdis_cdf(-sqrt( (thetahatLb(2)+1)*(1-thetahatLb(1))/(1+thetahatLb(1))), (thetahatLb(2)+1));
        
        % upper tail
        thetahatUb = fmincon('tail_tcopulaCL',theta0U,[],[],[],[],lower,upper,[],options,1-data,q);
        tauUb = 2*tdis_cdf(-sqrt( (thetahatUb(2)+1)*(1-thetahatUb(1))/(1+thetahatUb(1))), (thetahatUb(2)+1));
    end

end
