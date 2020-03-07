function [pvals,teststats,bootstats] = break_test_tv_copula_2(data,break_point,bootreps)
% function [pvals,teststats,bootstats] = break_test_tv_copula_1(data,break_point,bootreps)
%
%  Function to test for a time-varying copula through a break in rank correlation at a given point in the sample
%
% INPUTS:   data, a Tx2 matrix of (approx) Unif(0,1) data
%           break_point,   a scalar in (0,1), the fraction of the sample at which the break occurs
%                       OR a scalar in (-1,0), to *search* for a break anywhere in the interval [ -breakpoint*T, (1+breakpoint)*T ]
%           bootreps, a scalar, the number of bootstrap replications (default=1000)
%           QQ, a kx1 vector, some values of q to use in quantile dependence
%
% OUTPUTS:  pval, a scalar, the p-values for the tests of time-varying dependence,
%           teststats, a scalar, the test stats for the tests of time-varying dependence
%           bootstats, a bootreps x 1 matrix, the test stats on the bootstrap data
%
%  Andrew Patton
%
%  24 May 2011

T = size(data,1);

if nargin<2 || isempty(break_point)
    break_point = 0.5;
end
if nargin<3 || isempty(bootreps)
    bootreps=1000;
end

bootstats = nan(bootreps,1);
bootdates = randi(T,T,bootreps);  % these are the "dates" to use in the iid boostrap
if break_point>0  % then date of the break is known
    % rank correlation
    X = [[ones(floor(T*break_point),1);zeros(T-floor(T*break_point),1)],ones(T,1)];  % so just test for signif of beta 1
    UV = data(:,1).*data(:,2);
    Y = UV;
    xpxi = (X'*X)\eye(2);
    beta = xpxi*(X'*Y);
    e = Y-X*beta;
    vcv = (e'*e)/(T-2)*xpxi;
    teststat = beta(1)'*(vcv(1,1)\beta(1));
    
    for bb = 1:bootreps;
        % rank correlation
        UV = data(bootdates(:,bb),1).*data(bootdates(:,bb),2);
        Y = UV;
        xpxi = (X'*X)\eye(2);
        beta = xpxi*(X'*Y);
        e = Y-X*beta;
        vcv = (e'*e)/(T-2)*xpxi;
        bootstats(bb,1) = beta(1)'*(vcv(1,1)\beta(1));
    end
else  % we need to consider searching over the break
    
    if break_point>-0.5  % then use this as the start of the sample
        tmin = ceil(-break_point*T);
        tmax = floor((1--break_point)*T);
    else
        tmax = floor(-break_point*T);
        tmin = ceil((1--break_point)*T);
    end
    
    teststat = 0;
    UV = data(:,1).*data(:,2);
    for tt=tmin:tmax
        X = [[ones(tt,1);zeros(T-tt,1)],ones(T,1)];  % so just test for signif of beta 1
        Y = UV;
        xpxi = (X'*X)\eye(2);
        beta = xpxi*(X'*Y);
        e = Y-X*beta;
        vcv = (e'*e)/(T-2)*xpxi;
        teststat1 = beta(1)'*(vcv(1,1)\beta(1));
        teststat = max(teststat,teststat1);
    end
    
    for bb = 1:bootreps;
        UV = data(bootdates(:,bb),1).*data(bootdates(:,bb),2);
        Y = UV;
        bootstat = 0;
        for tt=tmin:tmax
            X = [[ones(tt,1);zeros(T-tt,1)],ones(T,1)];  % so just test for signif of beta 1
            xpxi = (X'*X)\eye(2);
            beta = xpxi*(X'*Y);
            e = Y-X*beta;
            vcv = (e'*e)/(T-2)*xpxi;
            bootstat1 = beta(1)'*(vcv(1,1)\beta(1));
            bootstat = max(bootstat,bootstat1);
        end
        bootstats(bb,1) = bootstat;
    end
end    
pvals = mean(bootstats(:,1)>teststat);  % proportion of simulations that generated a test stat bigger than the one observed in the data
    