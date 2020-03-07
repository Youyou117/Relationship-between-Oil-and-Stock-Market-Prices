function [out1] = dist_test_PEE(T,reps,ARparams,GARCHparams,SKEWTparams,GARCHstartvals);
%function [KSstats CvMstats] = dist_test_PEE(T,reps,ARparams,GARCHparams,SKEWTparams);
%
%  Function to simulate from an AR-GJRGARCH-SKEWt time series model, estimate the parameters of each, then compute the Kolmogorov-Smirnov and
%  Cramer-von Mises test statistics on the estimated probability integral transforms. This gives us a simulated distribution of these test stats
%  taking into account the estimation error from the unknown parameters.
%
%  INPUTS:  T, a scalar, the sample size
%           reps, a scalar, the number of replications to use
%           ARparams, a p+1 x 1 vector, [phi0,phi1,...phiP] the AR(p) parameters for the mean
%           GARCHparams, a 4x1 vector, [w,a,d,b] the parameters of the GJRGARCH model for the variance
%           SKEWTparams, a 2x1 vector, [nu,lam], the parameters of the skewed t distribution
%           GARCHstartvals, a 4x1 vector, starting values for GARCH estimation (default = [])
%
%  OUTPUTS: out1, a repsx2 matrix, [KSstats, CvMstats] the KS and CvM test stats on each of the simulated series
%
%  Andrew Patton
%
%  1 September 2011

ARp = length(ARparams)-1;
options = optimset('Display','off','TolCon',10^-12,'TolFun',10^-4,'TolX',10^-6,'Algorithm','interior-point');


out1 = nan(reps,2);
for rr=1:reps;
    eps = skewtdis_rnd(SKEWTparams(1),SKEWTparams(2),T);  % simulating the standardized residuals
    
    mu = nan(T,1);
    ht = nan(T,1);
    yt = nan(T,1);
    et = nan(T,1);
    
    % starting values for series
    if ARp>0
        mu(1) = ARparams(1)/(1-sum(ARparams(2:end)));
    else
        mu(1) = ARparams(1);
    end
    
    ht(1) = GARCHparams(1)/(1-GARCHparams(2)-0.5*GARCHparams(3)-GARCHparams(4));
    et(1) = sqrt(ht(1))*eps(1);
    yt(1) = mu(1) + et(1);

    % looping through rest of series
    for tt=2:T;
        ht(tt) = GARCHparams(1) + GARCHparams(2)*(et(tt-1)^2) + GARCHparams(3)*(et(tt-1)^2)*(et(tt-1)<0) + GARCHparams(4)*ht(tt-1);
        et(tt) = sqrt(ht(tt))*eps(tt);
        if tt<=ARp
            mu(tt) = mu(1);  % using unconditional mean for first ARp obs of conditional mean
        else
            mu(tt) = ARparams(1);
            if ARp>0
                mu(tt) = mu(tt) + ARparams(2:end)'*yt(tt-1:-1:tt-ARp);
            end
        end
        yt(tt) = mu(tt) + et(tt);
    end
    
    % estimating mean model
    temp = ols(yt(ARp+1:end),[ones(T-ARp,1),mlag(yt(ARp+1:end),ARp,mu(1))]);
    ehat = [zeros(ARp,1);temp.resid];
    
    % estimating the vol model
    [parameters4, likelihood1, ~, ~, hhat1] = multigarch(ehat,1,1,1,'GJRGARCH','NORMAL',[],GARCHstartvals);
    epshat = ehat./sqrt(hhat1);
    
    % estimating skew t model
    lower = [2.1, -0.99];
    upper = [100, 0.99 ];
    theta0 = [6;0];
    theta = fmincon('skewtdis_LL',theta0,[],[],[],[],lower,upper,[],options,epshat);
    Uhat = skewtdis_cdf(epshat,theta(1),theta(2));
    
    % the test stats
    KSstat = max(abs(sort(Uhat) - (1:T)'/T));
    CvMstat = sum( (sort(Uhat)-(1:T)'/T).^2 );
    out1(rr,:) = [KSstat,CvMstat];
end
    
    
