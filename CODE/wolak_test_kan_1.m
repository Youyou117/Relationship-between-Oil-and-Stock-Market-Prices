function pvalue = wolak_test_kan_1(data)
% function pvalue = wolak_test_kan_1(data,reps)
%
% Function to implement Wolak's (1989, JoE) test of
%
%       H0*: d1>=0,d2>=0,...,dK>=0
%
%  vs.  H1*: (d1,d2,...,dK) in R^K  (ie: general alternative)
%
%  ** USING SUGGESTION FROM RAYMOND KAN FOR SPEEDING THIS UP 
%
% Used in the papers:
%
% Patton, A.J., and A. Timmermann, 2010, Monotonicity in Asset Returns: New Tests with Applications to 
% the Term Structure, the CAPM and Portfolios Sorts, Journal of Financial Economics, 98(3), 605-625.
% Paper available at: http://public.econ.duke.edu/~ap172/Patton_Timmermann_sorts_JFE_Dec2010.pdf
%
% Patton, A.J., and A. Timmermann, 2012, Forecast Rationality Tests Based on Multi-Horizon Bounds, 
% Journal of Business and Economic Statistics, 30(1), 1-17.
% Paper available at: http://public.econ.duke.edu/~ap172/Patton_Timmermann_bounds_JBES_2012.pdf
%
%
% INPUTS:   data, a TxK matrix of data (return *differences* for the sorting application)
%
% OUTPUTS:  pvalue, a scalar, the p-value from Wolak's test 1 (null is weakly positive) 
%
%  Andrew Patton
%
%  5 October 2008

% Based on wolak_test_5.m

[T,K] = size(data);


muhat = nanmean(data)';  % unconstrained estimate of mu
omegahat1 = nannewey_west(data, floor(4*((T/100)^(2/9))) )/T; % HAC estimator of the covariance matrix of muhat
flag=0;
if isposdef(omegahat1)   % if HAC estimate is pos def, then use that in the tests
    omegahat = omegahat1;
    flag=1;
else  % then try HAC estimator of the covariance matrix of muhat, with lag length set to zero
    omegahat2 = nannewey_west(data, 0 )/T;
    if isposdef(omegahat2)   % if HAC estimate is not pos def, but one with zero lags is, then use that
        omegahat = omegahat2;
        flag=1;
    else
        omegahat3 = nannewey_west(data, 0 ,0)/T;  % trying the non-pairwise versino of this function
        if isposdef(omegahat3)   % if HAC estimate is not pos def, but one with zero lags is, then use that
            omegahat = omegahat3;
            flag=1;
        end
    end
end  % else, we will have to skip this test, so leave the "flag" variable at zero

if flag==1
    [~,pvalue] = ineqtest(muhat,omegahat);
end