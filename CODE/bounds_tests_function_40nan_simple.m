function [out1, out2] = bounds_tests_function_40nan_simple(actuals,forecasts,print_table)
%
% Function to compute a battery of tests of forecast rationality, described in the paper:
%
% Patton, A.J., and A. Timmermann, 2012, Forecast Rationality Tests Based on Multi-Horizon Bounds, 
% Journal of Business and Economic Statistics, 30(1), 1-17.
% Paper available at: http://public.econ.duke.edu/~ap172/Patton_Timmermann_bounds_JBES_2012.pdf
%
% ** Allows for NaN to fill in missing values, so that all possible data can be used (not just balanced panels of data)
%
% INPUTS:   actuals, a Tx1 vector of realizations for the target variable (may be left empty, and then only a sub-set of the tests will be run)
%           forecasts, a TxH matrix of forecasts, with shortest forecast horizon in first column and longest forecast horizon in last column
%           print_table, =1 if want a nice print out of table to screen, =0 (default) if do not want the printout
%
% OUTPUTS:  out1, a 13x1 vector of p-values for the various tests:
%                       BOUNDS TESTS USING DATA ON ACTUAL
%                           1)  INC MSE test
%                           2)  DEC COV test (covariance between actual and forecast)
%                           3)  Covariance bound test 1: V[ nu[t|1,h] ] <= 2*Cov[ Y[t], nu[t|1,h] ] for all h>=2
%
%                       BOUNDS TESTS USING DATA ONLY ON FORECASTS
%                           4)  DEC MSF test (mean-squared forecast)
%                           5)  INC MSFR test (mean-squared forecast revision)
%                           6)  DEC COV test, using shortest-horizon forecast as proxy for actual
%                           7)  Covariance bound test 1: V[ nu[t|1,h] ] <= 2*Cov[ Yhat[t|t-1], nu[t|1,h] ] for all h>=3, using shortest-horizon forecast as proxy for actual
%
%                       REGRESSION BASED TESTS
%                           8)  Univar MZ on shortest horizon
%                           9)  Univar opt revision regression
%                           10) Univar opt revision regression, with h=1 forecast as proxy
%
%                       BONFERRONI COMBINATIONS OF ABOVE TESTS
%                           11) Bonf using actuals, tests 1, 2, 3, 8, 9
%                           12) Bonf using forecasts, tests 4, 5, 6, 7, 10
%                           13) Bonf across all tests, tests 1-10
%
%           out2, a (6+2+3)xH matrix of results from simple univariate MZ regressions. Rows contain: (alpha, t-alpha, beta, t-beta (diff from 1), chi2stat, chi2pval, Bias, t-stat, MSFE, MSF,V[F]). Cols cover horizons 1 to H
%
% (Note: p-values for Bonferroni tests will lie between 0 and 1, but will not be uniform under the null. The decision rule is to reject rationality if the "p-value" is less than the desired size, eg, 0.05.)
%
%
%  Andrew Patton
%
%  1 December 2010.


[T,H] = size(forecasts);
T1 = size(actuals,1);  % if "actuals" is empty, then T1=0

if nargin<3 || isempty(print_table)
    print_table=0;
end

options = optimset('Display','off','TolCon',10^-6,'TolFun',10^-6,'TolX',10^-6);

out1 = nan(13,1);
if T1==T;  % then we have data on the actuals and so can run tests 1, 2, 3, 8, 9
    if H>2  % need more than one horizon for tests 1, 2, 3
        % 1) INC MSE test
        e2 = (actuals(1:T)*ones(1,H)-forecasts).^2;
        out1(1) = wolak_test_kan_1( e2(:,2:end)-e2(:,1:end-1) );
        
        % 2) DEC COV test (covariance between actual and forecast) <=> DEC E[y[t]*yhattilda[t|t-h]
        yyhat = (actuals*ones(1,H)).*forecasts;
        out1(2) = wolak_test_kan_1(  yyhat(:,1:end-1)-yyhat(:,2:end)  );
        
        % 3) Covariance bound test 1: V[ nu[t|1,h] ] <= 2*Cov[ Y[t], nu[t|1,h] ] for all h>=2
        % <=>  E[ nu[t|1,h]^2 ] <= 2*E[ Y[t]*nu[t|1,h] ]
        % <=>  2*E[ Y[t]*nu[t|1,h] ] - E[ nu[t|1,h]^2 ] >= 0
        fr = forecasts(:,1)*ones(1,H-1) - forecasts(:,2:H);
        temp2 = 2*(actuals*ones(1,H-1)).*fr - (fr.^2);
        out1(3) = wolak_test_kan_1( temp2  );
        
        % 9) Univar opt revision regression
        temp123 = ~isnan(sum([actuals,forecasts]')');  % a vector of 1 and 0 with a 1 if all forecasts and the actual are available, 0 if one or more of these are missing
        temp = ols(actuals(temp123),[ones(sum(temp123),1),forecasts(temp123,H),forecasts(temp123,1:H-1)-forecasts(temp123,2:H)]);
        chi2stat = (temp.beta-[0;ones(H,1)])'*inv(temp.vcv)*(temp.beta-[0;ones(H,1)]);
        out1(9) = 1-chi2cdf( chi2stat, 1+H );
    end
    
    % 8) Univar MZ regressions on shortest horizon
    temp124 = sum(isnan([actuals,forecasts(:,1)])')==0;  % all the rows where both actual and forecast are non-NAN
    temp = nwest(actuals(temp124),[ones(sum(temp124),1),forecasts(temp124,1)],3);  % newey-west truncation lag of 3
    chi2stat = (temp.beta-[0;1])'*inv(temp.vcv)*(temp.beta-[0;1]);
    out1(8) = 1-chi2cdf( chi2stat, 2 );
end
if H>1;  % then we have many horizons to play with, so can run some more tests. These tests don't require the actual so can run them even it T1=0
    % 10) Univar opt revision regression, with h=1 forecast as proxy
    temp123 = ~isnan(sum(forecasts')');  % a vector of 1 and 0 with a 1 if all forecasts are available, 0 if one or more of these are missing
    temp = ols(forecasts(temp123,1),[ones(sum(temp123),1),forecasts(temp123,H),forecasts(temp123,2:H-1)-forecasts(temp123,3:H)]);
    chi2stat = (temp.beta-[0;ones(H-1,1)])'*inv(temp.vcv)*(temp.beta-[0;ones(H-1,1)]);
    out1(10) = 1-chi2cdf( chi2stat, 1+H-1 );
    
    % 4) DEC MSF test (mean-squared forecast)
    out1(4) = wolak_test_kan_1( forecasts(:,1:H-1).^2 - forecasts(:,2:H).^2  );
    
    if H>2
        % 5) INC MSFR test (mean-squared forecast revision)
        fr = forecasts(:,1)*ones(1,H-1) - forecasts(:,2:H);
        out1(5) = wolak_test_kan_1( (fr(:,2:end).^2) - (fr(:,1:end-1).^2)  );
        
        % 6) DEC COV test, using shortest-horizon forecast as proxy for actual
        yyhat = (forecasts(:,1)*ones(1,H-1)).*forecasts(:,2:H);
        out1(6) = wolak_test_kan_1(  yyhat(:,1:end-1)-yyhat(:,2:end) );
        
        % 7) Covariance bound test 1: V[ nu[t|1,h] ] <= 2*Cov[ Yhat[t|t-1], nu[t|1,h] ] for all h>=3, using shortest-horizon forecast as proxy for actual
        fr2 = forecasts(:,2)*ones(1,H-2) - forecasts(:,3:H);
        temp2 = 2*(forecasts(1:T,1)*ones(1,H-2)).*fr2 - (fr2.^2);
        out1(7) = wolak_test_kan_1( temp2 );
    end
end

% Bonferroni tests
if T1==T;
    out1(11) = nanmin(out1([1, 2, 3, 8, 9]))*sum(~isnan(out1([1, 2, 3, 8, 9])));  % min p-value scaled by the number of tests that were ran
    out1(11) = min(out1(11),1);  % just making sure that this "pvalue" is less than one. (rejection requires pval to be less than 0.05 or 0.10, so this transformation does not change any conclusions, just makes the p-value look more "normal")
end
if H>1 
    out1(12) = nanmin(out1([4, 5, 6, 7, 10]))*sum(~isnan(out1([4, 5, 6, 7, 10])));  % min p-value scaled by the number of tests that were ran
    out1(12) = min(out1(12),1);  
end
out1(13) = nanmin(out1(1:10))*sum(~isnan(out1(1:10)));
out1(13) = min(out1(13),1);

out2 = nan(15,H);
if T1==T && (nargout>1  || print_table==1)  % then do the univariate MZ tests
    for hh=1:H;
        temp123 = ~isnan(sum([actuals,forecasts(:,hh)]')');  % a vector of 1 and 0 with a 1 if all forecasts and the actual are available, 0 if one or more of these are missing
        temp = nwest(actuals(temp123),[ones(sum(temp123),1),forecasts(temp123,hh)],3);
        out2([end-5;end-3],hh) = temp.beta;
        out2([end-4;end-2],hh) = (temp.beta-[0;1])./temp.se;
        out2(end-1,hh) = (temp.beta-[0;1])'*inv(temp.vcv)*(temp.beta-[0;1]);
        out2(end,hh) = 1-chi2cdf(out2(end-1,hh),2);
        
        temp = nwest(actuals(temp123)-forecasts(temp123,hh),[ones(sum(temp123),1)]);
        out2(1,hh) = temp.beta;
        out2(2,hh) = temp.tstat;
        
        temp = nwest((actuals(temp123)-forecasts(temp123,hh)).*forecasts(temp123,hh),[ones(sum(temp123),1)]);
        out2(3,hh) = temp.beta;
        out2(4,hh) = temp.tstat;
        
        
        out2(5,hh) = mean( (actuals(temp123)-forecasts(temp123,hh)).^2 );
        out2(6,hh) = mean( (forecasts(temp123,hh)).^2 );
        out2(7,hh) = cov( forecasts(temp123,hh) );
        out2(8,hh) = out2(5,hh)+out2(7,hh);
        out2(9,hh) = cov( actuals(temp123) );
        
    end
end

if print_table==1;
    clear info;
    info.fmt = '%10.3f';
    info.cnames = strvcat('p-value');
    info.rnames = strvcat('# Test',...
        '1  INC MSE', '2  DEC COV', '3  COV bound', ...
        '4  DEC MSF', '5  INC MSFR', '6  DEC COV with proxy', '7  COV bound with proxy', ...
        '8  Univar MZ short h','9  Univar opt revision','10 Univar opt revision with proxy',...
        '11 Bonferroni 1 (1,2,3,8,9)','12 Bonferroni 2 (4,5,6,7,10)','13 Bonferroni All (1-10)');
    sprintf(['\nP-values from tests of forecast rationality, T=',int2str(T),', H=',int2str(H),'\n'])
    mprint(out1,info)
    
    clear info;
    info.fmt = '%10.3f';
    info.cnames = 'h=1';
    for hh=2:H;
        info.cnames(hh,:) = ['h=',int2str(hh)];
    end
    info.rnames = strvcat('Forecast horizon:','Bias','t-stat','E[error*forecast]','t-stat','Mean Sq Error','Mean Sq Forecast','V[forecast]','MSE+V[forecast]','V[actual]',...
        'MZ alpha','MZ t-stat','MZ beta','MZ t-stat','MZ chi2stat','MZ chi2pval');
    sprintf(['\nResults from univariate MZ tests for each horizon, T=',int2str(T),', H=',int2str(H),'\n(Note: t-stat on beta compares estimate to 1 not 0)\n'])
    mprint(out2,info)
    
end
