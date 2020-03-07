function [out1, outboot] = rankcorrel(data,alpha,bootreps)
% function out1 = rankcorrel(data)
%
% Computes the sample estimator of Spearman's
% rank correlation for the given data
%
% USES THE STATIONARY BOOTSTRAP TO OBTAIN CRITICAL VALUES
%
% INPUTS:	data, a Tx2 vector of data
%           alpha = 0, if don't want confidence interval for correlation coefficient,
%                   else some number in (0,1) (eg 0.05) if want bootstrap
%                   confidence interval. block length for stationary bootstrap is
%                   obtained using algorithm in politis and white (2003)
%           bootreps, a scalar, the number of bootstrap replications to use
%
% OUTPUT:	out1, a scalar if alpha=0, a 4x1 vector else. Sample estimate of Spearman's rho
%                    with lower and upper confidence interval bounds, and bootstrap standard error
%
% Tuesday, 4 November, 2003
%
% Andrew Patton
%
% REVISED 24 March 2011 to include boostrap std errors, based on theory in Remillard (2010, wp)


if nargin<2 || isempty(alpha)
    alpha=0;
    bootreps=0;
end
if nargin<3 && alpha>0
    bootreps=1000;
end
if bootreps>0 && (alpha<=0 || isempty(alpha))
    alpha = 0.05;
end
T = size(data,1);

out1 = rankcorrel_calc(data(:,1),data(:,2));    % this is the rank correl coefficient


if bootreps>0
    bootdates = randi(T,T,bootreps);  % a set of dates to use in bootstrap (data are assumed iid so just shuffle dates without using blocks
    
    outboot = nan(bootreps,1);
    for bb=1:bootreps
        outboot(bb) = rankcorrel_calc(data(bootdates(:,bb),1),data(bootdates(:,bb),2));
    end
    out1(2:3) = quantile(outboot,[alpha/2,1-alpha/2]);
    out1(4) = std(outboot);
end




    function rhohat = rankcorrel_calc(X1,Y1);
        % helper function
   
        X1 = ranks(X1);
        Y1 = ranks(Y1);
        rhohat = corrcoef12(X1,Y1);    % this is the rank correl coefficient
    end


end