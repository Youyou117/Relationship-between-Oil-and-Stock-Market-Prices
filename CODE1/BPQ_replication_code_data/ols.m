function results=ols(y,x)
% PURPOSE: least-squares regression 
%---------------------------------------------------
% USAGE: results = ols(y,x)
% where: y = dependent variable vector (nobs x 1)
%        x = independent variables matrix (nobs x nvar)
%---------------------------------------------------
% RETURNS: a structure
%        results.meth  = 'ols'
%        results.beta  = bhat
%        results.tstat = t-stats
%        results.se = standard errors on coefficient estimates
%        results.yhat  = yhat
%        results.resid = residuals
%        results.sige  = e'*e/(n-k)
%        results.rsqr  = rsquared
%        results.rbar  = rbar-squared
%        results.dw    = Durbin-Watson Statistic
%        results.nobs  = nobs
%        results.nvar  = nvars
%        results.y     = y data vector
%   	 results.vcv	 = the variance-covariance matrix of bhat (AJP 4feb04)
%   	 results.Fstatpval	 = p-value on joint test that all parameters except the first (ie, the constant) equal zero (AJP 22nov04)
%---------------------------------------------------
% SEE ALSO: prt(results), plt(results)
%---------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

if (nargin ~= 2); error('Wrong # of arguments to ols');
else
    [nobs nvar] = size(x);
    [nobs2 junk] = size(y);
    if (nobs ~= nobs2); error('x and y must have same # obs in ols');
    end;
end;

results.meth = 'ols';
results.y = y;
results.nobs = nobs;
results.nvar = nvar;

[q r] = qr(x,0);
xpxi = (r'*r)\eye(nvar);


results.beta = xpxi*(x'*y);
results.yhat = x*results.beta;
results.resid = y - results.yhat;
sigu = results.resid'*results.resid;
results.sige = sigu/(nobs-nvar);
results.vcv = results.sige*xpxi;		% the vcv matrix of bhat - AJP inserted 4feb04
tmp = (results.sige)*(diag(xpxi));
results.se = sqrt(tmp);
results.tstat = results.beta./(sqrt(tmp));
ym = y - mean(y);
rsqr1 = sigu;
rsqr2 = ym'*ym;
results.rsqr = 1.0 - rsqr1/rsqr2; % r-squared
rsqr1 = rsqr1/(nobs-nvar);
rsqr2 = rsqr2/(nobs-1.0);
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
ediff = results.resid(2:nobs) - results.resid(1:nobs-1);
results.dw = (ediff'*ediff)/sigu; % durbin-watson
if length(results.beta)>1
    results.Fstatpval = 1-chi2cdf(results.beta(2:end)'*inv(results.vcv(2:end,2:end))*results.beta(2:end),length(results.beta-1));
else
    results.Fstatpval=[];
end