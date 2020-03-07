function [KSstat CVMstat] = copula_GOF_stats(data,copula_cdf_str,varargin)
% function [KSstat CVMstat] = copula_GOF_stats(data,copula_cdf_str,varagin)
%
%  Function to compute the Kolmogorov-Smirnov test statistic and the Cramer-von Mises test statistic on a given sample of data, for a given copula
%  model. (Note that no critical values are provided - these need to be obtained by simulation or bootstrap. This function can be used inside a
%  simulation or bootstrap to get critical values.)
%
%  INPUTS:  data, a Txk matrix of Unif(0,1) data 
%           copula_cdf_str, a string containing the name of the function that returns a Tx1 vector of the value of the copula CDF 
%           varargin, any other arguments to pass to the copula_cdf_str. (I assume that data is the first input to the CDF function)
%
%  OUTPUTS: KSstat, a scalar, the value of the KS test statistic
%           CVMstat, a scalar, the value of the CvM test statistic
%               
%  Andrew Patton
%
%  30 September 2011

[T,k] = size(data);

fitted_cdf = feval(copula_cdf_str,data(:,1),data(:,2),varargin{:});
%fitted_cdf = feval(copula_cdf_str,data,varargin{:});
empirical_cdf = empirical_copula(data);

KSstat  = max( abs(fitted_cdf-empirical_cdf) );
CVMstat = sum( (fitted_cdf-empirical_cdf).^2 );