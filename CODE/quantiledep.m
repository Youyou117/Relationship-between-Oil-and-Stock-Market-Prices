function out1 = quantiledep(X,Y,q)
% Computes the sample estimator of the quantile
% dependence measure ('near tail' dependence)
% 
% INPUTS:	X, a Tx1 vector of data
%				Y, a Tx1 vector of data
%				q, a kx1 vector of quantiles to evaluate the dependence measure at
%
% OUTPUT:	out1, a qx1 vector of the measure at each quantile
%
% Tuesday, 24 April, 2001
%
% Andrew Patton
%
% see also: quantiledepCOP.m


T = size(X,1);
k = size(q,1);

out1 = -999.99*ones(k,1);
for jj = 1:k;
   if q(jj)<=0.5;		% then this is a 'lower' quantile dependence measure
      out1(jj) = mean((X<=quantile(X,q(jj))).*(Y<=quantile(Y,q(jj))))/q(jj);
   else
      out1(jj) = (1-2*q(jj)+mean((X<=quantile(X,q(jj))).*(Y<=quantile(Y,q(jj)))))/(1-q(jj));
   end
end