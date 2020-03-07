function out1 = ang_chen1(X,Y,q)
% function out1 = ang_chen1(X,Y,q)
%
% Computes the 'exceedence correlations' discussed in
% Ang and Chen (2001):  Corr[X,Y|X>v,Y>v] if v>0 and X<v,Y<v if v<0
% 
% INPUTS:	X, a Tx1 vector of data
%				Y, a Tx1 vector of data
%				q, a kx1 vector of quantiles to estimate the correl 
%
% OUTPUT:	out1, a qx1 vector of the measure at each quantile
%
% Wednesday, 10 October, 2001
%
% Andrew Patton
%
% see also: quantiledep.m
%


T = size(X,1);
k = size(q,1);
out1 = nines(k+1,1);		% there's an extra one because we do both above and below 0.5

counter = 1;
for jj = 1:k;
   if q(jj)<0.5;		% then this is a 'lower' quantile dependence measure
      temp = find( (X<=quantile(X,q(jj))).*(Y<=quantile(Y,q(jj))) );
      if size(temp,1)<=1
         out1(counter,1) = -999.99;		% not enough obs to compute the correlation
      else
         out1(counter,1) = corrcoef12(X(temp),Y(temp));
      end
      counter=counter+1;
   elseif q(jj)==0.5
      temp = find( (X<=quantile(X,q(jj))).*(Y<=quantile(Y,q(jj))) );
      out1(counter,1) = corrcoef12(X(temp),Y(temp));
      temp = find( (X>quantile(X,q(jj))).*(Y>quantile(Y,q(jj))) );      
      out1(counter+1,1) = corrcoef12(X(temp),Y(temp));
      counter=counter+2;
   else
      temp = find( (X>quantile(X,q(jj))).*(Y>quantile(Y,q(jj))) );      
      if size(temp,1)<=1
         out1(counter,1) = -999.99;% not enough obs to compute the correlation
      else
         out1(counter,1) = corrcoef12(X(temp),Y(temp));
      end
      counter=counter+1;
   end
end