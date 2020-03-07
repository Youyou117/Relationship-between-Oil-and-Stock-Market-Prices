function out1 = empiricalCDF(data,xx);
% Computes empirical probability integral transform
%
%  INPUTS: 	data, a Txk matrix of data (organised in columns)
%					xx, (optional) a qxk or qx1 vector of points in support of data at which to compute the empricial CDF
%
%  OUTPUTS: a Txk (or qxk if xx entered) matrix of probability integral transforms of
%              the data, using the empirical CDF
%
%  Andrew Patton
%
%  9 May, 2001.

[T,k] = size(data);
if nargin<2
   out1 = -999.99*ones(T,k);
   for jj = 1:k
      temp = [data(:,jj),(1:1:T)'];
      temp2 = sortrows(temp,1);	
      temp3 = [temp2,(1:1:T)'/(T+1)];		% dividing by T+1 rather than T so that none equals 1
      temp4 = sortrows(temp3,2);
      out1(:,jj) = temp4(:,3);
   end
else
   [q,q2] = size(xx);
   if q2==1
      xx = xx*ones(1,k);
   end
   out1 = -999.99*ones(q,k);
   for jj = 1:k
      for ii = 1:q
         out1(ii,jj) = mean(data(:,jj)<=xx(ii,jj));
      end
   end
end
   