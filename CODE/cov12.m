function c = cov12(x,y)
% covariance coefficients.
%   COV(X) is the covariance between the two
%   series (each series is a column vector)
%   
%   See also COV, STD, CORRCOEF.

%   J. Little 5-5-86
%   Revised 6-9-88 LS 2-13-95 BJ
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.5 $  $Date: 1997/11/21 23:23:34 $
%
% 	  Modified by Andrew Patton
% 	  27 April, 2001.


switch nargin
case 1
   c = cov(x);
   c = c(1,2);
case 2
   c = cov(x,y);
   c = c(1,2);
otherwise
  error('Not enough input arguments.');
end
