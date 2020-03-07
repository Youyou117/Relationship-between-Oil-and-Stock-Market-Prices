function xy = corrcoef(x,y)
%CORRCOEF Correlation coefficients.
%   CORRCOEF(X) is a matrix of correlation coefficients formed
%   from array X whose each row is an observation, and each
%   column is a variable.
%   CORRCOEF(X,Y), where X and Y are column vectors is the same as
%   CORRCOEF([X Y]).
%   
%   If C is the covariance matrix, C = COV(X), then CORRCOEF(X) is
%   the matrix whose (i,j)'th element is
%
%          C(i,j)/SQRT(C(i,i)*C(j,j)).
%
%   See also COV, STD.

%   J. Little 5-5-86
%   Revised 6-9-88 LS 2-13-95 BJ
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.5 $  $Date: 1997/11/21 23:23:34 $

switch nargin
case 1
   c = cov(x);
case 2
   c = cov(x,y);
otherwise
  error('Not enough input arguments.');
end

d = diag(c);
xy = c./sqrt(d*d');
xy = xy(1,2);

