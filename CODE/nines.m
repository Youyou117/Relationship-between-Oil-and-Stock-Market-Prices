function out1 = nines(varargin)
%function out1 = nines(varargin)
%
%   	Returns a matrix with all elements -999.99
%	useful for pre-creating matrices and checking
%	that all elements have been used.
%
%   See also ZEROS, ONES
%
% Modified version of the ZEROS and ONES built-in function
%
%  Andrew Patton
%
%  August 2000

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.5 $  $Date: 1997/11/21 23:29:17 $
%   Built-in function.

if nargin==1	% only one input
   if size(varargin{1},1)==1 & size(varargin{1},2)==1
      out1 = -999.99*ones(varargin{:},1);	% if a scalar, then assume its a column vector
   else
      out1 = -999.99*ones(size(varargin{1}));	% if a matrix, then match its dimensions
   end
elseif nargin==2
   if size(varargin{1},1)==1 
      if size(varargin{2},2)==1
         out1 = -999.99*ones(varargin{:});
      else
         out1 = -999.99*ones(varargin{1},size(varargin{2},2));
      end
   else
      if size(varargin{2},2)==1
         out1 = -999.99*ones(size(varargin{1},1),varargin{2});
      else
         out1 = -999.99*ones(size(varargin{1},1),size(varargin{2},2));
      end
   end
else
   temp = zeros(1,nargin);
   for ii = 1:nargin
      if length(varargin{ii})==1
         temp(ii) = varargin{ii};
      else
         temp(ii) = length(varargin{ii});
      end
   end
   out1 = -999.99*ones(temp);
end