function out1 = GLquad(f_name,a,b,n,varargin)
%a function out1 = GLquad(f_name,a,b,n,varargin)
%
% Simple function to do Gauss-Legendre quadrature on a univariate function
% over the inteval [a,b]. The function "f_name(x)" must take a vector of
% values x at which it is to be evaluated
%
%  INPUTS:  f_name, a string, the function name
%           a, a scalar, the lower bound of the interval
%           b, a scalar, the upper bound of the interval
%           n, a scalar, the number of "nodes" to use in the quadrature
%           varargin, arguments to pass to the function
%
%  OUTPUT:  out1, a scalar, the estimated integral of the function "f_name" over the interval [a,b]
%
% Andrew Patton
%
% 21 Feb 2011

if nargin<4 || isempty(n)
    n=10;
end
if nargin<3 || isempty(b)
    b=1;
end
if nargin<2 || isempty(a)
    a=-1;
end

if isscalar(n) 
[x,w] = GLNodeWt(n); 
else
    x = n(:,1);
    w = n(:,2);
end

fval = feval(f_name,(x+1)*(b-a)/2+a,varargin{:})*(b-a)/2;  % have to re-scale nodes and function values if integral is over [a,b] rather than [-1,1]
out1 = fval'*w;
