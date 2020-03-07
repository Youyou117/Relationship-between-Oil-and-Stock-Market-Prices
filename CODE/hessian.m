function H = hessian(f,x,varargin)
% PURPOSE: Computes finite difference Hessian
% -------------------------------------------------------
% Usage:  H = hessian(func,x,varargin)
% Where: func = function name, fval = func(x,varargin)
%           x = vector of parameters (n x 1)
%    varargin = optional arguments passed to the function
% -------------------------------------------------------
% RETURNS:
%           H = finite differnce hessian
% -------------------------------------------------------

% Code from:
% COMPECON toolbox [www4.ncsu.edu/~pfackler]
% documentation modified to fit the format of the Ecoometrics Toolbox
% by James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

%eps = 1e-5;   % AJP commented out 5feb01, this is the mistake that Kevin found.

n = size(x,1);
fx = feval(f,x,varargin{:});

% Compute the stepsize (h)
%h = (1e-5).^(1/3)*max(abs(x),1e-2);
%h = eps.^(1/3)*max(abs(x),1e-4);	% trying to get rid of neg on diag of inverse hessian (4may01)
h = -eps.^(1/3)*max(abs(x),1e-2);
h = -eps.^(1/3)*max(abs(x),1e-4);
%h = 0.001;
%sh = -0.0001;
h;


xh = x+h;
h = xh-x;
ee = sparse(1:n,1:n,h,n,n);

% Compute forward step
g = zeros(n,1);
for i=1:n
    g(i) = feval(f,x+ee(:,i),varargin{:});
end

H=h*h';
% Compute "double" forward step
for i=1:n
    for j=i:n
        H(i,j) = (feval(f,x+ee(:,i)+ee(:,j),varargin{:})-g(i)-g(j)+fx)/H(i,j);
        H(j,i) = H(i,j);
    end
end
