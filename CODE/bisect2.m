function [x,fval] = bisect(func,a,b,info,varargin)
% PURPOSE: Uses bisection to find root x of a univariate function on [a,b]
%          to tolerance tol, provided f(a) and f(b) have different signs
% ------------------------------------------------
% Usage:    [x,fval] = bisect(func,a,b,info,varargin)
%        or [x,fval] = bisect(func,a,b,[],varargin)
% Where : func = function of the form fval = func(x)
%            a = left endpoint of interval
%            b = right endpoint of interval
% info is a structure variable with:
%         .maxit = maximum # of iterations (default = 100)
%           .tol = convergence tolerance (default = sqrt(eps))   
%        .pflag  = 1 for printing during search (default = 0)    
%    varargin = arguments passed to func()
% ------------------------------------------------
% RETURNS:
%        x  = value at the root of func
%     fval  = function value at the root
%
% MODIFIED: 15 November, 2000 by Andrew Patton
% Changing the code so that it doesn't abort if
% the function has the same sign at the end points
% i.e., if I haven't chosen a wide enough interval
% to bound the zero...
 
% set defaults
maxit = 100;
tol = sqrt(eps);
pflag = 0;
if length(info) > 0
  if ~isstruct(info)
    error('newton: options should be in a structure variable');
  end;
% parse options
fields = fieldnames(info);
nf = length(fields); ycheck = 0; xcheck = 0;
  for i=1:nf
    if strcmp(fields{i},'maxit')
        maxit = info.maxit; 
    elseif strcmp(fields{i},'tol')
        tol = info.tol;
    elseif strcmp(fields{i},'pflag')
        pflag = info.pflag;
    end;
  end;
else % case where we have no options
% no input option, so we use default options

end; % end of case where we have options

% these are for the case of printing intermediate results during iteration
input.cnames = strvcat('function value','converge','argument');


if a>b, error('b must exceed a'), end
sa = sign(feval(func,a,varargin{:}));
sb = sign(feval(func,b,varargin{:}));
%if sa==sb, error ('f has same sign at endpoints'), end   % original code

%modification:
if sa==sb
   if sa==-1   % then the zero is above the upper end-point
      x = b;	 % setting sol'n to upper end-point
      fval = +999.99;		% this is my 'error-code' telling me that it hasn't converged
   else
      x = a;
      fval = -999.99;		% this is my 'error-code' telling me that it hasn't converged
   end
else	% everything is OK so use original code
   x = (a+b)/2;
   d = (b-a)/2;
   while d>tol
      d = d/2;
      if sa == sign(feval(func,x,varargin{:}))
         x = x+d;
      else
         x = x-d;
      end 
      if pflag == 1
         mprint([sa d x],input);
      end;
   end
   fval = sa;
end