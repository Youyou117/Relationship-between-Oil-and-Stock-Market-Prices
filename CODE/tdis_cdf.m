function FF = tdis_cdf (x, n)
% PURPOSE: returns cdf at x of the t(n) distribution
%---------------------------------------------------
% USAGE: cdf = tdis_cdf(x,n)
% where: x = a maxtrix (modified to allow for matrix inputs: AJP, 26aug03)
%        n = a scalar parameter with dof
%---------------------------------------------------
% RETURNS:
%        a vector of cdf at each element of x of the t(n) distribution      
% --------------------------------------------------
% SEE ALSO: tdis_inv, tdis_rnd, tdis_pdf, tdis_prb
%---------------------------------------------------

%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg
% modified by J.P. LeSage

if nargin ~= 2
error('Wrong # of arguments to tdis_cdf');
end;

if any(any(n<=0))
   error('tdis_cdf dof is wrong');
end
[nobs k] = size(x);
FF = [];
data=x;
for ii=1:k;
   x=data(:,ii);
   neg = x<0;
%   F = fdis_cdf(x.^2,1,n);
   F = fdis_cdf(x.^2,ones(size(x,1),size(x,2)),n.*ones(size(x,1),size(x,2)));  % 26jan08: modified by AJP to allow for vector inputs of n
   iota = ones(nobs,1);
   out = iota-(iota-F)/2;
   F = out + (iota-2*out).*neg;
   FF = [FF,F];
end

   
   
   