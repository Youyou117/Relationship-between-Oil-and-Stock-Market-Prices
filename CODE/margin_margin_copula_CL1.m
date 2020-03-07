% This function returns the negative copula likelihood
% of a general bivariate distribution constructed
% using copula theory
%
% INPUTS: theta = [phi;gamma;kappa], the parameters of the two margins and the copula
%				depvar = [X,Y] a Tx2 matrix of data
%				margin1_name = a string containing the name of the first marginal distribution
%				margin2_name = a string containing the name of the second marginal distribution
%				copula_name  = a string containing the name of the copula
%				pp1, a scalar indicating the number of parameters in phi
%				qq1, a scalar indicating the number of parameters in gamma
%				pp2, a scalar indicating the number of arguments to pass to margin1_name
%				qq2, a scalar indicating the number of arguments to pass to margin2_name
%				varargin, arguments to pass to margin1_name, margin2_name and copula_name
%
%  OUTPUTS: CL, the negative copula likelihood
%
%  Andrew Patton
%
%  Friday, 21 Sep, 2001

% this code is used for finding the bottom rows of the hessian. 

function CL = margin_margin_copula_CL1(theta,depvar,margin1_name,margin2_name,copula_name,pp1,qq1,pp2,qq2,varargin);

phi = theta(1:pp1);
gam = theta(pp1+1:pp1+qq1);
kap = theta(pp1+qq1+1:end);

[ tempL1 temp temp uu] = feval(margin1_name,phi,depvar(:,1),varargin{1:pp2});
[ tempL2 temp temp vv] = feval(margin2_name,gam,depvar(:,2),varargin{pp2+1:pp2+qq2});

%[tempL1,tempL2]
%[[min(uu),mean(uu),max(uu)];[min(vv),mean(vv),max(vv)]]

if size(uu,1)>size(vv,1);		% matching up uu and vv, assuming that if different, then one started later than other
   uu = uu(size(uu,1)-size(vv,1)+1:size(uu,1));
else
   vv = vv(size(vv,1)-size(uu,1)+1:size(vv,1));
end

nn = nargin - 9;
CL = feval(copula_name,kap,[uu,vv],varargin{pp2+qq2+1:nn});

