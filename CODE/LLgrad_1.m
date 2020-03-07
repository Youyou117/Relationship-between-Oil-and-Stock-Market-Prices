% This program calculates the scores of the 
% a given function
%
% INPUTS: theta, a kx1 vector of parameters
%				function_name, a string containing the name of a function that returns the 
%					value of the log-likelihood at each point in time
%				varargin, arguments to pass to function_name
%
% OUTPUTS; grad, a Txk matrix of scores at each point in time
%
% Wednesday, 19 September, 2001.
%
%  Andrew Patton

function grad = LLgrad_1(function_name,theta,varargin);

k = size(theta,1);
temp = feval(function_name,theta,varargin{:});
T = size(temp,1);
B1 = -999.99*ones(T,k+1);
B2 = -999.99*ones(T,k);
B1(:,1) = temp;
eye1 = eye(k);
for jj = 1:k;
   h = 2.2204e-016^(1/3)*max(abs(theta(jj)),0.01);		% following the same idea in the 'hessian.m' code
   thetajj = theta - h*eye1(:,jj);
   B1(:,jj+1) = feval(function_name,thetajj,varargin{:});
   B2(:,jj) = (B1(:,1)-B1(:,jj+1))/h;
end
grad=B2;