% The negative copula log-likelihood of a 
% member of Clayton's family
% From Joe(1997), p141
%
% Friday, 29 Sep, 2000
%
% Andrew Patton

% INPUTS: theta 
%				data = [U V];
% 
% OUTPUTS: LL, a Tx1 vector of the log-likelihood evaluated
%		at each point in the sample

function CL = claytonCL(theta,data)

if abs(theta)<=0.00001    % 19oct04: Clayton cop is well defined at theta=0, but it is a limit. setting theta to zero in the following
    theta = 0.00001;     % formula involves 1/0. so just use theta=0.00001 for theta at or near zero.
end
CL = log(data(:,1).*data(:,2))*(1+theta);
CL = CL + (2+1/theta).*log( (data(:,1).^(-theta)) + (data(:,2).^(-theta)) -1);
CL = log(1+theta) - CL;
CL = -CL;