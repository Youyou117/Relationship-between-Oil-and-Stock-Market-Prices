% The negative copula log-likelihood of a 
% member of Frank's family
% Taken from Joe (1997), p141.
%
% Thursday, 28 Sep, 2000
%
% Andrew Patton

% INPUTS: theta ;
%				data = [U V];

function CL = frankCL(theta,data)

eta = 1-exp(-theta);

CL = 2*log(eta - (1-exp(-theta*data(:,1))).*(1-exp(-theta*data(:,2))));
CL = CL + theta*(data(:,1)+data(:,2));
CL = size(data,1)*log(theta*eta) - sum(CL);
CL = -CL;
