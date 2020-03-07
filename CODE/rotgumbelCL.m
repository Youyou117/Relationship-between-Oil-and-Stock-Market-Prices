% The negative copula log-likelihood of a 
% member of ROTATED Gumbel's family
% Taken from Joe (1997), p142.
%
% Thursday, 27 July, 2000
%
% Andrew Patton

% INPUTS: theta ;
%				data = [U V];

function CL = rotgumbelCL(theta,data)

data = 1-data;			% this is the rotation.


ut = -log(data(:,1));
vt = -log(data(:,2));

CL = log(gumbel_cdf(data(:,1),data(:,2),theta)) - log(data(:,1)) - log(data(:,2));
CL = CL + (theta-1)*(log(ut)+log(vt)) - (2-1/theta)*(log(ut.^theta+vt.^theta));
CL = CL + log((ut.^theta + vt.^theta).^(1/theta) + theta - 1);
CL = sum(CL);
CL = -CL;
