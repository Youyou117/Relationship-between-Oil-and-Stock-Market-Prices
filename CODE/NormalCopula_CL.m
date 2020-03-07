% The negative copula log-likelihood of a 
% bivariate Normal distribution
% with constant correlation
%
% Monday, 4 Sep, 2000
%
% Andrew Patton

% INPUTS: theta ;
%				data = [U V];

function CL = NormalCopula_CL(theta,Zdata)

x = norminv(Zdata(:,1),0,1);
y = norminv(Zdata(:,2),0,1);

CL = -1*(2*(1-theta^2))^(-1)*(x.^2+y.^2-2*theta*x.*y);
CL = CL + 0.5*(x.^2+y.^2);  
CL = sum(CL) - size(x,1)/2*log(1-theta^2);
CL = -CL;

