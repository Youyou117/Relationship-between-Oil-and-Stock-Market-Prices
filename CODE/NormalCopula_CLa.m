% The negative copula log-likelihood of a 
% bivariate Normal distribution
% with constant correlation
%
% Monday, 4 Sep, 2000
%
% Andrew Patton

% INPUTS: theta ;
%				data = [U V];

function CLa = NormalCopula_CLa(theta,Zdata)

x = norminv(Zdata(:,1),0,1);
y = norminv(Zdata(:,2),0,1);

CLa = -1*(2*(1-theta^2))^(-1)*(x.^2+y.^2-2*theta*x.*y);
CLa = CLa + 0.5*(x.^2+y.^2);  
CLa = CLa - 1/2*log(1-theta^2);
CLa = -CLa;

