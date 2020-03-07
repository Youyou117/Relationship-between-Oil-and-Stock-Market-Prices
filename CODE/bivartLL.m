function LL = bivartLL(theta,data);
% Computes the negative log-likelihood of
% a standardised bivariate Student's t distribution

mu1 = theta(1);
mu2 = theta(2);
v1  = theta(3);
cov = theta(4);
v2  = theta(5);
nu  = theta(6);
HHAT = [v1,cov;cov,v2];

data = [data(:,1)-mu1,data(:,2)-mu2];
T = size(data,1);
LL = 0.5*2*(log(nu)-log(nu-2)) + gammaln(0.5*(2+nu)) - log(nu*pi) - gammaln(0.5*nu)- 0.5*log(det(HHAT));
for tt = 1:T;
   LL = LL  - 0.5*(2+nu)*log(1+ 1/(nu-2)*data(tt,:)*inv(HHAT)*data(tt,:)');
end
LL = sum(LL);
LL = -LL;
