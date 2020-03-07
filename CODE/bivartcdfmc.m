function out1 = bivartcdfmc(b1,b2,rho,nu,nobs);
% Numerically approximates the double integral
% of a bivariate Student's distribution, with zero
% means and unit variances, correlation rho and degrees of freedom nu
% VIA MONTE CARLO
%
% INPUTS:	b1, an Nx1 vector or a scalar, the upper bounds on the outer integral
%				b2, an Nx1 vector or a scalar, the upper bounds on the inner integral
%				rho, an NxM matrix or a 1xM vector, the correlation coefficients of the joint densities
%				nu,  an NxM matrix or a 1xM vector, the degrees of freedom parameter
%				nobs, a scalar: the number of Monte Carlo samples to draw (actually 4 times this effectively)
%						default is 25,000.
%
% OUTPUTS: out1, a NxM matrix of the cdf's of the joint standard Normal dist'ns 
%								at (b1,b2) with correlation coefficient rho
%
%  Andrew Patton
%
% Friday, 27 April, 2001.
%
% corrected to deal with standardisation (thanks to comment from Miguel Segoviano): 18 August, 2003.

if nargin<5
   nobs = 25000;
end


N = max([size(b1,1);size(b2,1);size(rho,1);size(nu,1)]);
M = max(size(rho,2),size(nu,2));

if size(b1,1)~=N;			% stretching the input arguments
   b1 = b1*ones(N,1);
end
if size(b2,1)~=N;
   b2 = b2*ones(N,1);
end
if size(rho,1)~=N;
   rho = ones(N,1)*rho;
end
if size(nu,1)~=N;
   nu = ones(N,1)*nu;
end
if size(rho,2)~=M;
   rho = rho*ones(1,M);
end
if size(nu,2)~=M;
   nu = nu*ones(1,M);
end

out1 = -999.99*ones(N,M);
for n = 1:N;
   for m = 1:M;
      
      % In this step we make use of the fact that the mv t dist'n is symmetric about both
      % the main and secondary diagonals. Thus with each draw we effectively get four
      % observations, by forcing the sample to be symmetric about these two lines.
      % This appears to speed things up by roughly 4 times (as you might guess).
      temp = mvtrnd([1,rho(n,m);rho(n,m),1],nu(n,m),nobs);
      temp = temp*sqrt((nu(n,m)-2)/nu(n,m));											% this is the new line 18aug03
      temp = [temp;-temp;[temp(:,2),temp(:,1)];-[temp(:,2),temp(:,1)]];  
      out1(n,m) = mean((temp(:,1)<=b1(n,m)).*(temp(:,2)<=b2(n,m)));
   end
end

   
