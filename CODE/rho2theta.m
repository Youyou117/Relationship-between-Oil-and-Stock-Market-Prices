function out1 = rho2theta(rho)
%function out1 = rho2theta(rho)
%
% Takes a kxk correlation matrix and returns
% a k*(k-1)/2 x 1 vector of the correlation coeffs
%
% INPUTS: 	rho, a kxk correlation matrix
%
% OUTPUTS:	out1, a kx(k-1)/2 x 1 vector of correlations [rho12,rho13,...rho(k-1)k]
%
%  Andrew Patton
%
%  Monday, 5 August, 2002.


k = size(rho,2);
out1 = nines(k*(k-1)/2,1);
counter=1;
for ii=1:k;
   for jj=ii+1:k
      out1(counter) = rho(ii,jj);
      counter=counter+1;
   end
end
