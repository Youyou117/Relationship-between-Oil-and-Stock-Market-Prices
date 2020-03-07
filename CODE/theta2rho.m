function out1 = theta2rho(theta)
%function out1 = theta2rho(theta)
%
% Takes a k*(k-1)/2 x 1 vector and transforms it into
% a kxk correlation matrix
%
% INPUTS: 	theta, a kx(k-1)/2 x 1 vector of correlations [rho12,rho13,...rho(k-1)k]
%
% OUTPUTS:	out1, a kxk correlation matrix
%
%  Andrew Patton
%
%  Monday, 5 August, 2002.

m = length(theta);
k = (1 + sqrt(1+8*m))/2;
if mod(k,1)~=0;	% then k is not an integer
   'vector is not of correct length'
end
out1 = nines(k,k);
counter=1;
for ii=1:k;
   for jj=ii:k;
      if ii==jj;
         out1(ii,jj,1)=1;
      else
         out1(ii,jj,1)=theta(counter);
         out1(jj,ii,1)=theta(counter);
         counter = counter+1;
      end
   end
end

