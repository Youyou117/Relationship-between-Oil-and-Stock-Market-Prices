function out1 = bivarnormcdf(b1,b2,rho);
% Numerically approximates the double integral
% of a bivariate Normal distribution, with zero
% means and variance-covariance matrix: [ 1, rho ; rho , 1]
%
% INPUTS:	b1, an Nx1 vector or a scalar, the upper bounds on the outer integral
%				b2, an Nx1 vector or a scalar, the upper bounds on the inner integral
%				rho, an NxM matrix or a 1xM vector, the correlation coefficients of the joint densities
%
% OUTPUTS: out1, a NxM matrix of the cdf's of the joint standard Normal dist'ns 
%								at (b1,b2) with correlation coefficient rho
%
%  Andrew Patton
%
% Thursday, 14 December, 2000
%
% Compared with bvncdf.src code in GAUSS: seems OK (max difference ~1e-06) which
% is about the error that 'quadl.m' is aiming for.
% (Much slower than Gauss code, though...)
% 
% see bivarnormcdf2.m if want values in extremes (<-5)

if sum(size(b1)==size(b2))+sum(size(b1)==size(rho))~=4	% then they're not all the same size
   
   
   N = max(size(b1,1),size(b2,1));
   M = size(rho,2);
   
   %stretching the input arguments
   if size(b1,1)~=N;
      b1 = b1*ones(N,1);
   end
   if size(b2,1)~=N;
      b2 = b2*ones(N,1);
   end
   if size(rho,1)~=N;
      rho = ones(N,1)*rho;
   end
   
%   b1
%   b2
%   rho
   
   out1 = -999.99*ones(N,M);
   for n = 1:N;
      for m = 1:M;
         if b1(n)>0.01
            out1(n,m) = quadl('bivarnormcdf_arg',-5,b1(n),[],[],b2(n),rho(n,m));
         else
            out1(n,m) = quadl('bivarnormcdf_arg',-10,b1(n),[],[],b2(n),rho(n,m));
         end
      end
   end
else
   N = size(b1,1);
   M = size(b2,2);
   out1 = -999.99*ones(N,M);
   for n = 1:N;
      for m = 1:M;
         if b1(n,m)>0.01
            out1(n,m) = quadl('bivarnormcdf_arg',-5,b1(n,m),[],[],b2(n,m),rho(n,m));
         else
            out1(n,m) = quadl('bivarnormcdf_arg',-10,b1(n,m),[],[],b2(n,m),rho(n,m));
         end
      end
   end
end
   
