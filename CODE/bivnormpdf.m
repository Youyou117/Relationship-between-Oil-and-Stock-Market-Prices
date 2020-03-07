% Computes the value of the multivariate Normal pdf at a specified point
% 
% INPUTS:	X,   a Tx1 vector
%				Y,		a Tx1 vector
%				MU,  a 1x2 vector of means of each of the x[i]'s in X
%				VCV, a 2x2 variance-covariance matrix
%
% Tuesday, 10 April, 2001.
%
% Andrew Patton

% From Hamilton (1994, p748).

function out1 = bivnormpdf(X,Y,MU,VCV)

k = 2;
detVCV = det(VCV);
invVCV = inv(VCV);

T = max(size(X,1),size(Y,1));
N = max(size(X,2),size(Y,2));	% can only stretch in one direction: choose whichever is bigger

if T>=N
   if size(X,1)<T
      X = X(1,1)*ones(T,1);
   end
   if size(Y,1)<T
      Y = Y(1,1)*ones(T,1);
   end
   X = [X,Y];
   out1 = -999.99*ones(T,1);
   for tt = 1:T;
      out1(tt) = (2*pi)^(-k/2)*detVCV^(-0.5)*exp(-0.5*(X(tt,:)-MU)*invVCV*((X(tt,:)-MU)'));
   end
else
   if size(X,2)<N
      X = X(1,1)*ones(1,N);
   end
   if size(Y,2)<N
      Y = Y(1,1)*ones(1,N);
   end
   X = [X;Y]';
   out1 = -999.99*ones(1,N);
   for tt = 1:N;
      out1(tt) = (2*pi)^(-k/2)*detVCV^(-0.5)*exp(-0.5*(X(tt,:)-MU)*invVCV*((X(tt,:)-MU)'));
   end
end   

