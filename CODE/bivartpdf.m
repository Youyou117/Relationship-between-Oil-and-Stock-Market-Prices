% Computes the value of the bivariate standardised Student's t pdf at a specified point
% 
% INPUTS:	X1,  a Tx1 vector or a scalar
%				X2,  a Tx1 vector or a scalar
%				RHO, a Tx1 vector or a scaler, correlation coefficient
%				NU, 	a Tx1 vector or a scalar, the degrees of freedom parameter
%
% Thursday, 26 April, 2001.
%
% Andrew Patton

% working OK (but double-check if any weird results are obtained)
function out1 = bivartpdf(X1,X2,RHO,NU)

T = max(size(X1,1),max(size(X2,1),max(size(RHO,1),size(NU,1))));
if size(X1,1)~=T;
   X1 = ones(T,1)*X1;
end
if size(X2,1)~=T;
   X2 = ones(T,1)*X2;
end
if size(RHO,1)~=T;
   RHO = ones(T,1)*RHO;
end
if size(NU,1)~=T;
   NU = ones(T,1)*NU;
end

out1 = -999.99*ones(T,1);
for tt = 1:T;
   if NU(tt)<300
      out1(tt) = gamma((NU(tt)+2)/2)/gamma(NU(tt)/2)*((NU(tt)*pi)^(-1))/sqrt(1-RHO(tt)^2)*((1 + (X1(tt)^2+X2(tt)^2-2*RHO(tt)*X1(tt)*X2(tt))/(NU(tt)*(1-RHO(tt)^2)))^(-(NU(tt)+2)/2));
   else
      out1(tt) = bivnormpdf(X1(tt),X2(tt),[0,0],[1,RHO(tt);RHO(tt),1]);
   end
end
   