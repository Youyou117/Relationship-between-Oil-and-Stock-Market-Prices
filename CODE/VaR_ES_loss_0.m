function [Eloss, loss] = VaR_ES_loss_0(ve,y,alpha)
% function loss = VaR_ES_loss_0(ve,y,alpha)
%
%  Function to compute the FZ loss function value. This version is the
%  homogeneous of degree zero version. Only works when VaR and ES are strictly negative.
%
%       G1(x) = 0 
%       G2(x)= -1/x
%       caliG2(x) = -log(-x)
%
% INPUTS:   [ve], a 2x1 vector or a Tx2 matrix, = [VaR, ES] forecasts 
%           y, a scalar or a Tx1 vector, realizations of the target variable
%           alpha, a scalar or a kx1 vector, in (0,1), the probability level for the VaR and ES forecasts
%
%  OUTPUTS: Eloss, a scalar or a kx1 vector, the mean of "loss"
%           loss, a scalar or a Tx1 vector, values of the loss function at all points in time
%
%  Andrew Patton
%
%  7 November 2015
%
% This code was used in: 
% Patton, A.J., J.F. Ziegel, and R. Chen, 2018, Dynamic Semiparametric
% Models for Expected Shortfall (and Value-at-Risk), Journal of Econometrics, forthcoming. 

% Loss function general form is:
% L(y,v,e;a,G1,G2) = (y<=v).*(G1(v)-G1(y))  +  1/a*G2(e)*(y<=v)*(v-y)  +  G2(e)*(e-v) - caliG2(e)

if max(size(ve))==2  
    v = ve(1);
    e = ve(2);
else
    v = ve(:,1);
    e = ve(:,2);
end

T = max([size(y,1),size(ve,1)]);
k = length(alpha);

if size(y,1)<T
    y = y(1)*ones(T,1);
end
if size(v,1)<T
    v = v(1)*ones(T,1);
end
if size(e,1)<T
    e = e(1)*ones(T,1);
end

G1v = 0;
G1y = 0;
G2 = -1./e;
caliG2 = -log(-e);

loss = nan(T,k);
for aa=1:length(alpha)
    a = alpha(aa);
    loss(:,aa) = ( (y<=v)-a ).*(G1v-G1y) ...
        + 1/a*G2.*(y<=v).*(v-y) ...
        + G2.*(e-v) - caliG2;
end
Eloss = mean(loss)';
