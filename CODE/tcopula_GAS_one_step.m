function [rhot,ft] = tcopula_GAS_one_step(theta,lag_rho,lag_U,rho0,RR,NN,HESSstudt)
%function [rhot,ft] = tcopula_GAS_one_step(theta,lag_rho,lag_U,rho0,RR,NN,HESSstudt)
%
% Function to compute the one-step ahead parameter of the GAS model for a t-coupla with time-varying rho (useful for OOS forecasting studies of this
% model, faster than having to call the likelihood and compute rhot for the whole time series)
%
% INPUTS:   theta = [w,a,b,nuinv], the parameters of the GAS specification, and the *inverse* DoF parameter
%           lag_rho, a scalar, the lagged rho
%           lag_U, a 1x2 vector, the lagged value of [U1,U2]
%           rho0, a scalar, the value of the correlation parameter to use as the starting value (perhaps the estimate from a constant version of this model)
%           RR, a k1x1 vector of values of rho at which the hessian was computed
%           NN, a k2x1 vector of values of nu at which the hessian was computed
%           HESSstudt, a k1 x k2 x 2 x 2 matrix, containing the 2x2 hessian for each combination of values of [rho,nu]
%
% OUTPUTS:  rhot, a scalar, the next-period value of rho
%           ft, a scalar, the next period value of f
%
%  Andrew Patton
%
%  1 Nov 2011
%
% (This is basically a stripped down version of the likelihood function, "tcopula_GAS_CL.m"

% will model  f[t] = w + b*f[t-1] + a*DELTA[t-1]*S[t-1]
% where f[t] = log( (1+rho[t]) / (1-rho[t]) ),  so that rho is always inside (-1,1)

% for numerical stability, will use a modified transformation to keep rho inside (-RBAR,+RBAR),  where RBAR = 0.999

RBAR = 0.9999;  % can make this equal to 1, in theory, but issues occur very close to that boundary

w = theta(1);
a = theta(2);
b = theta(3);
nuinv = theta(4);
nu = 1/theta(4);

h = 0.00001;  % step size to use when computing score

lag_f = log( (RBAR+lag_rho)/(RBAR-lag_rho) );  % used below.

% now computing the one-step ahead value of rho and f
It = interp2(RR',NN',squeeze(HESSstudt(:,:,1,1))',lag_rho,nu,'nearest')    ;  % hessian for rho
DELTAt = (-tcopulaCL([lag_rho+h;nu],lag_U)--tcopulaCL([lag_rho;nu],lag_U))/h    ;         % estimated score for rho. NOTE: my tcopulaCL2 code returns the *neg* log-like, so need to undo that here

drhodf = 2*RBAR*exp(-lag_f)/((1+exp(-lag_f))^2);  % used below
Itildat = It / ( drhodf^2) ;                         % estimated hessian for f
DELTAtildat = DELTAt / (  drhodf  );                    % estimated score for f

Itildat = max(Itildat,1e-6);                            % imposing a min value here to stop the likelihood blowing up when Rho is very close to the boundary
DELTAtildat = max(min(DELTAtildat,1e4),-1e4);           % imposing that this is inside (-1e6,1e6)

ft = w + b*lag_f + a*DELTAtildat/sqrt(Itildat);

ft = max(min(ft,100),-100);                     % imposing that this is inside (-100,100)
rhot = RBAR*(1-exp(-ft))/(1+exp(-ft));

