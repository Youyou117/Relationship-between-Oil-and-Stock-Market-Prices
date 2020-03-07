function [kappat,ft] = rotgumbel_GAS_one_step(theta,lag_k,lag_U,kappa0,hessINFO)
%function CL = rotgumbel_GAS_CL(theta,data)
%
% Function to compute the one-step ahead parameter of the GAS model for a Rotated Gumbel with time-varying kappa (useful for OOS forecasting studies of this
% model, faster than having to call the likelihood and compute rhot for the whole time series)
%
% INPUTS:   theta = [w,a,b], the parameters of the GAS specification
%           lag_k, a scalar, the lagged kappa
%           lag_U, a 1x2 vector, the lagged value of [U1,U2]
%           kappa0, a scalar, the value of the Rotated Gumbel copula to use as the starting value (perhaps the estimate from a constant version of this model)
%           hessINFO, a kx2 matrix containing the grid of values at which the hessian was evaluated, and the value of the hessian at those points
%
% OUTPUTS:  kappat, a scalar, the next-period value of kappa
%           ft, a scalar, the next period value of f
%
%  Andrew Patton
%
%  1 Nov 2011
%
% (This is basically a stripped down version of the likelihood function, "Rotgumbel_GAS_CL.m"


lag_U = 1-lag_U;			% this is the rotation, and below we just use the Gumbel log-likelihood everywhere

w = theta(1);
a = theta(2);
b = theta(3);

% will model  f[t] = w + b*f[t-1] + a*DELTA[t-1]*S[t-1]
% where f[t] = log(kappa[t]-1),  so that kappa[t] = 1+exp{ f[t] }, ensuring kappa is always above 1

h = 0.0001;  % step size to use when computing score

lag_f = log(lag_k -1);


It = interp1(hessINFO(:,1),hessINFO(:,2),lag_k);     % estimated hessian for kappa
DELTAt = (-gumbelCL(lag_k  +h,lag_U)--gumbelCL(lag_k,lag_U))/h;         % estimated score for kappa. NOTE: my gumbelCL code returns the *neg* log-like, so need to undo that here

dkappadf = exp(lag_k);  % used below
Itildat = It / ( dkappadf^2) ;                         % estimated hessian for f

DELTAtildat = DELTAt / (  dkappadf  );                    % estimated score for f

ft = w + b*lag_f + a*DELTAtildat/sqrt(Itildat);  % dkappadf drops out in this univar case (as Creal et al. note) but I keep it here for completeness
kappat = 1+exp(ft);
