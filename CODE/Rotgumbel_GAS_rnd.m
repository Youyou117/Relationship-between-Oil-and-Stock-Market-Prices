function [data,kappat,ft] = Rotgumbel_GAS_rnd(theta,T,kappa0,hessINFO)
%function [data,kappat,ft] = Rotgumbel_GAS_rnd(theta,T,kappa0,hessINFO)
% 
% Function to simulate from a BIVARIATE member of ROTATED Gumbel's family
% (from Joe (1997), p142), where the parameter varies through time according to the "GAS" model for Creal, Koopman and Lucas (2011, JAE)
%
% INPUTS:   theta = [w,a,b], the parameters of the GAS specification
%           T, a scalar, the number of observations to simulate
%           kappa0, a scalar, the value of the Rotated Gumbel copula to use as the starting value (perhaps the estimate from a constant version of this model)
%           hessINFO, a kx2 matrix containing the grid of values at which the hessian was evaluated, and the value of the hessian at those points
%
% OUTPUTS:  data, a Tx2 matrix of simulated Unif(0,1) data
%           kappat, a Tx1 vector, the time series of the copula parameter
%           ft, a Tx1 vector, the time series of the transformed copula parameter
%
%  Andrew Patton
%
%  15 Sep 2011

% Below I will simulate the data as a time-varying GUMBEL, and then "rotate" it at the end

% will model  f[t] = w + b*f[t-1] + a*DELTA[t-1]*S[t-1]
% where f[t] = log(kappa[t]-1),  so that kappa[t] = 1+exp{ f[t] }, ensuring kappa is always above 1


w = theta(1);
a = theta(2);
b = theta(3);

h = 0.0001;  % step size to use when computing score

ft = nan(T,1);
kappat = nan(T,1);
data = nan(T,2);

kappat(1) = kappa0;
ft(1) = log(kappat(1)-1);
data(1,:) = Gumbel_rnd(kappat(1),1);

for tt=2:T
    It = interp1(hessINFO(:,1),hessINFO(:,2),kappat(tt-1));     % estimated hessian for kappa
    DELTAt = (-gumbelCL(kappat(tt-1)+h,data(tt-1,:))--gumbelCL(kappat(tt-1),data(tt-1,:)))/h;         % estimated score for kappa. NOTE: my gumbelCL code returns the *neg* log-like, so need to undo that here
    
    dkappadf = exp(kappat(tt-1));  % used below
    Itildat = It / ( dkappadf^2) ;                         % estimated hessian for f
    
    DELTAtildat = DELTAt / (  dkappadf  );                    % estimated score for f
    
    ft(tt) = w + b*ft(tt-1) + a*DELTAtildat/sqrt(Itildat);  % dkappadf drops out in this univar case (as Creal et al. note) but I keep it here for completeness
    kappat(tt) = 1+exp(ft(tt));
    data(tt,:) = Gumbel_rnd(kappat(tt),1);
end

data = 1-data;			% this is the rotation so that this data comes from the Rotated Gumbel copula
