function [data,rhot,ft] = tcopula_GAS_rnd(theta,T,rho0,RR,NN,HESSstudt)
%function [data,rhot,ft] = tcopula_GAS_rnd(theta,T,rho0,RR,NN,HESSstudt)
%
% Function to simulate from a BIVARIATE Student's t copula,
% where the CORRELATION parameter varies through time according to the "GAS" model for Creal, Koopman and Lucas (2011, JAE)
% but the DEGREES OF FREEDOM parameter is constant.
%
% INPUTS:   theta = [w,a,b,nuinv], the parameters of the GAS specification, and the *inverse* DoF parameter
%           T, a scalar, the number of observations to simulate
%           rho0, a scalar, the value of the correlation parameter to use as the starting value (perhaps the estimate from a constant version of this model)
%           RR, a k1x1 vector of values of rho at which the hessian was computed
%           NN, a k2x1 vector of values of nu at which the hessian was computed
%           HESSstudt, a k1 x k2 x 2 x 2 matrix, containing the 2x2 hessian for each combination of values of [rho,nu]
%
% OUTPUTS:  data, a Tx2 matrix of simulated Unif(0,1) data
%           rhot, a Tx1 vector, the time series of the correlation parameter
%           ft, a Tx1 vector, the time series of the transformed copula parameter
%
%  Andrew Patton
%
%  8 Sep 2011

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

% generating the time series of rho
ft = nan(T,1);
rhot = nan(T,1);
data = nan(T,2);

rhot(1) = rho0;
ft(1) = log( (RBAR+rhot(1))/(RBAR-rhot(1)) );
data(1,:) = tdis_cdf(mvtrnd([[1,rhot(1)];[rhot(1),1]],nu),nu);

for tt=2:T
    It = interp2(RR',NN',squeeze(HESSstudt(:,:,1,1))',rhot(tt-1),nu,'nearest')    ;  % hessian for rho
    DELTAt = (-tcopulaCL([rhot(tt-1)+h;nu],data(tt-1,:))--tcopulaCL([rhot(tt-1);nu],data(tt-1,:)))/h    ;         % estimated score for rho. NOTE: my tcopulaCL2 code returns the *neg* log-like, so need to undo that here
    
    drhodf = 2*RBAR*exp(-ft(tt-1))/((1+exp(-ft(tt-1)))^2);  % used below
    Itildat = It / ( drhodf^2) ;                         % estimated hessian for f
    DELTAtildat = DELTAt / (  drhodf  );                    % estimated score for f
    
    ft(tt) = w + b*ft(tt-1) + a*DELTAtildat/sqrt(Itildat);
    
    rhot(tt) = RBAR*(1-exp(-ft(tt)))/(1+exp(-ft(tt)));
    data(tt,:) = tdis_cdf(mvtrnd([[1,rhot(tt)];[rhot(tt),1]],nu),nu);
    
end

