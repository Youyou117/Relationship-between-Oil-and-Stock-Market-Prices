function [CL,kappat,ft] = rotgumbel_GAS_CL(theta,data,kappa0,hessINFO)
%function CL = rotgumbel_GAS_CL(theta,data)
% The negative copula log-likelihood of a member of ROTATED Gumbel's family
% (from Joe (1997), p142), where the parameter varies through time according to the "GAS" model for Creal, Koopman and Lucas (2011, JAE)
%
% INPUTS:   theta = [w,a,b], the parameters of the GAS specification
%           data, a Tx2 matrix of Unif(0,1) random variables.
%           kappa0, a scalar, the value of the Rotated Gumbel copula to use as the starting value (perhaps the estimate from a constant version of this model)
%           hessINFO, a kx2 matrix containing the grid of values at which the hessian was evaluated, and the value of the hessian at those points
%
% OUTPUTS:  CL, a scalar, the negative log-likelihood
%
%  Andrew Patton
%
%  7 Sep 2011

%theta'

if sum(isnan(theta))==0
    
    data = 1-data;			% this is the rotation, and below we just use the Gumbel log-likelihood everywhere
    
    T = size(data,1);
    
    ut = -log(data(:,1));
    vt = -log(data(:,2));
    
    w = theta(1);
    a = theta(2);
    b = theta(3);
    
    % will model  f[t] = w + b*f[t-1] + a*DELTA[t-1]*S[t-1]
    % where f[t] = log(kappa[t]-1),  so that kappa[t] = 1+exp{ f[t] }, ensuring kappa is always above 1
    
    h = 0.0001;  % step size to use when computing score
    
    ft = nan(T,1);
    kappat = nan(T,1);
    kappat(1) = kappa0;
    ft(1) = log(kappat(1)-1);
    
    for tt=2:T
        It = interp1(hessINFO(:,1),hessINFO(:,2),kappat(tt-1));     % estimated hessian for kappa
        DELTAt = (-gumbelCL(kappat(tt-1)+h,data(tt-1,:))--gumbelCL(kappat(tt-1),data(tt-1,:)))/h;         % estimated score for kappa. NOTE: my gumbelCL code returns the *neg* log-like, so need to undo that here

        dkappadf = exp(kappat(tt-1));  % used below
        Itildat = It / ( dkappadf^2) ;                         % estimated hessian for f
        
        DELTAtildat = DELTAt / (  dkappadf  );                    % estimated score for f
        
        ft(tt) = w + b*ft(tt-1) + a*DELTAtildat/sqrt(Itildat);  % dkappadf drops out in this univar case (as Creal et al. note) but I keep it here for completeness
        kappat(tt) = 1+exp(ft(tt));
    end
    
    
    CL = log(gumbel_cdf(data(:,1),data(:,2),kappat)) - log(data(:,1)) - log(data(:,2));
    CL = CL + (kappat-1).*(log(ut)+log(vt)) - (2-1./kappat).*(log(ut.^kappat+vt.^kappat));
    CL = CL + log((ut.^kappat + vt.^kappat).^(1./kappat) + kappat - 1);
    CL = sum(CL);
    CL = -CL;
    
    if isnan(CL) || isinf(CL)
        CL = 1e8;
    end
    
else
    CL = 1e7;
end