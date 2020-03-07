function [Eloss,VEhat,loss,VCV,outT] = garch_FZ_LL(theta,data,alpha,tau,omega,h0,cT)
%function [Eloss,VEhat,loss,VCV,outT] = garch_FZ_LL(theta,data,alpha,tau,omega,h0,cT)
%
% Function to compute the average FZ loss of the GARCH model:
%
%  y(t) = sig(t)*eta(t)
%  sig2(t) = omega + beta*sig2(t-1) + gamma*(y(t-1)^2)
%   e(t) = b*sig(t), for b<0  (since alpha<<0.5)
%   v(t) = c*e(t), for 0<c<1
%
% ** Note: omega is not identified in this model. If it is not provided as an input, it will be set to 1. 
%
%  INPUTS:  theta = [beta,gamma,b,c]', a vector of parameters. We set beta=normcdf(beta) and gamma = exp(gamma) to impose that beta lies in (0,1) and that alpha is positive.
%           data, a Tx1 vector of returns (assumed in this model to be zero mean, although that is not required)
%           alpha, a scalar inside (0,1), the probability level
%           tau, a scalar, a positive number if want to smooth the hit variable (eg, tau=10), use -1 (default) if want NO smoothing and want to use the true objective function
%                   tau=5 is quite a bit of smoothing, tau=10 is some, tau=20 is not much. tau=50 is almost none.
%           omega, a positive scalar, the fixed value for omega to use.  Default=1.
%           h0, starting value for variance series. Default=cov(data)
%           cT, a scalar, the bandwidth parameter to use in computing the VCV matrix. Default cT=T^(-1/3)
%
%  OUTPUTS: Eloss, a scalar, the average FZ loss
%           VEhat, a Tx3 vector, the model-implied [VaR, ES, Variance] series
%           loss, a Tx1 vector, the FZ loss values at each point in time
%           VCV, a 4x4 matrix, the asymptotic covariance matrix of the estimated parameters (assuming that we are evaluating the function at a local min)
%           outT, a 4x3 matrix, the transformed parameters, std errors for the transformed parameters, and t-stats
%
%  Andrew Patton
%
%  4 January 2017
%
% This code was used in: 
% Patton, A.J., J.F. Ziegel, and R. Chen, 2018, Dynamic Semiparametric
% Models for Expected Shortfall (and Value-at-Risk), Journal of Econometrics, forthcoming. 


beta = normcdf(theta(1));
gamma = exp(theta(2));
b = -exp(theta(3));
c = normcdf(theta(4));
a = c*b;
T = length(data);

if nargin<4 || isempty(tau)
    tau = -1;  % default option is to NOT use any smoothing of the indicator function
end
if nargin<5 || isempty(omega)
    omega = 1;  % use the true value of omega here if this is for a simulation, so that [b,c] are centered on the right values. If wrong omega is used the fit will be identical, but we cannot compare bias and efficiency so easily.
end
if nargin<6 || isempty(h0)
    h0 = cov(data);
end
if nargin<7 || isempty(cT)
    cT = T^(-1/3);
end


hhat = nan(T,1);
LLt = nan(T,1);
VEhat = nan(T,2);

hhat(1) = h0;
VEhat(1,:) = [a,b]*sqrt(hhat(1));

if tau==-1
    hitS = (data(1)<VEhat(1,1));  % then use the actual indicator function
else
    hitS = 1./(1+exp( tau*(data(1)-VEhat(1,1)) ));
end

loss = zeros(T,1);
loss(1) = -1/alpha/VEhat(1,2)*hitS*(VEhat(1,1)-data(1)) - 1/VEhat(1,2)*(VEhat(1,2)-VEhat(1,1)) + log(-VEhat(1,2));
for tt=2:T
    hitSL = hitS;
    hhat(tt) = omega + beta*hhat(tt-1) + gamma*(data(tt-1)^2);
    VEhat(tt,:) = [a,b]*sqrt(hhat(tt));
    
    if tau==-1
        hitS = (data(tt)<VEhat(tt,1));  % then use the actual indicator function
    else
        hitS =  1./(1+exp( tau*(data(tt)-VEhat(tt,1)) ));
    end
    loss(tt) = -1/alpha/VEhat(tt,2)*hitS*(VEhat(tt,1)-data(tt)) - 1/VEhat(tt,2)*(VEhat(tt,2)-VEhat(tt,1)) + log(-VEhat(tt,2)) ;
end

Eloss = mean(loss);
if sum(VEhat(:,1)<VEhat(:,2))>0 || max(max(VEhat))>0 || min(VEhat(:,2))<5*min(data) % require ES<VaR, and since I'm using Hom deg 0 loss, also require both to be below 0. Also, don't want ES (and VaR) to go crazy, so limit how far it can be below observed min
    Eloss = 1e6;
end
if sum(sum(isnan(VEhat)))>0 || ~isreal(VEhat) || isinf(Eloss) || ~isreal(theta) || sum(data<=VEhat(:,1))==0  % adding rule that we must see *some* hits, otherwise parameter vector is not identified
    Eloss = 1e7;
end
VEhat(:,3) = hhat;  

if nargout>3  % then compute the VCV matrix
    p = length(theta);  % number of parameters in this model
    
    % first: the derivative of L w.r.t V and E. Same formula for all models for (V,E).
    dLdVE = nan(T,2);      % Note:in this model lamE (the RHS variable in the GAS model) is *not* the same as dL/dE, as we reparameterized the forcing variable. So compute it here.
    dLdVE (:,1) = -1/alpha./VEhat(:,2).*(data<=VEhat(:,1) - alpha);  % dL/dV
    dLdVE (:,2) = +1/alpha./(VEhat(:,2).^2).*(data<=VEhat(:,1)).*(VEhat(:,1)-data) - VEhat(:,1)./(VEhat(:,2).^2) + 1./VEhat(:,2);   % dL/dE
    
    % second: deriv of V and E w.r.t theta, which depends on model
    delV = zeros(T,p);
    delE = zeros(T,p);
    delF = zeros(T,p);  % delV and delE are both driven by delF, so will create this here too
    
    delF(1,1) = omega/((1-beta-gamma)^2);  % getting impact on f(1) from the assumed starting value for f
    delF(1,2) = omega/((1-beta-gamma)^2);  % getting impact on f(1) from the assumed starting value for f
    delF(1,3) = 0;
    delF(1,4) = 0;
    for tt=2:T
        delF(tt,1) = VEhat(tt-1,3) + beta*delF(tt-1,1) ;  % dF/dbeta
        delF(tt,2) = data(tt-1)^2     + beta*delF(tt-1,2) ;  % dF/dgamma
        delF(tt,3) = 0;  % dF/db
        delF(tt,4) = 0;  % dF/dc
    end
    delV(:,1) = ( VEhat(:,1)/2./VEhat(:,3).*delF(:,1) )    *normpdf(theta(1));  % dV/dtheta1 = dV/dbeta * dbeta/dtheta1
    delV(:,2) = ( VEhat(:,1)/2./VEhat(:,3).*delF(:,2) )    *(+exp(theta(2)));   % dV/dtheta2 = dV/dgamma * dgamma/dtheta2
    delV(:,3) = ( VEhat(:,1)/2./VEhat(:,3).*delF(:,3)+1/b )*(-exp(theta(3)));   % dV/dtheta3 = dV/db * db/dtheta3
    delV(:,4) = ( VEhat(:,1)/2./VEhat(:,3).*delF(:,1)+1/c )*normpdf(theta(4));  % dV/dtheta4 = dV/dc * dc/dtheta4
    
    delE(:,1) = ( VEhat(:,2)/2./VEhat(:,3).*delF(:,1) )    *normpdf(theta(1));  % dV/dtheta1 = dV/dbeta * dbeta/dtheta1
    delE(:,2) = ( VEhat(:,2)/2./VEhat(:,3).*delF(:,2) )    *(+exp(theta(2)));   % dV/dtheta2 = dV/dgamma * dgamma/dtheta2
    delE(:,3) = ( VEhat(:,2)/2./VEhat(:,3).*delF(:,3)+1/b )*(-exp(theta(3)));   % dV/dtheta3 = dV/db * db/dtheta3
    delE(:,4) = ( VEhat(:,2)/2./VEhat(:,3).*delF(:,1) )    *normpdf(theta(4));  % dV/dtheta4 = dV/dc * dc/dtheta4
    
    dLdtheta = (dLdVE(:,1)*ones(1,p)).*delV + (dLdVE(:,2)*ones(1,p)).*delE;
    
    % information matrix part of the variance-covariance matrix
    OPG = cov(dLdtheta);  % using the fact that mean is zero, so OPG = cov
    
    % part of Hessian matrix contributed by V
    HessV = nan(T,p,p);  % time period; variance-covariance matrix of parameters(p x p)
    HessE = nan(T,p,p);  % time period; variance-covariance matrix of parameters(p x p)
    VCV = nan(p,p);      % variance-covariance matrix of parameters(p x p)
    for tt=1:T
        HessE(tt,:,:) = 1/(VEhat(tt,2)^2)*delE(tt,:)'*delE(tt,:);
    end
    
    kT = sum((abs(data(:)-VEhat(:,1))<cT));  % check this to make sure that the bandwidth is not "too small" -- we need at least a few returns that are close to the value at risk
    % Engle and Manganelli set this directly (see p373). They set kT=(40,60) for alpha=(.01,.05), when T=2892.
    for tt=1:T
        HessV(tt,:,:) = -1/alpha/(2*cT)*(abs(data(tt)-VEhat(tt,1))<cT)/VEhat(tt,2)*delV(tt,:)'*delV(tt,:);
    end
    Hess = squeeze(mean(HessE+HessV));
    temp = inv(Hess)*OPG*inv(Hess);
    VCV = (temp+temp')/2/T;  % To make sure the matrix is symmetric
    
    % now transforming the parameter vector back to the original units, and
    % using the delta method to get standard errors and t-stats
    thetaT = [normcdf(theta(1)) ; exp(theta(2)) ; -exp(theta(3))*normcdf(theta(4)) ; -exp(theta(3))  ];
    stderr = sqrt(diag(VCV));
    a34 = [-exp(theta(3))*normcdf(theta(4)) , -exp(theta(3))*normpdf(theta(4)) ];
    stderrT = sqrt([normpdf(theta(1))^2*(stderr(1)^2) ; (exp(theta(2))*stderr(2))^2 ;  a34*VCV(3:4,3:4)*a34' ; (-exp(theta(3))*stderr(3))^2 ]);
    outT = [thetaT,stderrT,thetaT./stderrT];
    
end

