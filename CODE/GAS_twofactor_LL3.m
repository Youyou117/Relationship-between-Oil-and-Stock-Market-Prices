function [Eloss,VEhat,lamV,lamE,VCV,outT] = GAS_twofactor_LL3(theta,Y,alpha,VE0,tau,cT)
% function [Eloss,VEhat,lamV,lamE,VCV,outT] = GAS_twofactor_LL3(theta,Y,alpha,VE0,tau,cT)
%
% Function to return the average loss using the VaR_ES_loss function,
% assuming GAS dynamics for VaR and ES. 
%
%  Normalizing the score by the inverse expected Hessian. This matrix is
%  diagonal, so ends up just multiplying lamV by V*E*alpha, and then lamE  by E^2.
%
% INPUTS:   theta = [w1,w2,b1,b2,a11,a12,a21,a22], the parameters of the GAMM dynamics
%           Y, a Tx1 vector, returns on the target variable
%           alpha, a scalar in (0,1), the probability level
%           VE0, a 2x1 vector, starting value for VE time series. Default is unconditional VE from Y
%           tau, a scalar, a positive number if want to smooth the hit variable (eg, tau=10), use -1 (default) if want NO smoothing and want to use the true objective function
%           cT, a scalar, the bandwidth parameter to use in computing the VCV matrix. Default ccT=floor(T^(1/3))
%
% OUTPUTS:  Eloss, a scalar, the average loss using these parameter values
%           VEhat, a Tx2 matrix, the fitted values of VaR and ES at each point in time
%           VCV, a 8x8 matrix, the asymptotic covariance matrix of the estimated parameters (assuming that we are evaluating the function at a local min)
%           outT, a 8x3 matrix, the transformed parameters, std errors for the transformed parameters, and t-stats
%
%  Andrew Patton
%
%  28 Feb 2016
%
% This code was used in: 
% Patton, A.J., J.F. Ziegel, and R. Chen, 2018, Dynamic Semiparametric
% Models for Expected Shortfall (and Value-at-Risk), Journal of Econometrics, forthcoming. 


w1 = theta(1);
w2 = theta(2);
b1  = normcdf(theta(3));  % imposing this lies inside (0,1))
b2  = normcdf(theta(4));  % imposing this lies inside (0,1))
a11 = theta(5);
a12 = theta(6);
a21 = theta(7);
a22 = theta(8);

b = diag([b1,b2]);
a = [[a11,a12];[a21,a22]];
w = [w1;w2];

T = length(Y);
VEhat = nan(T,2);
if nargin<5 || isempty(tau)
    tau = -1;  % default option is to NOT use any smoothing of the indicator function
end

if nargin<4  || isempty(VE0)
    VE1(1,:) = VEhatTAU(Y,alpha,tau);
else
    VE1(1,:) = VE0';
end
VEhat(1,:) = VE1;


if nargin<6 || isempty(cT)
    cT = T^(-1/3);
end

if tau==-1
    hitS = (Y(1)<VEhat(1,1));  % then use the actual indicator function
else
    hitS = 1./(1+exp( tau*(Y(1)-VEhat(1,1)) ));
end

lamE = nan(T,1);
lamV = nan(T,1);
lamE(1)=0;
lamV(1)=0;
loss = zeros(T,1);
loss(1) = -1/alpha/VEhat(1,2)*hitS*(VEhat(1,1)-Y(1)) - 1/VEhat(1,2)*(VEhat(1,2)-VEhat(1,1)) + log(-VEhat(1,2));
for tt=2:T
    hitSL = hitS; 
    lamV(tt) = -VEhat(tt-1,1)*(hitSL - alpha);           
    lamE(tt) = 1/alpha*hitSL*Y(tt-1) - VEhat(tt-1,2);
    
    VEhat(tt,:) = ( w + b*VEhat(tt-1,:)' + a*[lamV(tt);lamE(tt)] )';
    if tau==-1
        hitS = (Y(tt)<VEhat(tt,1));  % then use the actual indicator function
    else
        hitS =  1./(1+exp( tau*(Y(tt)-VEhat(tt,1)) ));
    end
    loss(tt) = -1/alpha/VEhat(tt,2)*hitS*(VEhat(tt,1)-Y(tt)) - 1/VEhat(tt,2)*(VEhat(tt,2)-VEhat(tt,1)) + log(-VEhat(tt,2)) ;
end

Eloss = mean(loss);
if sum(VEhat(:,1)<VEhat(:,2))>0 || max(max(VEhat))>=0  % require ES<VaR, and since I'm using Hom deg 0 loss, also require both to be below 0
    Eloss = 1e6;
end
if sum(sum(isnan(VEhat)))>0 || ~isreal(VEhat) || isinf(Eloss) || ~isreal(theta) || sum(Y<=VEhat(:,1))==0  % adding rule that we must see *some* hits, otherwise parameter vector is not identified
    Eloss = 1e7;
end

if nargout>4  % then compute the VCV matrix
    p = length(theta);  % number of parameters in this model
    
    % first: the derivative of L w.r.t V and E. Same formula for all models for (V,E).
    dLdVE = nan(T,2);      % Note:in this model lamE (the RHS variable in the GAS model) is *not* the same as dL/dE, as we reparameterized the forcing variable. So compute it here.
    dLdVE (:,1) = -1/alpha./VEhat(:,2).*(Y<=VEhat(:,1) - alpha);  % dL/dV
    dLdVE (:,2) = +1/alpha./(VEhat(:,2).^2).*(Y<=VEhat(:,1)).*(VEhat(:,1)-Y) - VEhat(:,1)./(VEhat(:,2).^2) + 1./VEhat(:,2);   % dL/dE
    
    % second: deriv of V and E w.r.t theta, which depends on model
    delV = zeros(T,p);
    delE = zeros(T,p);
    for tt=2:T
        delV(tt,:) = [1,0,normpdf(theta(3))*VEhat(tt-1,1),0,  lamV(tt-1),lamE(tt-1),0,0] + (b1+a11*((Y(tt-1)<VEhat(tt-1,1))-alpha))*delV(tt-1,:) + -a12*delE(tt-1,:);
        delE(tt,:) = [0,1,0,normpdf(theta(4))*VEhat(tt-1,2),  0,0,lamV(tt-1),lamE(tt-1)] + a21*((Y(tt-1)<VEhat(tt-1,1))-alpha)*delV(tt-1,:) + (b2-a22)*delE(tt-1,:);
    end
    dLdtheta = (dLdVE(:,1)*ones(1,p)).*delV + (dLdVE(:,2)*ones(1,p)).*delE;
    
    % information matrix part of the variance-covariance matrix
    OPG = cov(dLdtheta);  % using the fact that mean is zero, so OPG = cov
    
    % part of Hessian matrix contributed by V
    HessV = nan(T,p,p);% time period; variance-covariance matrix of parameters(p x p)
    HessE = nan(T,p,p);% time period; variance-covariance matrix of parameters(p x p)
    VCV = nan(p,p);      % variance-covariance matrix of parameters(p x p)
    for tt=1:T
        HessE(tt,:,:) = 1/(VEhat(tt,2)^2)*delE(tt,:)'*delE(tt,:);
    end
    for tt=1:T
        HessV(tt,:,:) = -1/alpha/(2*cT)*(abs(Y(tt)-VEhat(tt,1))<cT)/VEhat(tt,2)*delV(tt,:)'*delV(tt,:);
    end
    Hess = squeeze(mean(HessE+HessV));
    temp = inv(Hess)*OPG*inv(Hess);
    VCV = (temp+temp')/2/T;  % To make sure the matrix is symmetric
    
    thetaT = [w1;w2;b1;b2;a11;a12;a21;a22];
    stderr = sqrt(diag(VCV));
    stderrT = [ stderr(1:2);   sqrt(normpdf(theta(3))^2*(stderr(3)^2)) ; sqrt(normpdf(theta(4))^2*(stderr(4)^2)) ; stderr(5:8) ];
    outT = [thetaT,stderrT,thetaT./stderrT];
    
end
