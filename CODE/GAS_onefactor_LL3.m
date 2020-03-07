function [Eloss,VEhat,loss,VCV,outT] = GAS_onefactor_LL3(theta,Y,alpha,tau,cT)
% function [Eloss,VEhat,loss,VCV,outT] = GAS_onefactor_LL3(theta,Y,alpha,tau,cT)
%
% Function to return the average loss for the one-factor GAS model on the FZ loss function.
%
% Using the hom deg 0 FZ loss and the inverse Hessian as scaling matrix, we get the evolution equation as follows:
%
%  y(t) = exp{f(t)}*eta(t),  eta(t)~iid F(0,1)
%  f(t) = omega + beta*f(t) + gamma/(-e(tt-1))*( 1/alpha*hit(tt-1)*y(tt-1) - e(tt-1) )
% 
%  and v(t) = a*exp{f(t)}, e(t) = b*exp{f(t)}
%
%  INPUTS:  theta = [beta,gamma,b,c], We set beta=normcdf(beta) and gamma=exp(gamma), b=-exp(b) and c=normcdf(c). We set omega=0 as this parameter is not identified.
%           Y, a Tx1 vector of returns
%           alpha, a scalar inside (0,1), the probability level
%           tau, a scalar, a positive number if want to smooth the hit variable (eg, tau=10), use -1 (default) if want NO smoothing and want to use the true objective function
%                   tau=5 is quite a bit of smoothing, tau=10 is some, tau=20 is not much. tau=50 is almost none.
%           cT, a scalar, the bandwidth parameter to use in computing the VCV matrix. Default ccT=floor(T^(1/3))
%
%  OUTPUTS: Eloss, a scalar, the average FZ loss
%           VEhat, a Tx3 vector, the model-implied [VaR, ES, factor] series
%           loss, a Tx1 vector, the FZ loss values at each point in time
%           VCV, a 4x4 matrix, the asymptotic covariance matrix of the estimated parameters (assuming that we are evaluating the function at a local min)
%           outT, a 4x3 matrix, the transformed parameters, std errors for the transformed parameters, and t-stats
%
%  Andrew Patton
%
%  6 Dec 2016
%
% This code was used in: 
% Patton, A.J., J.F. Ziegel, and R. Chen, 2018, Dynamic Semiparametric
% Models for Expected Shortfall (and Value-at-Risk), Journal of Econometrics, forthcoming. 


beta = normcdf(theta(1));
gamma = exp(theta(2));
b = -exp(theta(3));
c = normcdf(theta(4));
a = c*b;

omega = 0;  % use the true value of omega here, if this is for a simulation, so that [b,c] are centered on the right values. If wrong omega is used the fit will be identical, but we cannot compare bias and efficiency so easily.

kalpha = 1;  % this is just a normalizing constant. Can set it to 1, or 1/sqrt(alpha), or any other positive value. May want to choose a smaller value than 1 when parameter is close to zero (numerical optimization issues).

T = length(Y);

if nargin<4 || isempty(tau)
    tau = -1;  % default option is to NOT use any smoothing of the indicator function
end
if nargin<5 || isempty(cT)
    cT = T^(-1/3);
end


fhat = nan(T,1);  % this is the "one factor" driving the dynamics. It is labelled "kappa" in the paper.
VEhat = nan(T,2);
fhat(1) = omega/(1-beta);  % model-implied unconditional average for f
VEhat(1,:) = [a,b]*exp(fhat(1));

if tau==-1
    hitS = (Y(1)<VEhat(1,1));  % then use the actual indicator function
else
    hitS = 1./(1+exp( tau*(Y(1)-VEhat(1,1)) ));
end

loss = zeros(T,1);
loss(1) = -1/alpha/VEhat(1,2)*hitS*(VEhat(1,1)-Y(1)) - 1/VEhat(1,2)*(VEhat(1,2)-VEhat(1,1)) + log(-VEhat(1,2));
for tt=2:T
    hitSL = hitS;  % lagged hit
    fhat(tt) = omega + beta*fhat(tt-1) + gamma/kalpha/VEhat(tt-1,2)*( 1/alpha*hitSL*Y(tt-1) - VEhat(tt-1,2) );
    VEhat(tt,:) = [a,b]*exp(fhat(tt));
        
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
VEhat(:,3) = fhat;  % put this here so it doesn't mess up the above set of insanity checks

if nargout>3  % then compute the VCV matrix
    p = length(theta);  % number of parameters in this model
    
    % first: the derivative of L w.r.t V and E. Same formula for all models for (V,E).
    dLdVE = nan(T,2);      % Note:in this model lamE (the RHS variable in the GAS model) is *not* the same as dL/dE, as we reparameterized the forcing variable. So compute it here.
    dLdVE (:,1) = -1/alpha./VEhat(:,2).*(Y<=VEhat(:,1) - alpha);  % dL/dV
    dLdVE (:,2) = +1/alpha./(VEhat(:,2).^2).*(Y<=VEhat(:,1)).*(VEhat(:,1)-Y) - VEhat(:,1)./(VEhat(:,2).^2) + 1./VEhat(:,2);   % dL/dE
    
    % second: deriv of V and E w.r.t theta, which depends on model
    delV = zeros(T,p);
    delE = zeros(T,p);
    delF = zeros(T,p);  % delV and delE are both driven by delF, so will create this here too. (f=kappa=log(sigma))
    delF(1,:) = [normpdf(theta(1))/(1-normcdf(theta(1))),exp(theta(2))/kalpha,0,0];  % starting value for delF
    delE(1,:) = VEhat(1,2)*(delF(1,:)+[0,0,1,0]);
    for tt=2:T
        psiL = -exp(theta(2))/kalpha/alpha*(  Y(tt-1)<VEhat(tt-1,1))*Y(tt-1)/(VEhat(tt-1,2)^2  );
        psi2 = +exp(theta(2))/kalpha/alpha*( (Y(tt-1)<VEhat(tt-1,1))*Y(tt-1)/VEhat(tt-1,2) - 1 );
        delF(tt,:) = normcdf(theta(1))*delF(tt-1,:) + psiL*delE(tt-1,:) + [normpdf(theta(1))*VEhat(tt-1,3) , psi2, 0, 0];
        delE(tt,:) = VEhat(tt,2)*(delF(tt,:) + [0,0,1,0]);
    end
    delV = (VEhat(:,1)*ones(1,4)).*(delF + ones(T,1)*[0,0,1,normpdf(theta(4))/normcdf(theta(4))]);
    
    dLdtheta = (dLdVE(:,1)*ones(1,p)).*delV + (dLdVE(:,2)*ones(1,p)).*delE;
    
    % information matrix part of the variance-covariance matrix
    OPG = cov(dLdtheta);  % using the fact that mean is zero, so OPG = cov
    
    % part of Hessian matrix contributed by V
    HessV = nan(T,p,p);  % time period; V part of Hessian matrix of parameters(p x p)
    HessE = nan(T,p,p);  % time period; E part of Hessian matrix of parameters(p x p)
    for tt=1:T
        HessE(tt,:,:) = 1/(VEhat(tt,2)^2)*delE(tt,:)'*delE(tt,:);
    end
    
    kT = sum((abs(Y-VEhat(:,1))<cT))  % check this to make sure that the bandwidth is not "too small" -- we need at least a few returns that are close to the value at risk
    % Engle and Manganelli set this directly (see p373). They set kT=(40,60) for alpha=(.01,.05), when T=2892.
    for tt=1:T
        HessV(tt,:,:) = -1/alpha/(2*cT)*(abs(Y(tt)-VEhat(tt,1))<cT)/VEhat(tt,2)*delV(tt,:)'*delV(tt,:);
    end
    Hess = squeeze(mean(HessE+HessV));
    temp = inv(Hess)*OPG*inv(Hess);
    VCV = (temp+temp')/2/T;  % To make sure the matrix is symmetric (solves numerical issues sometimes)
    
    thetaT = [normcdf(theta(1)) ; exp(theta(2)) ; -exp(theta(3))*normcdf(theta(4)) ; -exp(theta(3))  ];
    stderr = sqrt(diag(VCV));
    a34 = [-exp(theta(3))*normcdf(theta(4)) , -exp(theta(3))*normpdf(theta(4)) ];
    stderrT = sqrt([normpdf(theta(1))^2*(stderr(1)^2) ; (exp(theta(2))*stderr(2))^2 ;  a34*VCV(3:4,3:4)*a34' ; (-exp(theta(3))*stderr(3))^2 ]);
    outT = [thetaT,stderrT,thetaT./stderrT];
    
end
