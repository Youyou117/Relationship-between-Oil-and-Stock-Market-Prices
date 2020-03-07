function [sigmahat1,sigmahat2,out1,scoresS,Qjt,Qhatat] = copula_deriv2_chen_fan(function_name,theta,data,varargin)
% function out1 = copula_deriv2_chen_fan(function_name,theta,varargin)
%
% This program computes d2(log copula likelihood)/( dUj * dtheta ), which
% is used in getting standard errors for the class of semiparametric copula
% models considered by Chen and Fan (2006, JoE)
%
% INPUTS:   theta, a kx1 vector of parameters
%           function_name, a string containing the name of a function that returns the value of the log-likelihood *at each point in time*
%           data, a Txd matrix, the matrix of Unif(0,1) variables that are distributed according to this copula
%           varargin, any other arguments to pass to function_name
%
% OUTPUTS:  sigmahat1, the central part of the asymptotic covariance matrix of Chen and Fan (2006) (ie, the simghat hat in  V=inv(B)*sigmahat*inv(B)). Imposing that scores are mean zero.
%           sigmahat2, the central part of the asymptotic covariance matrix of Chen and Fan (2006) (ie, the simghat hat in  V=inv(B)*sigmahat*inv(B)). NOT imposing that scores are mean zero
%           out1, a Txkxd matrix of second derivatives at each point in time
%
%  1 April 2011
%
%  Andrew Patton

k = size(theta,1);
[T,d] = size(data);

out1 = nan(T,k,d);


% Compute the stepsize (h)
hU = 0.0001;  % making step size for U smaller, so that i don't go outside [0,1]
hA = 0.0001*ones(k,1);
eU = hU*eye(d);
eA = diag(hA);

% hA = eps.^(1/3)*max(abs(theta),1e-4)
% hA = 1e-4*abs(theta)
% hU = 0.0001;
% eU = hU*eye(d);
% eA = diag(hA);

fx = feval(function_name,theta,data,varargin{:});

% computing the two first derivative terms
dfx1 = nan(T,d);
dfdU = nan(T,d);
for jj=1:d;
    dataj = data;
    dataj(:,jj) = data(:,jj) + hU;  
    dataj(:,jj) = max(min(dataj(:,jj),1-1e-6),1e-6);  % making sure adjusted data is insed [1e-6 , 1-1e-6]
    dfx1(:,jj) = feval(function_name,theta,dataj,varargin{:});  % log like at each point in time, when Uj is slightly larger
    dfdU(:,jj) = (dfx1(:,jj)-fx)/hU;  % first derivative of the log-like w.r.t. U[j,t]
end

dfx2 = nan(T,k);
scores = nan(T,k);
for ii=1:k;
    dfx2(:,ii) = feval(function_name,theta+eA(:,ii),data,varargin{:});  % log like at each point in time, when thetai is slightly larger
    scores(:,ii) = (dfx2(:,ii)-fx)/hA(ii);
end


% computing the second deriv term
d2fx = nan(T,k,d);
for jj=1:d;
    dataj = data;
    dataj(:,jj) = data(:,jj) + hU;  
    dataj(:,jj) = max(min(dataj(:,jj),1-1e-6),1e-6);  % making sure adjusted data is insed [1e-6 , 1-1e-6]
    for ii=1:k;
        d2fx(:,ii,jj) = feval(function_name,theta+eA(:,ii),dataj,varargin{:});  % log like at each point in time, when Uj is slightly larger AND theati is slightly larger
    end
end

% finally creating the second derivative
for jj=1:d;
    for ii=1:k;
        out1(:,ii,jj) = ( d2fx(:,ii,jj) - dfx1(:,jj) - dfx2(:,ii) + fx ) / (hU*hA(ii));
    end
end
   
Qhatat = nan(T,k,d);  % used in Chen-Fan std error calculations
Qjt = nan(T,d);       % used in Chen-Fan PLR test  
for tt=1:T;
    temp124 = setdiff((1:T)',tt);  % all dates except for this one
    for jj=1:d;
        for ii=1:k;
            Qhatat(tt,ii,jj) = sum( d2fx(temp124,ii,jj) .* ( (data(tt,jj)<=data(temp124,jj)) - data(temp124,jj) ) )/T;
        end
        Qjt(tt,jj) = sum( dfdU(temp124,jj) .* ( (data(tt,jj)<=data(temp124,jj)) - data(temp124,jj) ) )/T;  % eq 4.3 on p137
    end
end

scoresS = [scores + sum(Qhatat,3)];
sigmahat2 = scoresS'*scoresS/T;  % not imposing that the scores are mean zero (which holds in the limit)
sigmahat1 = cov(scoresS);        % imposing that scores are mean zero

