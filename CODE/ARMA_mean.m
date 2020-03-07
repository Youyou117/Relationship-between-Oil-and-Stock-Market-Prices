function muhat = ARMA_mean(theta,p,q,data);
%function muhat = ARMA_mean(theta,p,q,data);
%
% Function to return the conditional mean from an ARMA(p,q) model
%
%  INPUTS:  theta, a kx1 vector = [phi0,phi1,...,phip,theta1,...,thetaq]
%           p, a scalar, the order of the AR part
%           q, a scalar, the order of the MA part
%           data, a Tx1 vector of data
%
%  OUTPUS:  muhat, a Tx1 vector of the conditional mean
%
%  Andrew Patton
%
%  31 March 2011

T = length(data);

resids = zeros(T+max(p,q),1);
y = [mean(data)*ones(max(p,q),1);data];  % using sample mean of y for y(t) with t<1
muhat = nan(T+max(p,q),1);              
for tt=max(p,q)+1:max(p,q)+T;
    muhat(tt) = theta(1);
    for pp=1:p;
        muhat(tt) = muhat(tt) + theta(1+pp)*y(tt-pp);
    end
    for qq=1:q;
        muhat(tt) = muhat(tt) + theta(1+p+qq)*resids(tt-qq);
    end
    resids(tt) = y(tt)-muhat(tt);
end

muhat = muhat(max(p,q)+1:max(p,q)+T);  % dropping pre-sample obs