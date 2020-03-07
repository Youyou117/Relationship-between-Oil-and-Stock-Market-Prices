function hhat = GJR_GARCH_vol(theta,p,o,q,data);
%function hhat = GJR_GARCH_vol(theta,p,o,q,data);
%
% Function to return the conditional variance from a GJR-GARCH(p,o,q) model
%
%  INPUTS:  theta, a kx1 vector = [omega,arch1,...,archp,tarch1,...,tarcho,garch1,...,garchq]
%           p, a scalar, the order of the ARCH part (coeff on eps[t-j]^2
%           o, a scalar, the order of the TARCH part (coeff on (eps[t-j]^2)*1(eps[t-j]<0)
%           q, a scalar, the order of the MA part (coeff on sig2[t-j]
%           data, a Tx1 vector of data (assumed to be mean zero, as is std in GARCH models)
%
%  OUTPUS:  hhat, a Tx1 vector of the conditional mean
%
%  Andrew Patton
%
%  31 March 2011

T = length(data);

eps = [zeros(max([p,o,q]),1);data];  % using sample mean of y for y(t) with t<1
hhat = [cov(data)*ones(max([p,o,q]),1);nan(T,1)];
for tt=max([p,o,q])+1:max([p,o,q])+T;
    hhat(tt) = theta(1);
    for pp=1:p;
        hhat(tt) = hhat(tt) + theta(1+pp)*(eps(tt-pp)^2);
    end
    for oo=1:o;
        hhat(tt) = hhat(tt) + theta(1+p+oo)*( (eps(tt-oo)^2).*(eps(tt-oo)<0) );
    end
    for qq=1:q;
        hhat(tt) = hhat(tt) + theta(1+p+o+qq)*hhat(tt-qq);
    end
end

hhat = hhat(max([p,o,q])+1:max([p,o,q])+T);  % dropping pre-sample obs