function [out1 out2 ] = robust_loss_1(b,volproxy,hhat);
% function [out1 out2 ] = robust_loss_1(b,volproxy,hhat)
%
% This family of loss functions is useful for ranking volatility forecasts
% when the volatility proxy available (eg: daily squared returns) can be
% assumed to be conditionally unbiased for the latent object of interest
% (eg: conditional variance) but measures the object of interest with
% noise.
%
% INPUTS:       b, a Mx1 vector, the loss function shape parameter. (b=0 ->MSE, b=-2 ->QLIKE)
%               r2, a Tx1 vector of squared returns
%				hhat, a Txk vector of variance forecasts
%
% OUTPUTS:	out1, a TxkxM matrix of the loss function evaluated at each point in time
%               (out1 will be "squeezed" so all singleton dimensions are removed)
%           out2, a kxkxM matrix of t-statistics from Diebold-Mariano-West tests of equal predictive ability
%               (diagonals are left empty, element (i,j) is tstat on d(t) = L( volproxy(t),hhat_i(t) ) - L( volproxy(t),hhat_j(t) )
%
%
%
% This is the parametric family of loss functions that was proposed in:
%
% Patton, A.J., 2006, "Volatility Forecast Comparison using Imperfect Volatility Proxies", 
% manuscript, London School of Economics.
%
%  L(sig2,h;b) = 1/((b+1)*(b+2))*( sig2^(b+2) - hhat^(b+2) ) - 1/(b+1)*(hhat^(b+1))*(sig2-hhat),    for b not equal to -1 or -2
%  L(sig2,h;b) = hhat - sig2 + sig2*log(sig2/hhat),                                                 for b = -1
%  L(sig2,h;b) = sig2/hhat - log(sig2/hhat) - 1,                                                    for b = -2
%
% MSE is obtained for b=0, QLIKE is obtained for b=-2 (up to additive and multiplicative constants)
%
%
%  Andrew Patton
%
%  Thursday, 13 April, 2006.


[T,k] = size(hhat);
M = length(b);

out1 = nines(T,k,M);
out2 = nines(k,k,M);

if M==1  % trying to avoid using a loop
    if b==-1
        out1 = hhat - volproxy + volproxy.*log(volproxy./hhat);
    elseif b==-2;
        out1 = volproxy./hhat - log(volproxy./hhat) - 1;
    else
        out1 = 1/(b+1)/(b+2)*( (volproxy.^(b+2)) - (hhat.^(b+2)) ) - 1./(b+1)*(hhat.^(b+1)).*(volproxy-hhat);
    end
else  % M>1
    temp1 = find(b==-1);
    temp2 = find(b==-2);
    temp3 = find((b~=-1).*(b~=-2));
    if k==1;  % then there is only one series of forecasts to analyse.
        b = ones(T,1)*b';  % making size conformable with volproxy and hhat
        hhat = hhat*ones(1,M);
        volproxy = volproxy*ones(1,M);
        if size(volproxy,1)<T;
            volproxy = ones(T,1)*volproxy;
        end
        out1(:,1,temp1) = hhat(:,temp1) - volproxy(:,temp1) + volproxy(:,temp1).*log(volproxy(:,temp1)./hhat(:,temp1));
        out1(:,1,temp2) = volproxy(:,temp2)./hhat(:,temp2) - log(volproxy(:,temp2)./hhat(:,temp2)) - 1;
        out1(:,1,temp3) = 1./(b(:,temp3)+1)./(b(:,temp3)+2).*((volproxy(:,temp3).^(b(:,temp3)+2))-(hhat(:,temp3).^(b(:,temp3)+2))) - ...
            1./(b(:,temp3)+1).*(hhat(:,temp3).^(b(:,temp3)+1)).*(volproxy(:,temp3)-hhat(:,temp3));
    else  % in this case k>1 and M>1 (and I won't bother checking the special case that T=1 since that is not common) so use a loop
        volproxy = volproxy*ones(1,k);  % making volproxy conformable with hhat
        for ii=1:M;
            if b(ii)==-1
                out1(:,:,ii) = hhat - volproxy + volproxy.*log(volproxy./hhat);
            elseif b(ii)==-2;
                out1(:,:,ii) = volproxy./hhat - log(volproxy./hhat) - 1;
            else
                out1(:,:,ii) = 1/(b(ii)+1)/(b(ii)+2)*( (volproxy.^(b(ii)+2)) - (hhat.^(b(ii)+2)) ) - 1/(b(ii)+1)*(hhat.^(b(ii)+1)).*(volproxy-hhat);
            end
        end
    end
end
if nargout>1
    for kk=1:M;
        for ii=1:k-1;
            for jj=ii+1:k;
                dt = squeeze(out1(:,ii,kk) - out1(:,jj,kk));
                temp = nwest(dt,ones(T,1));
                out2(ii,jj,kk) = temp.tstat;
                out2(jj,ii,kk) = -temp.tstat;
            end
        end
    end
end
out1 = squeeze(out1);  % dropping all "singleton" dimensions
out2 = squeeze(out2);
