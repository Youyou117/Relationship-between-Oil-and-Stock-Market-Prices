function out1 = nannewey_west(data,lag,pairwise);
%function out1 = nannewey_west(data,lag);
%
% Newey-West estimator of V[ n^(-1/2)*sum(data) ] 
% (equals, asymptotically, cov(data) if the data are uncorrelated)
%
%  ADJUSTED TO ALLOW FOR MISSING OBSERVATIONS (REPRESENTED BY NaN). NOTE
%  THAT THIS MATRIX MIGHT NOT BE POS DEF, SO CAVEAT EMPTOR!
%  (If no observations are missing, this function gives the same answer as
%  "newey_west.m", but is slightly slower due to the loops)
%
%
%  Andrew Patton
%
%  Tuesday 11 nov, 2003
%
% adapted: 26 Nov 2009

[T,K] = size(data);

if nargin<2 || isempty(lag)
    lag = floor(4*((T/100)^(2/9))); % this is the rule used by EViews
end
if nargin<3 || isempty(pairwise)
    pairwise=1;
end


if pairwise==0  % then only use rows where all variables are non-missing
    temp123 = ~isnan(sum(data')');  % indicator for rows where NONE of the variables are missing
    out1 = newey_west(data(temp123,:),lag);  % then just get the newey-West variance estimate for these rows
else  % then use all available data on each variable
    data = data - ones(T,1)*nanmean(data);
    
    B0 = data'*data/T;
    for ii=1:K;
        for jj=ii:K;
            B0(ii,jj) = nanmean(data(:,ii).*data(:,jj));
            B0(jj,ii) = B0(ii,jj);
        end
    end
    
    for kk=1:lag;
        for ii=1:K;
            for jj=1:K;
                B1(ii,jj) = nansum(data(1+kk:end,ii).*data(1:end-kk,jj))/T;
            end
        end
        B0 = B0 + (1-kk/(lag+1))*(B1+B1');
    end
    out1 = 1/2*(B0+B0');  % making sure this matrix is symmetric (the NaN's might affect this??)
end