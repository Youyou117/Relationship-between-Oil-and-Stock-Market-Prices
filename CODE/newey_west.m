function out1 = newey_west(data,lag);
%function out1 = newey_west(data,lag);
%
% Newey-West estimator of V[ n^(-1/2)*sum(data) ] 
% (equals, asymptotically, cov(data) if the data are uncorrelated)
%
%  Andrew Patton
%
%  Tuesday 11 nov, 2003

[T,K] = size(data);

if nargin<2
    lag = floor(4*((T/100)^(2/9))); % this is the rule used by EViews
end

data = data - ones(T,1)*mean(data);

B0 = data'*data/T;
for ii=1:lag;
    B1 = data(1+ii:end,:)'*data(1:end-ii,:)/T;
    B0 = B0 + (1-ii/(lag+1))*(B1+B1');
end
out1 = B0;
