function out1 = SB_variance(data,b);
% function out1 = SB_variance(data,b);
%
%  Function to compute the stationary bootstrap estimate of 
%  V[ sqrt(T)*Xbar ] (ie: the long-run variance) ANALYTICALLY
%
%  Using Lemma 1 of Politis and Romano (1994)
%
% INPUTS:   data, a TxN matrix of data
%           b, a scalar, the average block length in the stationary bootstrap
%
% OUTPUTS,  out1, a NxN matrix, the estimated long-run covariance matrix of the data
%
% 
%  Andrew Patton
%
%  4 May, 2008

[T,N] = size(data);

data = data - ones(T,1)*mean(data);  % de-meaning the data to start off with. now I can drop the mean parts below

data = [data;data];  % stacking the data on top of itself - PR's Lemma 1 uses "circular" autocovariances and this is useful below.

out1 = data(1:T,:)'*data(1:T,:)/T;
%for ii=1:N-1;
for ii=1:min(N,T)-1;  % dealing with case that N>T
        out1 = out1 + (1-ii/N)*((1-1/b)^ii) * ...
            ( data(1:T,:)'*data(1+ii:T+ii,:)/T  + data(1+ii:T+ii,:)'*data(1:T,:)/T  )    ;  % 4apr09: bug corrected
end
