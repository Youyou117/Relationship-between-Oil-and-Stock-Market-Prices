function [out1,out2] = VEforecast_filter(VEforecast,mode,VE0)
% function out = VEforecast_filter(VEforecast,mode,VE0);
%
% Function to take in a time series of VaR and ES forecasts, and make sure
% they are consistent with each other.
%
%  INPUTS:  VEforecast, either at Tx2 matrix, or a TxKx2 matrix, with VaR forecasts in first column and ES forecasts in second
%           mode, =1 (default) if want to use previous filtered forecast when current is "insane", =0  if want to use uncond forecast
%           VE0, either a 1x2 vector, or a kx2 matrix, the value to use if mode=0
%
%  OUTPUTS: out1, a Tx2 matrix or a TxKx2 matrix, filtered VaR and ES forecasts
%           out2, a TxK matrix flagging the time periods and models where a filter was used
%
%  Andrew Patton
% 
% 5 Dec 2015
%
% This code was used in: 
% Patton, A.J., J.F. Ziegel, and R. Chen, 2018, Dynamic Semiparametric
% Models for Expected Shortfall (and Value-at-Risk), Journal of Econometrics, forthcoming. 



T = size(VEforecast,1);
if size(VEforecast,3)==1
    K=1;
    temp = nan(T,1,2);
    temp(:,1,:) = VEforecast;
    VEforecast = temp;
else
    K = size(VEforecast,2);
end

if nargin<2 || isempty(mode)
    mode = 1;
end

out1 = nan(T,K,2); 
out2 = nan(T,K);
for kk=1:K
    if VEforecast(1,kk,1)<VEforecast(1,kk,2) || max(VEforecast(1,kk,:))>0   % then VaR<ES which makes no sense, or one of these are positive, which is a problem for the FZ0 loss function
        if ~isempty(VE0)  % for t=1, use the uncond for both "modes", but only if it is provided.
            out1(1,kk,:) = VE0(kk,:);
        end
    else
        out1(1,kk,:) = VEforecast(1,kk,:);
    end
    for tt=2:T
        if VEforecast(1,kk,1)<VEforecast(1,kk,2)  || max(VEforecast(tt,kk,:))>0    % then VaR<ES which makes no sense, or one of these are positive, which is a problem for the FZ0 loss function
            if mode==0
                out1(tt,kk,:) = VE0(kk,:);
            elseif mode==-1
                out1(tt,kk,:) = out1(tt-1,kk,:);
            end
            out2(tt,kk) = 1;
        else
            out1(tt,kk,:) = VEforecast(tt,kk,:);
            out2(tt,kk) = 0;
        end
    end
end

out1 = squeeze(out1);  % squeezing out the second dimension if it equals 1