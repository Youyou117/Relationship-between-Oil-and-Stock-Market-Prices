function [VEhat,vcv] = sample_VE(data,alpha)
% function [VEhat,vcv] = sample_VE(data,alpha)
%
% Function to compute the sample VaR and ES from a data set
%
%  INPUTS:  data, a TxN matrix, data on N time series
%           alpha, a px1 vector, the values of alpha to consider
%
%  OUTPUTS: VEhat, a Nx2xp matrix, the sample VaR and ES for each series
%           vcv, a 2x2xp matrix, the covariance matrix of the estimated VaR and ES
%
%  Andrew Patton
%
%  7 November 2015
%
% This code was used in: 
% Patton, A.J., J.F. Ziegel, and R. Chen, 2018, Dynamic Semiparametric
% Models for Expected Shortfall (and Value-at-Risk), Journal of Econometrics, forthcoming. 

[T,N] = size(data);
p = length(alpha);
VEhat = nan(N,2,p);
for ii=1:N
    data1 = data;
    for aa=1:p
        VEhat(ii,1,aa) = quantile(data,alpha(aa));
        if sum( (data1<=VEhat(ii,1,aa) ) )>0
            VEhat(ii,2,aa) = mean(data1(data1<=VEhat(ii,1,aa)));
        end
    end
    
    if nargout>1  % then return VCV matrix
        scores = [ -1/alpha/VEhat(ii,2,aa)*( (data1<=VEhat(ii,1,aa))-alpha ) , ...
            1/alpha/(VEhat(ii,2,aa)^2)*( (data1<=VEhat(ii,1,aa)).*(VEhat(ii,1,aa)-data1) + alpha*(VEhat(ii,2,aa)-VEhat(ii,1,aa)) ) ];
        OPG = cov(scores);  % using the fact that mean is zero, so OPG = cov
        
        cT = 1/(T^(1/3));  % bandwidth parameter for density estimation, needs to be smaller order than 1/sqrt(T)
        pdfv = 1/(2*T*cT)*sum( (abs(data1-VEhat(ii,1,aa))<cT)  );  % density of y at VaR
        
        hess11 =  -pdfv/alpha/VEhat(ii,2,aa);
        hess12 = 0;
        hess22 = 1/(VEhat(ii,2,aa)^2);
        
        hess = [[hess11,hess12];[hess12,hess22]];
        vcv = 1/T*inv(hess)*OPG*inv(hess);
    end
end
VEhat = squeeze(VEhat);  % just removing any redundant dimensions
if nargout>1
    vcv = squeeze(vcv);
end