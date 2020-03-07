function VEhat = skewt_VE(theta,alpha)
% function VEhat = skewt_VE(theta,alpha)
%
% Function to compute the theoretical VaR and ES for Hansen's skew t distribution
%
%  INPUTS:  theta, a Nx2 matrix containing [nu,lambda], the two parameters of Hansen's skew t distribution
%              (Note that setting lam=0 leads to the standardized (i.e., unit variance) Student's t distribution)
%           alpha, a px1 vector, the values of alpha to consider
%
%  OUTPUTS: VEhat, a Nx2xp matrix, the sample VaR and ES for each series
%
%  Andrew Patton
%
%  11 October 2018
%
% The formula below appeared in:
%
% Patton, A.J., J.F. Ziegel, and R. Chen, 2018, Dynamic Semiparametric
% Models for Expected Shortfall (and Value-at-Risk), Journal of Econometrics, forthcoming. 
%
% and is based on the following paper, which presented results for the Student's t distribution (among other results):
%
% Dobrev, D., T.D. Nesmith, D.H. Oh, 2017, Accurate Evaluation of Expected Shortfall for Linear Portfolios with 
% Elliptically Distributed Risk Factors, Journal of Risk and Financial Management, 10(5), 1-14.

if prod(size(theta))==2 && size(theta,1)==2 % ie, just a single value of (nu,lam) is entered but it is a column vector
    theta = theta';
end

    
N = size(theta,1);
p = length(alpha);
VEhat = nan(N,2,p);
for ii=1:N
    nu = theta(ii,1);
    lam = theta(ii,2);
    
    c = exp(gammaln((nu+1)/2)-gammaln(nu/2))/sqrt(pi*(nu-2));      % parameters used in Hansen's skew t PDF
    a = 4*lam*c*(nu-2)/(nu-1);
    b = sqrt(1+3*lam^2-a^2);

    if max(skewtdis_inv(alpha,nu,lam))>(-a/b)  % then at least one of the alpha values requires "special" treatment, so use a loop
        for aa=1:length(alpha)
            alpha1 = alpha(aa);
            VVZ = skewtdis_inv(alpha1,nu,lam);
            VEhat(ii,1,aa) = VVZ;  % theoretical VaR for the skew t
            if VVZ < (-a/b)  % then quantile is to the left of the mode of the distribution, and the calcs are a little easier)
                AAY = skewtdis_cdf( b/(1-lam)*(VVZ+a/b) ,nu,0);  % this is the "effective alpha" for the ES of the t variable
                esY = -sqrt((nu-2)/nu).*(nu.^(nu/2))./(2*AAY*sqrt(pi)).*exp(gammaln((nu-1)/2)-gammaln(nu/2)).*(( tinv(AAY,nu).^2 + nu).^(-(nu-1)/2));  % this is the analytical expression for ES of a studentised Student's t
                VEhat(ii,2,aa) = AAY./alpha1*(1-lam).*(-a/b + (1-lam)/b*esY);  % theoretical ES for the skew t
            else  % need to find ES to the *right* of the quantile, using the above simple formula, and then flip it around
                % the equations below use the fact that if Z~skew t(nu,lam) then -Z~skew t(nu,-lam)
                lam2 = -lam;    % this is part of the "flip" step
                nu2 = nu;       % unchanged
                c2 = c;
                a2 = 4*lam2*c2*(nu2-2)/(nu2-1);
                b2 = sqrt(1+3*lam2^2-a2^2);
                VVZ2 = skewtdis_inv(1-alpha1,nu2,lam2);

                AAY = skewtdis_cdf( b2/(1-lam2)*(VVZ2+a2/b2) ,nu2,0);  % this is the "effective alpha" for the ES of the t variable
                esY = -sqrt((nu2-2)/nu2).*(nu2.^(nu2/2))./(2*AAY*sqrt(pi)).*exp(gammaln((nu2-1)/2)-gammaln(nu2/2)).*(( tinv(AAY,nu2).^2 + nu2).^(-(nu2-1)/2));  % this is the analytical expression for ES of a studentised Student's t
                VEhat(ii,2,aa) = (1-alpha1)/alpha1 * AAY./(1-alpha1)*(1-lam2).*(-a2/b2 + (1-lam2)/b2*esY);  % theoretical ES for the skew t
            end
        end
    else  % no alpha value requires special treatment, and so I can vectorize
        VVZ = skewtdis_inv(alpha,nu,lam);
        VEhat(ii,1,:) = VVZ;  % theoretical VaR for the skew t
        AAY = skewtdis_cdf( b/(1-lam)*(VVZ+a/b) ,nu,0);  % this is the "effective alpha" for the ES of the t variable
        esY = -sqrt((nu-2)/nu).*(nu.^(nu/2))./(2*AAY*sqrt(pi)).*exp(gammaln((nu-1)/2)-gammaln(nu/2)).*(( tinv(AAY,nu).^2 + nu).^(-(nu-1)/2));  % this is the analytical expression for ES of a studentised Student's t
        VEhat(ii,2,:) = AAY./alpha*(1-lam).*(-a/b + (1-lam)/b*esY);  % theoretical ES for the skew t
    end
end
VEhat = squeeze(VEhat);  % removing any redundant dimensions
