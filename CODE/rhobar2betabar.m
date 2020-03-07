function betabar = rhobar2betabar(rhobar)
% function betabar = rhobar2betabar(rhobar);
%
% Function to take in a sample rank correlation matrix and return a vector of implied loadings on the common factor
% (assuming a simple one-factor structure for the data).
%
% NOTE: This function maps the sample rank correlation to the model-implied *linear* correlation (which has a nice closed-form expression). The difference
% between rank and linear correlation is small for rho inside about [-0.75,0.75]. For rho above 0.8 or so we find rank correlation is below linear
% correlation, which leads to estimated betabar being below the true beta. So be wary applying this function to highly correlated series.
%
%  INPUTS:  rhobar, a NxN matrix, the sample rank correlation matrix
%
%  OUTPUTS: betabar, a Nx1 vector, the estimated implied loadings on the common factor
%
%  Dong Hwan Oh and Andrew Patton
%
%  17 October 2012
%
%  This code is to accompany the paper:
%  Oh, D.H. and A.J. Patton, 2015, Time-Varying Systemic Risk: Evidence 
%  from a Dynamic Copula Model of CDS Spreads, working paper, Duke University. 



N = size(rhobar,1);
if N<3
    'this mapping requires there to be at least 3 variables'
end

theta0 = ones(N,1);
options = optimset('Display','off','TolCon',10^-12,'TolFun',10^-4,'TolX',10^-6);
rhobar2betabar_calc(theta0,rhobar);
warning off;
betabar = fminunc(@(beta) rhobar2betabar_calc(beta,rhobar),theta0,options);
warning on;


    function out1 = rho2theta(rho)
        % this function takes in a NxN correlation matrix and returns a vector of just the upper-triangle of correlations
        k = size(rho,2);
        out1 = nines(k*(k-1)/2,1);
        counter=1;
        for ii=1:k;
            for jj=ii+1:k
                out1(counter) = rho(ii,jj);
                counter=counter+1;
            end
        end
    end


    function LL = rhobar2betabar_calc(beta,rhobar)
        % objective function for estimating betabar
        Nb = length(beta);
        rho = nan(Nb,Nb);
        for ii=1:Nb;
            for jj=ii+1:Nb;
                rho(ii,jj) = beta(ii)*beta(jj)/sqrt( (1+beta(ii)^2)*(1+beta(jj)^2) );
                rho(jj,ii) = rho(ii,jj);
            end
        end
        LL = sum(rho2theta( (rho-rhobar).^2) );  % summing the squared diff between the model implied rho and the sample rhobar
    end


end
