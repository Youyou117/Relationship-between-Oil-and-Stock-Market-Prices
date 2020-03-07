function [out, LL, lam, log_lam, s ] = LL_GASFacCop_Skewtt_VT(theta, data_u, GLweight, lam_bar, lam_ini)
% [out, LL, lam, log_lam, s ]  = LL_GASFacCop_Skewtt_NGroup(theta, data_u , GLweight, group_code, lam_ini)
%
% This function computes the sum of (negative) log likelihoods of factor copula (skew t - t) with
% GAS recursion
%
% INPUTS:  theta,        a vector of parameters, [alpha, beta, nuinv_z, nuinv_eps, psi_z]
%          data_u,       a TxN matrix, the matrix of Unif(0,1) to be modeled with this copula
%          GLweight,     a [vec of nodes x vec of weights] matrix,  nodes and weight for Gauss-Legendre quadrature
%          lam_bar,      a (N x 1) vector of (average) factor loadings given from the separate stage     
%          lam_ini,      a (N x 1) vector of factor loadings at t=1
%
% OUTPUTS: out,     a scalar, the sum of (negative) log-likelihoods of factor copula evaluated at each time t
%          LL,      a Tx1 vector of log-likelihoods of factor copula evaluated at each time t
%          lam,     a (T x N) matrix of (time-varying) factor loadings
%          log_lam, a (T x N) matrix of (time-varying) log of factor loadings
%          s,       a (T x N) matrix of (time-varying) score
%
%  Dong Hwan Oh and Andrew Patton
%
%  16 February 2016
%
%  This code is to accompany the paper:
%  Oh, D.H. and A.J. Patton, 2015, Time-Varying Systemic Risk: Evidence
%  from a Dynamic Copula Model of CDS Spreads, working paper, Duke University.

if size(theta,2)>size(theta,1)
    theta = theta';
end

[TT,NN] = size(data_u);

group_code = (1:1:NN);
Ngroup = max(group_code);

epsi = 0.001 ; % step size of the numerical derivative for score

alpha = theta(1);
beta  = theta(2);
nuinv_z = theta(3);
nuinv_eps = theta(4);
psi_z   = theta(5);

omega = (1-beta) * log(lam_bar) ;

LL = nan(TT,1);

lam_ini(lam_ini>4) = 4;
lam_ini(lam_ini<0.01) = 0.01;

lam = nan(TT,Ngroup);
lam(1,:) = lam_ini';

log_lam = nan(TT,Ngroup);
log_lam(1,:) = log(lam_ini);

s  = nan(TT,Ngroup);

%%% Evaluate the marginal cdf (G) and pdf (g) of skew t - t factor model at various x and factor loadings fixing other
% parameters such as nuinv_z, nuinv_eps, and psi_z.
x1 = -15;
x2 = 15;
Npoints = 100;
x_grid = (x1:(x2-x1)/(Npoints-1):x2)';
x_grid = [-30;x_grid;30];

lam_grid = [0.001,0.01,(0.05:0.05:2.5)]' ;

[Gcdf_ini, Gpdf_ini] = Gcdfpdf_Skewtt([nuinv_z; nuinv_eps; psi_z], GLweight, x_grid, lam_grid);

%%% Interpolate the marginal cdf and pdf along with finer grids of factor loadings (lam)
dense_lam_grid = [0.001;(0.01:0.001:2.5)'];
Gcdf = nan(length(x_grid),length(dense_lam_grid));
Gpdf = nan(length(x_grid),length(dense_lam_grid));

for i = 1:length(dense_lam_grid)
    Gcdf(:,i) = interp2(lam_grid, x_grid, Gcdf_ini, dense_lam_grid(i), x_grid);
    Gpdf(:,i) = interp2(lam_grid, x_grid, Gpdf_ini, dense_lam_grid(i), x_grid);
end

%%%
%%% Evaluate log density of skew t - t factor copula with GAS recursion at each time t
for tt = 1:TT
    
    if tt ~= 1
        log_lam(tt,:) = omega'+ alpha*s(tt-1,:)+ beta*log_lam(tt-1,:);
        lam(tt,:) = exp(log_lam(tt,:));
    end
    
    lam(tt, (lam(tt,:) < 0.01)) = 0.01;
    lam(tt, (lam(tt,:) > 2.5 )) = 2.5;
    
    [L_temp, N_derivative] = LL_eacht_GASFacCop_Skewtt_Ngroup([lam(tt,:)'; nuinv_z; nuinv_eps; psi_z], data_u(tt,:)', GLweight, Gcdf, Gpdf, x_grid, dense_lam_grid, group_code, epsi) ;
    
    LL(tt) = L_temp;
    s(tt,:) = N_derivative'.*lam(tt,:);   
    
end

out = -sum(LL);


