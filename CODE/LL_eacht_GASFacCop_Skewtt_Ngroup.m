function [LL_t , N_derivative] = LL_eacht_GASFacCop_Skewtt_Ngroup(theta_t, u_t, GLweight, Gcdf, Gpdf, x_grid, lam_grid, group_code, epsi)
% [LL_t , N_derivative] = LL_eacht_GASFacCop_Skewtt_Ngroup(theta_t, u_t, GLweight, Gcdf, Gpdf, x_grid, lam_grid, group_code, epsi)
%
% This function computes the log likelihood of (skew t - t) factor copula
% evaluated at time t. It also generates the numerical derivatives of log copula
% density with respect to each group's factor loading at time t
%
% INPUTS: theta,    a vector of parameters, [factor loadings; nuinv_z; nuinv_eps; psi_z]
%         u_t,      a Nx1 vector Unif(0,1) that are distributed according to this copula at time t
%         GLweight, a matrix, nodes and weight for Gauss-Legendre quadrature
%         Gcdf,     a [Num_x_grid by Num_lam_grid] matrix of the marginal cdfs of skew t-t factor model at x and factor loading (lam)
%         Gpdf,     a [Num_x_grid by Num_lam_grid] matrix of the marginal pdfs of skew t-t factor model at x and factor loading (lam)
%         x_grid,   a vector of x that Gcdf and Gpdf are evaluated at
%         lam_grid, a vector of factor loading that Gcdf and Gpdf are evaluated at
%         group_code, a Nx1 vector of group codes into which each firm is classified
%         epsi,     a scalar, the step size for numerical derivative for score calculation
%
% OUTPUTS:
%          LL_t, a scalar, the log-likelihood of (skew t - t) factor copula evaluated at each time t
%          N_derivative, a (N_group x 1) vector, each element is the numerical derivative of log copula density
%                        with respect to each group's factor loading
%
%  Dong Hwan Oh and Andrew Patton
%
%  16 April 2015
%
%  This code is to accompany the paper:
%  Oh, D.H. and A.J. Patton, 2015, Time-Varying Systemic Risk: Evidence
%  from a Dynamic Copula Model of CDS Spreads, working paper, Duke University.


if size(theta_t,2)>size(theta_t,1)
    theta_t = theta_t';
end

if size(u_t,2)>size(u_t,1)
    u_t = u_t';
end

if length(group_code) ~= length(u_t)
    error('The length of group_code does not equal the length of u');
end

if length(theta_t) ~= (3+max(group_code))
    error('theta or group_code has incorrect size');
end

Ginv_u_t        = nan(size(u_t,1),1);
Ginv_u_t_ateps  = nan(size(u_t,1),1);
deno_t          = nan(size(u_t,1),1);
deno_t_ateps    = nan(size(u_t,1),1);

for i = 1: max(group_code)
    inx = find((group_code == i));
    u = u_t(inx,1);
    
    theta_group       = [theta_t(i)        ; theta_t(end-2); theta_t(end-1); theta_t(end)];
    theta_group_eps   = [theta_t(i) + epsi ; theta_t(end-2); theta_t(end-1); theta_t(end)];
    
    Ginv_u_t(inx,1)       = factor_cop_Gcdfinv_spline(u, theta_group,     Gcdf, x_grid, lam_grid);
    Ginv_u_t_ateps(inx,1) = factor_cop_Gcdfinv_spline(u, theta_group_eps, Gcdf, x_grid, lam_grid);
    
    deno_t(inx,1)         = factor_cop_Gpdfvec_Skewtt(Ginv_u_t(inx,1),       theta_group,     Gpdf, x_grid, lam_grid);
    deno_t_ateps(inx,1)   = factor_cop_Gpdfvec_Skewtt(Ginv_u_t_ateps(inx,1), theta_group_eps, Gpdf, x_grid, lam_grid);
    
end

numer_t = factor_cop_FXpdf_GL_Skewtt([Ginv_u_t, Ginv_u_t_ateps], theta_t, GLweight, group_code, epsi);

LL_t = log(numer_t(1))- sum(log(deno_t));

%%% The numerical derivative of log copula density w.r.t. each group's factor loading
N_derivative = nan(max(group_code),1);
if nargout>=2
    for i = 1: max(group_code)
        
        inx = find((group_code == i));
        
        deno_temp = deno_t;
        deno_temp(inx,1) = deno_t_ateps(inx,1);
        
        LL_eps = log(numer_t(1+i)) - sum(log(deno_temp));
        
        N_derivative(i) = (LL_eps - LL_t)/epsi ;
    end
end

end
