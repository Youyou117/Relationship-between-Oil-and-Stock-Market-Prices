function out1 = factor_cop_Gcdf_GL_Skewtt(x, theta, GLweight)
% function out1 = factor_cop_Gcdf_GL_Skewtt(x, theta, GLweight)
%
% The marginal CDF of X_i associated with skew t - t factor model
%
% INPUTS:   x, a scalar, the value of x that we will evaluate G at
%           theta, =[lam, nuinv_z, nuinv_esp, psi_z], the parameters of Fz and Feps
%           GLweight, a matrix, nodes and weights for Gauss-Legendre quadrature
%
%  OUTPUTS: out1, a scalar, the value of the marginal CDF at x
%
%  Dong Hwan Oh and Andrew Patton
%
%  16 April 2015
%
%  This code is to accompany the paper:
%  Oh, D.H. and A.J. Patton, 2015, Time-Varying Systemic Risk: Evidence 
%  from a Dynamic Copula Model of CDS Spreads, working paper, Duke University. 

out1 = GLquad('factor_cop_Gcdf_calc1_Skewtt_loading', 1e-5, 1-1e-5, GLweight, x, theta);

end