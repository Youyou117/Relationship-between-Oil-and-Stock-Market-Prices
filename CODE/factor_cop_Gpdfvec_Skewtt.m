function out = factor_cop_Gpdfvec_Skewtt(q, theta, Gpdf, x_grid, lam_grid)
% function out = factor_cop_Gpdfvec_Skewtt(q, theta, Gpdf, x_grid, lam_grid)
%
% The vector of (approximate) marginal densities at the vector q associated with skew t - t factor model 
%
% INPUTS: q,        a vector, the values that we will (approximately) evaluate marginal pdfs at
%         theta,  = [factor loading, nuinv_z, nuinv_esp, psi_z], the parameters of Fz and Feps
%         Gpdf,     a matrix of the marginal pdfs of skew t-t factor model at x and factor loading (lam)
%         x_grid,   a vector of x that Gpdf is evaluated at
%         lam_grid, a vector of factor loading that Gpdf is evaluated at
%
%  OUTPUTS: out, a vector, the marginal pdfs approximated at q
%
%  Dong Hwan Oh and Andrew Patton
%
%  16 April 2015
%
%  This code is to accompany the paper:
%  Oh, D.H. and A.J. Patton, 2015, Time-Varying Systemic Risk: Evidence 
%  from a Dynamic Copula Model of CDS Spreads, working paper, Duke University. 

[~, inx_a] = min( abs( lam_grid - theta(1) ) );
Gpdf_appro = Gpdf(:,inx_a);

xx = [x_grid, Gpdf_appro];

x = spline(xx(:,1),xx(:,2),q);

%%% switch to numerical integral if spline doesn't work properly
inx = find(x<0);
if isempty(inx) == 0
    disp('spline failed, switch to numerical integral');
    [Wa, Wb] = GLNodeWt(50);
    GLweight = [Wa Wb];
    for i = 1:length(inx)
        x(inx(i)) = GLquad('factor_cop_Gpdf_calc1_Skewtt_loading',1e-5,1-1e-5,GLweight,q(inx(i)),theta);
    end
    
end

out = x;
