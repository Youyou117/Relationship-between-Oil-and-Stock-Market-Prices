function x = factor_cop_Gcdfinv_spline(q, theta, Gcdf, x_grid, lam_grid)
% function x = factor_cop_Gcdfinv_spline(q, theta, Gcdf, x_grid, lam_grid)
%
% The inverse cdf of x
%
% INPUTS:   q,        a vector in (0,1), the value of q that we will evaluate Ginv at
%           theta,    = [factor loadings; nuinv_z; nuinv_eps; psi_z], the parameters of Fz and Feps
%           Gcdf,     a [Num_x_grid by Num_lam_grid] matrix of the marginal cdfs of skew t-t factor model at x and factor loading (lam)
%           x_grid,   a vector of x that Gcdf is evaluated at
%           lam_grid, a vector of factor loading that Gcdf is evaluated at
%
%  OUTPUTS: x, a vector, the value of the inverse CDF at q
%
%  Dong Hwan Oh and Andrew Patton
%
%  16 April 2015
%
%  This code is to accompany the paper:
%  Oh, D.H. and A.J. Patton, 2015, Time-Varying Systemic Risk: Evidence
%  from a Dynamic Copula Model of CDS Spreads, working paper, Duke University.

N = length(q);

[~, inx_a] = min( abs( lam_grid - theta(1) ) );
Gcdf_appro = Gcdf(:,inx_a);

xx = [x_grid, Gcdf_appro];
[~, b, ~]= unique(xx(:,2),'first'); % just in case xx(:,2) has identical values
xx = xx(b,:);

aa_temp = (q>max(xx(:,2)));
bb_temp = (q<min(xx(:,2)));
q(aa_temp) = max(xx(:,2));
q(bb_temp) = min(xx(:,2));

max_xx = find(xx(:,2)>max(q), 1 ) +5;
min_xx = find(xx(:,2)<min(q), 1, 'last' ) -5;
inx_max = min([max_xx,size(xx,1)]);
inx_min = max([min_xx,1]);
xx = xx(inx_min:inx_max,:);

x4 = spline(xx(:,2),xx(:,1),q);

x2a = x4;

% adding a check that the spline quantiles are monotonic in q (linear interpolation guarantees this, but spline has possibility to do weird
% stuff). will "switch" any non-monotonic cases, in the spirit of Chernozhukov's paper on quantile estimation
temp = [(1:N)',q];
temp = sortrows(temp,2);  % sorting by q
x2b = sort(x2a);  % sorting the estimated quantiles
temp = sortrows([temp,x2b],1);  % putting things back in original order
sum(abs(x2b-x2a));  % just seeing whether this re-ordering changed anything
x(:,1) = temp(:,3);  % using the sorted quantiles as spline-based estimates

end
