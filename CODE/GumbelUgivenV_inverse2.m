function u = GumbelUgivenV_inverse2(v,t,k)
% function u = GumbelUgivenV_inverse2(v,t,k)
%
% Computes the ivnerse value of the conditional Gumbel (U given V)
% copula cdf at a specified point
% 
% INPUTS:	t, a scalar uniform r.v.
%				v, a scalar r.v. G(Y[t])
%				k, a scalar of kappa
%
% Saturday, 28 July, 2001
%
% Andrew Patton


info=struct('maxit',{100},'tol',{10^-6},'pflag',0);
[u,junk] = bisect2('GumbelUgivenV_t',10^-4,1-10^-4,info,v,t,k);


% dealing with the extreme values that don't converge
if junk==-999.99 		% then the zero of the function is very near to 0
   [u,junk] = bisect2('GumbelUgivenV_t',10^-18,10^-4,info,v,t,k);
elseif junk==999.99		% then the zero of the funciton is very near to 1
   [u,junk] = bisect2('GumbelUgivenV_t',1-10^-4,1-10^-18,info,v,t,k);
end
