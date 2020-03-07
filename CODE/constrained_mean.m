function  out1 = constrained_mean(mu,muhat,omega)
% function  out1 = constrained_mean(mu,muhat,omega)
%
%  Function used to find the estimate of the mean, subject to the
%  constraint that it is weakly positive. Used in implementing Wolak's test
%
%  INPUTS:  mu, a kx1 vector, the parameter to be estimated
%           muhat, a kx1 vector, the unconstrained estimate of the means
%           omegahat, a kxk matrix, the covariance matrix of the unconstrained estimator
%
%  OUPUTS:  out1, a scalar, the objective function: (muhat-mu)'*inv(omega)*(muhat-mu)
%
%  Andrew Patton
%
%  4 May 2008

out1 = (muhat-mu)'*inv(omega)*(muhat-mu);