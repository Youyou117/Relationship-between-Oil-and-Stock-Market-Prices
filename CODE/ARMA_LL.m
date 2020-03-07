function [LL] = ARMA_LL(theta,p,q,data);
%function [LL,LLa] = ARMA_LL(theta,p,q,data);
%
% Function to return the negative log-likelihood (total and individual) for
% an ARMA(p,q) with normal residuals
%
%  INPUTS:  theta, a kx1 vector = [phi0,phi1,...,phip,theta1,...,thetaq]
%           p, a scalar, the order of the AR part
%           q, a scalar, the order of the MA part
%           data, a Tx1 vector of data
%
%  OUTPUS:  LL, a scalar, the negative log-likelihood
%
%  Andrew Patton
%
%  31 March 2011

muhat = ARMA_mean(theta,p,q,data);
resids = data - muhat;

LLa = -1/2*log(2*pi) - 1/2*log(cov(resids)) - 1/2/cov(resids)*(resids.^2);
LLa = -LLa;
LL = sum(LLa);
