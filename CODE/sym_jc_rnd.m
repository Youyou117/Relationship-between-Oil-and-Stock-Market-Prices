function out1 = sym_jc_rnd(tauU,tauL,T,state)
% function out1 = sym_jc_rnd(tauU,tauL,T,state)
% This program generates observations from
% a Symmetrised Joe-Clayton (BB7) copula
%
%  INPUTS:	tauU, a scalar or a Tx1 vector, the upper tail dependence parameter
%					tauL, a scalar or a Tx1 vector, the lower tail dependence parameter
%					T,	a scalar, the number of obs in each sample
%				    state, an integer to use to seed the random number generator
%
%  OUTPUT:	out1, a Tx2 matrix of data.
%
% Wednesday, 22 May, 2002
%
%	Andrew Patton
%
% On my Pentium II 266MHz laptop: T=1000 => ~20 seconds

% Published in:
%
% Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence,
% International Economic Review, 47(2), 527-556. 
%
% See also:
% Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of
% Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
%
% Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and
% Asymmetric Dependence for Asset Allocation, Journal of Financial
% Econometrics, 2(1), 130-168. 
%
% http://fmg.lse.ac.uk/~patton

% 1aug07: Important bug fixed following comment from Byoung Kang, Monash University, Australia.

if nargin<4
   rand('state',sum(1234*clock));	% setting RNG to new seed according to computer clock time.
else
   rand('state',state);
end
if nargin<5;
    direct=0;
end


if size(tauU,1)~=T;		% stretching the input vectors to match each other
   tauU = tauU*ones(T,1);
end
if size(tauL,1)~=T;
   tauL = tauL*ones(T,1);
end

% drawing from the joe-clayton copula 
temp1 = tau2kappa([tauU,tauL]);
temp2 = tau2kappa([tauL,tauU]);
temp3 = (rand(T,1)<=0.5);	% randomly picking BB7_1 or BB7_2
out1a = bb7_rnd(temp1(find(temp3==0),1),temp1(find(temp3==0),2),sum(temp3==0));
out1b = 1-bb7_rnd(temp2(find(temp3==1),1),temp2(find(temp3==1),2),sum(temp3==1));
out1  = [rand(T,1),[out1a;out1b]];
out1 = sortrows(out1,1);
out1 = out1(:,2:3);
